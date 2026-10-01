#include "bfv_ptct.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using namespace lbcrypto;
using namespace std;

typedef std::complex<double> ComplexLD;

namespace {

struct Preset {
    int id;
    uint32_t ringDim;
    uint64_t ptm;
    uint32_t depth;
    const char* name;
};

struct Stats {
    double minVal;
    double q50;
    double q90;
    double q99;
    double maxVal;
    double mean;
};

Stats ComputeStats(vector<double> vals) {
    Stats out{0, 0, 0, 0, 0, 0};
    if (vals.empty())
        return out;

    sort(vals.begin(), vals.end());
    auto atQuantile = [&vals](double q) {
        if (vals.size() == 1)
            return vals[0];
        const double pos = q * static_cast<double>(vals.size() - 1);
        const size_t lo = static_cast<size_t>(floor(pos));
        const size_t hi = static_cast<size_t>(ceil(pos));
        if (lo == hi)
            return vals[lo];
        const double w = pos - static_cast<double>(lo);
        return vals[lo] * (1.0 - w) + vals[hi] * w;
    };

    out.minVal = vals.front();
    out.q50 = atQuantile(0.50);
    out.q90 = atQuantile(0.90);
    out.q99 = atQuantile(0.99);
    out.maxVal = vals.back();
    out.mean = accumulate(vals.begin(), vals.end(), 0.0) / static_cast<double>(vals.size());
    return out;
}

double CenteredInfNorm(DCRTPoly poly) {
    poly.SetFormat(Format::COEFFICIENT);
    auto polyInt = poly.CRTInterpolate();
    auto Q = polyInt.GetModulus();
    auto Qhalf = Q >> 1;

    double mx = 0.0;
    const size_t n = polyInt.GetLength();
    for (size_t i = 0; i < n; ++i) {
        auto v = polyInt[i];
        const double centered = (v > Qhalf)
                                         ? -static_cast<double>((Q - v).ConvertToDouble())
                                         : static_cast<double>(v.ConvertToDouble());
        mx = std::max(mx, std::abs(centered));
    }
    return static_cast<double>(mx);
}

DCRTPoly ExtractNoisePoly(PrivateKey<DCRTPoly> sk, ConstCiphertext<DCRTPoly> ct) {
    const vector<DCRTPoly>& cv = ct->GetElements();
    DCRTPoly s = sk->GetPrivateElement();

    size_t sizeQl = cv[0].GetParams()->GetParams().size();
    size_t sizeQs = s.GetParams()->GetParams().size();
    if (sizeQs > sizeQl) {
        s.DropLastElements(sizeQs - sizeQl);
    }

    DCRTPoly sPower(s);
    DCRTPoly e = cv[0];
    e.SetFormat(Format::EVALUATION);

    for (size_t i = 1; i < cv.size(); ++i) {
        DCRTPoly ci = cv[i];
        ci.SetFormat(Format::EVALUATION);
        e += (sPower * ci);
        if (i < cv.size() - 1) {
            sPower *= s;
        }
    }

    return e;
}

ComplexLD EvalComplexValue(DCRTPoly poly, size_t kIdx) {
    poly.SetFormat(Format::COEFFICIENT);
    auto intPoly = poly.CRTInterpolate();
    auto Q = intPoly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = intPoly.GetLength();

    double theta = (2.0 * kIdx - 1.0) * M_PI / static_cast<double>(N);
    ComplexLD omega(std::cos(theta), std::sin(theta));
    ComplexLD acc(0.0, 0.0);

    for (int j = static_cast<int>(N) - 1; j >= 0; --j) {
        auto v = intPoly[j];
        double c = (v > Qhalf)
                            ? -static_cast<double>((Q - v).ConvertToDouble())
                            : static_cast<double>(v.ConvertToDouble());
        acc = (acc * omega) + c;
    }

    return acc;
}

double GetVnorm(CryptoContext<DCRTPoly> cc) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersBFVRNS>(cc->GetCryptoParameters());
    if (!cryptoParams) {
        throw std::runtime_error("Invalid crypto parameters: Expected BFVRNS.");
    }

    uint32_t n = cc->GetRingDimension();
    double sigma = cryptoParams->GetDistributionParameter();
    double alpha = cryptoParams->GetAssuranceMeasure();
    double Berr = sigma * std::sqrt(alpha);

    double Bkey = 0;
    uint32_t thresholdParties = cryptoParams->GetThresholdNumOfParties();
    if (cryptoParams->GetSecretKeyDist() == GAUSSIAN) {
        Bkey = std::sqrt(static_cast<double>(thresholdParties)) * Berr;
    }
    else {
        Bkey = static_cast<double>(thresholdParties);
    }

    auto delta = [](uint32_t nVal) -> double {
        return 4.0 * std::sqrt(static_cast<double>(nVal));
    };

    auto deltaMS = [](uint32_t nVal) -> double {
        return 4.0 * std::sqrt(static_cast<double>(nVal));
    };

    EncryptionTechnique encTech = cryptoParams->GetEncryptionTechnique();
    if (encTech == EXTENDED) {
        return (1.0 + deltaMS(n) * Bkey) / 2.0;
    }
    return Berr * (1.0 + 2.0 * delta(n) * Bkey);
}

bool GetIntArg(const vector<string>& args, const string& key, int& out) {
    for (const auto& a : args) {
        if (a.rfind(key + "=", 0) == 0) {
            out = stoi(a.substr(key.size() + 1));
            return true;
        }
    }
    return false;
}

vector<Preset> BuildPresets() {
    return {
        {12, 4096, 65537, 2, "N=2^12"},
        {13, 8192, 65537, 2, "N=2^13"},
        {14, 16384, 65537, 2, "N=2^14"},
        {15, 32768, 65537, 2, "N=2^15"},
        {16, 65536, 65537, 2, "N=2^16"},
    };
}

CryptoContext<DCRTPoly> BuildContext(const Preset& p) {
    CCParams<CryptoContextBFVRNS> params;
    params.SetRingDim(p.ringDim);
    params.SetMultiplicationTechnique(HPSPOVERQLEVELED);
    params.SetSecurityLevel(SecurityLevel::HEStd_NotSet);
    params.SetPlaintextModulus(p.ptm);
    params.SetMultiplicativeDepth(p.depth);
    return GenCryptoContext(params);
}

void PrintStatsLine(const string& name, const Stats& s) {
    cout << name << " | min=" << setw(9) << s.minVal << "  p50=" << setw(9) << s.q50 << "  p90=" << setw(9)
         << s.q90 << "  p99=" << setw(9) << s.q99 << "  max=" << setw(9) << s.maxVal << "  mean=" << setw(9)
         << s.mean << '\n';
}

}  // namespace

int main(int argc, char* argv[]) {
    vector<string> args(argv + 1, argv + argc);

    int numKeys = 100;
    int ctPerKey = 1000;
    int kIdx = 1; // why not try all?
    int onlyParam = -1;

    GetIntArg(args, "--keys", numKeys);
    GetIntArg(args, "--cts", ctPerKey);
    GetIntArg(args, "--idx", kIdx);
    GetIntArg(args, "--only", onlyParam);

    if (numKeys <= 0 || ctPerKey <= 0) {
        cerr << "--keys and --cts must be > 0" << endl;
        return 1;
    }

    const auto presets = BuildPresets();

    cout << "==== BFV ||e0||_inf and ||e0*pt||_inf Profiling ====" << endl;
    cout << "Args: keys=" << numKeys << ", cts-per-key=" << ctPerKey << ", idx=" << kIdx;
    if (onlyParam != -1)
        cout << ", only=" << onlyParam;
    cout << "\n" << string(72, '=') << endl;

    for (const auto& p : presets) {
        if (onlyParam != -1 && p.id != onlyParam) {
            continue;
        }

        auto runStart = chrono::high_resolution_clock::now();

        CryptoContext<DCRTPoly> cc = BuildContext(p);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);

        const auto cryptoParams = std::static_pointer_cast<CryptoParametersBFVRNS>(cc->GetCryptoParameters());
        double logQ = 0.0;
        for (const auto& prm : cryptoParams->GetElementParams()->GetParams()) {
            logQ += log2(prm->GetModulus().ConvertToDouble());
        }

        uint32_t N = cc->GetRingDimension();
        if (kIdx < 1 || static_cast<uint32_t>(kIdx) > N) {
            cerr << "idx must satisfy 1 <= idx <= N for " << p.name << endl;
            return 1;
        }

        const double vNorm = GetVnorm(cc);
        const double logVnorm = log2(vNorm);

        Plaintext pt = LargePtxtApp(cc, static_cast<size_t>(kIdx), 1);
        DCRTPoly ptPoly = pt->GetElement<DCRTPoly>();
        ptPoly.SetFormat(Format::EVALUATION);
        const ComplexLD ptW = EvalComplexValue(pt->GetElement<DCRTPoly>(), static_cast<size_t>(kIdx));
        const double logPtW = log2(std::abs(ptW));

        vector<int64_t> zeros(N, 0);

        vector<double> log2E0;
        vector<double> log2E0Pt;
        vector<double> e0VsVnormBits;
        vector<double> e0PtGrowthBits;
        log2E0.reserve(static_cast<size_t>(numKeys) * static_cast<size_t>(ctPerKey));
        log2E0Pt.reserve(static_cast<size_t>(numKeys) * static_cast<size_t>(ctPerKey));
        e0VsVnormBits.reserve(static_cast<size_t>(numKeys) * static_cast<size_t>(ctPerKey));
        e0PtGrowthBits.reserve(static_cast<size_t>(numKeys) * static_cast<size_t>(ctPerKey));

        const size_t totalSamplesTarget = static_cast<size_t>(numKeys) * static_cast<size_t>(ctPerKey);
        size_t sampleCount = 0;
        int nextProgressPct = 10;

        cout << "\n[" << p.name << "]";
        for (int keyIdx = 0; keyIdx < numKeys; ++keyIdx) {
            KeyPair<DCRTPoly> kp = cc->KeyGen();

            for (int c = 0; c < ctPerKey; ++c) {
                Plaintext zeroPt = cc->MakeCoefPackedPlaintext(zeros);
                auto ct0 = cc->Encrypt(kp.publicKey, zeroPt);

                DCRTPoly e0 = ExtractNoisePoly(kp.secretKey, ct0);
                const double e0Inf = CenteredInfNorm(e0);

                DCRTPoly e0Mul = e0;
                e0Mul.SetFormat(Format::EVALUATION);
                e0Mul *= ptPoly;
                const double e0PtInf = CenteredInfNorm(e0Mul);

                const double lE0 = log2(max(e0Inf, std::numeric_limits<double>::min()));
                const double lE0Pt = log2(max(e0PtInf, std::numeric_limits<double>::min()));

                log2E0.push_back(lE0);
                log2E0Pt.push_back(lE0Pt);
                e0VsVnormBits.push_back(lE0 - logVnorm);
                e0PtGrowthBits.push_back(lE0Pt - lE0);
                sampleCount++;

                while (nextProgressPct <= 100 &&
                       sampleCount * static_cast<size_t>(100) >=
                           totalSamplesTarget * static_cast<size_t>(nextProgressPct)) {
                    cout << "..." << nextProgressPct << "%";
                    nextProgressPct += 10;
                }
            }
        }
        cout << endl;

        Stats sLogE0 = ComputeStats(log2E0);
        Stats sLogE0Pt = ComputeStats(log2E0Pt);
        Stats sE0VsVnorm = ComputeStats(e0VsVnormBits);
        Stats sGrowth = ComputeStats(e0PtGrowthBits);

        auto runEnd = chrono::high_resolution_clock::now();
        const double sec = chrono::duration<double>(runEnd - runStart).count();

        // cout << "\n" << string(72, '-') << endl;
        cout << "Preset: " << p.name << " (paramChoice=" << p.id << ")" << endl;
        cout << "N=" << N << ", log2Q=" << fixed << setprecision(2) << logQ << ", t=" << p.ptm
             << ", depth=" << p.depth << endl;
        cout << "Vnorm: log2(Vnorm)=" << logVnorm << endl;
        cout << "LargePtxtApp(idx=" << kIdx << ", constant=1): log2|pt(w)|=" << logPtW << endl;
        cout << "Samples: " << sampleCount << " (= keys * cts)" << endl;

        cout << setprecision(4);
        PrintStatsLine("log2 ||e0||_inf          ", sLogE0);
        PrintStatsLine("log2 ||e0*pt||_inf       ", sLogE0Pt);
        PrintStatsLine("log2(||e0||/Vnorm) bits  ", sE0VsVnorm);
        PrintStatsLine("log2(||e0*pt||/||e0||)   ", sGrowth);

        cout << "Elapsed: " << sec << " s" << endl;
        cout << string(72, '-') << endl;
    }

    cout << "Done." << endl;
    return 0;
}
