#include "bfv_ptct.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

using namespace lbcrypto;
using namespace std;

namespace {

struct CandidateResult {
    int logN;
    uint32_t ringDim;
    uint64_t ptm;
    int setDepth;
    uint32_t depthEstimate;
    double logQ;
    double log2Vnorm;
    double log2CMin;
    double log2CMaxWrap;
    double log2CMaxCap;
    double log2CMaxFeasible;
    double cWindowBits;
    double log2CSuggest;
    string rating;
};

double GetNormExpE1Inf(int logN) {
    switch (logN) {
        case 12:
            return 16.1807;
        case 13:
            return 17.3890;
        case 14:
            return 18.1876;
        case 15:
            return 19.2704;
        case 16:
            return 20.2343;
        default:
            throw runtime_error("Unsupported logN. Expected 12 <= logN <= 16.");
    }
}

double GetVnorm(CryptoContext<DCRTPoly> cc) {
    const auto cryptoParams = dynamic_pointer_cast<CryptoParametersBFVRNS>(cc->GetCryptoParameters());
    if (!cryptoParams) {
        throw runtime_error("Invalid crypto parameters: Expected BFVRNS.");
    }

    uint32_t n = cc->GetRingDimension();
    double sigma = cryptoParams->GetDistributionParameter();
    double alpha = cryptoParams->GetAssuranceMeasure();
    double Berr = sigma * sqrt(alpha);

    double Bkey = 0.0;
    uint32_t thresholdParties = cryptoParams->GetThresholdNumOfParties();
    if (cryptoParams->GetSecretKeyDist() == GAUSSIAN) {
        Bkey = sqrt(static_cast<double>(thresholdParties)) * Berr;
    }
    else {
        Bkey = static_cast<double>(thresholdParties);
    }

    auto delta = [](uint32_t nVal) -> double {
        return 4.0 * sqrt(static_cast<double>(nVal));
    };

    auto deltaMS = [](uint32_t nVal) -> double {
        return 4.0 * sqrt(static_cast<double>(nVal));
    };

    EncryptionTechnique encTech = cryptoParams->GetEncryptionTechnique();
    if (encTech == EXTENDED) {
        return (1.0 + deltaMS(n) * Bkey) / 2.0;
    }
    return Berr * (1.0 + 2.0 * delta(n) * Bkey);
}

bool ParseIntArg(const vector<string>& args, const string& key, int& out) {
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == key && i + 1 < args.size()) {
            out = stoi(args[i + 1]);
            return true;
        }
    }
    return false;
}

bool ParseDoubleArg(const vector<string>& args, const string& key, double& out) {
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == key && i + 1 < args.size()) {
            out = stod(args[i + 1]);
            return true;
        }
    }
    return false;
}

bool HasFlag(const vector<string>& args, const string& key) {
    return find(args.begin(), args.end(), key) != args.end();
}

double GetLogQ(CryptoContext<DCRTPoly> cc) {
    const auto cryptoParams = dynamic_pointer_cast<CryptoParametersBFVRNS>(cc->GetCryptoParameters());
    if (!cryptoParams) {
        throw runtime_error("Invalid crypto parameters: Expected BFVRNS.");
    }

    double logQ = 0.0;
    for (const auto& param : cryptoParams->GetElementParams()->GetParams()) {
        logQ += log2(param->GetModulus().ConvertToDouble());
    }
    return logQ;
}

CandidateResult EvaluateCandidate(int logN, uint64_t ptm, int setDepth, double lowerMarginBits,
                                  double upperMarginBits) {
    CCParams<CryptoContextBFVRNS> params;
    params.SetRingDim(1U << logN);
    params.SetPlaintextModulus(ptm);
    params.SetMultiplicativeDepth(setDepth);
    params.SetMultiplicationTechnique(HPSPOVERQLEVELED);
    params.SetSecurityLevel(SecurityLevel::HEStd_NotSet);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(params);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    CandidateResult out{};
    out.logN = logN;
    out.ringDim = cc->GetRingDimension();
    out.ptm = ptm;
    out.setDepth = setDepth;
    out.logQ = GetLogQ(cc);

    const double vNorm = GetVnorm(cc);
    out.log2Vnorm = log2(vNorm);

    const double log2t = log2(static_cast<double>(ptm));
    const double log2N = log2(static_cast<double>(out.ringDim));
    const double delta = 2.0 * sqrt(static_cast<double>(out.ringDim));
    const double denominator = log2t - 1.0 + log2(delta);
    const double numerator = out.logQ - 1.0 - log2t - out.log2Vnorm;

    if (denominator <= 0.0 || numerator <= 0.0) {
        out.depthEstimate = 0;
        out.log2CMin = numeric_limits<double>::infinity();
        out.log2CMaxWrap = -numeric_limits<double>::infinity();
        out.log2CMaxCap = log2t - 1.0;
        out.log2CMaxFeasible = -numeric_limits<double>::infinity();
        out.cWindowBits = -numeric_limits<double>::infinity();
        out.log2CSuggest = numeric_limits<double>::infinity();
        out.rating = "discard";
        return out;
    }

    out.depthEstimate = static_cast<uint32_t>(floor(numerator / denominator));
    out.log2CMaxCap = log2t - 1.0;

    if (out.depthEstimate == 0) {
        out.log2CMin = numeric_limits<double>::infinity();
        out.log2CMaxWrap = -numeric_limits<double>::infinity();
        out.log2CMaxFeasible = -numeric_limits<double>::infinity();
        out.cWindowBits = -numeric_limits<double>::infinity();
        out.log2CSuggest = numeric_limits<double>::infinity();
        out.rating = "discard";
        return out;
    }

    const double e1Bits = GetNormExpE1Inf(logN);
    const double geomBits = log2N - log2(static_cast<double>(M_PI));
    const double k = static_cast<double>(out.depthEstimate);

    // Model:
    // log2||e_K|| ~= e1Bits + K*log2(c) + (K-1)*(log2(N)-log2(pi))
    // We require Q/(2t) * 2^{lowerMargin} <= ||e_K|| <= Q/2 * 2^{-upperMargin}.
    // This creates a feasible window [log2CMin, log2CMaxFeasible].
    const double log2LowerTarget = out.logQ - 1.0 - log2t + lowerMarginBits;
    const double log2UpperTarget = out.logQ - 1.0 - upperMarginBits;
    out.log2CMin = (log2LowerTarget - e1Bits - (k - 1.0) * geomBits) / k;
    out.log2CMaxWrap = (log2UpperTarget - e1Bits - (k - 1.0) * geomBits) / k;
    out.log2CMaxFeasible = min(out.log2CMaxWrap, out.log2CMaxCap);
    out.cWindowBits = out.log2CMaxFeasible - out.log2CMin;
    out.log2CSuggest = (out.log2CMin + out.log2CMaxFeasible) / 2.0;

    if (out.cWindowBits >= 2.0) {
        out.rating = "recommended";
    }
    else if (out.cWindowBits >= 1.0) {
        out.rating = "viable";
    }
    else if (out.cWindowBits > 0.0) {
        out.rating = "stretch";
    }
    else {
        out.rating = "poor";
    }

    return out;
}

string BuildCommand(const CandidateResult& result, double customFactor, double retryFactor) {
    ostringstream oss;
    oss << "./build/bin/examples/pke/bfv_attack_recovery_gen"
        << " --logN " << result.logN
        << " --ptm " << result.ptm
        << " --depth " << result.setDepth
        << " --cF " << fixed << setprecision(2) << customFactor
        << " --retry " << fixed << setprecision(2) << retryFactor;
    return oss.str();
}

bool BetterCandidate(const CandidateResult& lhs, const CandidateResult& rhs) {
    if (lhs.setDepth != rhs.setDepth) {
        return lhs.setDepth < rhs.setDepth;
    }
    if (lhs.ptm != rhs.ptm) {
        return lhs.ptm < rhs.ptm;
    }
    if (lhs.cWindowBits != rhs.cWindowBits) {
        return lhs.cWindowBits > rhs.cWindowBits;
    }
    if (lhs.log2CMin != rhs.log2CMin) {
        return lhs.log2CMin < rhs.log2CMin;
    }
    return lhs.depthEstimate > rhs.depthEstimate;
}

void PrintUsage(const char* exeName) {
    cout << "Usage: " << exeName << " [--logN N] [--all] [--minDepth D] [--maxDepth D]"
         << " [--top K] [--lower-margin bits] [--upper-margin bits] [--retry r] [--cF x]" << endl;
}

}  // namespace

int main(int argc, char* argv[]) {
    vector<string> args(argv + 1, argv + argc);

    const vector<int> defaultLogNs = {15, 16};
    const vector<uint64_t> candidatePtms = {
        65537,     262145,    786433,    1048577,   2097153,   4194305,
        8388609,   16777217,  33423361,  67108865,  134086657, 270794753,
    };

    int onlyLogN = 0;
    int minDepth = 4;
    int maxDepth = 10;
    int topK = 6;
    double lowerMarginBits = 0.0;
    double upperMarginBits = 1.0;
    double retryFactor = 1.02;
    double customFactor = 1.0;
    const bool includeAll = HasFlag(args, "--all");

    if (HasFlag(args, "--help")) {
        PrintUsage(argv[0]);
        return 0;
    }

    ParseIntArg(args, "--logN", onlyLogN);
    ParseIntArg(args, "--minDepth", minDepth);
    ParseIntArg(args, "--maxDepth", maxDepth);
    ParseIntArg(args, "--top", topK);
    ParseDoubleArg(args, "--lower-margin", lowerMarginBits);
    ParseDoubleArg(args, "--upper-margin", upperMarginBits);
    ParseDoubleArg(args, "--retry", retryFactor);
    ParseDoubleArg(args, "--cF", customFactor);

    if (minDepth <= 0 || maxDepth < minDepth || topK <= 0) {
        cerr << "Invalid depth or top-K configuration." << endl;
        return 1;
    }
    if (onlyLogN != 0 && (onlyLogN < 12 || onlyLogN > 16)) {
        cerr << "--logN must be between 12 and 16." << endl;
        return 1;
    }

    vector<int> targetLogNs;
    if (onlyLogN != 0) {
        targetLogNs.push_back(onlyLogN);
    }
    else if (includeAll) {
        targetLogNs = {12, 13, 14, 15, 16};
    }
    else {
        targetLogNs = defaultLogNs;
    }

    cout << "==== BFV Attack Parameter Prescreener ====" << endl;
    cout << "Heuristic goals:" << endl;
    cout << "  1. lower bound: ||e_K|| >= Q/(2t) * 2^{lower-margin}" << endl;
    cout << "  2. upper bound: ||e_K|| <= Q/2 * 2^{-upper-margin}" << endl;
    cout << "Search range: depths in [" << minDepth << ", " << maxDepth << "]" << endl;
    cout << "Lower margin bits: " << lowerMarginBits << endl;
    cout << "Upper margin bits: " << upperMarginBits << endl;
    cout << "Displayed retry      : " << retryFactor << endl;
    cout << "Displayed customFactor: " << customFactor << endl;
    cout << endl;

    for (int logN : targetLogNs) {
        vector<CandidateResult> allResults;
        for (int depth = minDepth; depth <= maxDepth; ++depth) {
            for (uint64_t ptm : candidatePtms) {
                allResults.push_back(EvaluateCandidate(logN, ptm, depth, lowerMarginBits, upperMarginBits));
            }
        }

        sort(allResults.begin(), allResults.end(), BetterCandidate);

        vector<CandidateResult> filtered;
        for (const auto& result : allResults) {
            if (result.rating == "recommended" || result.rating == "viable") {
                filtered.push_back(result);
            }
        }
        if (filtered.empty()) {
            for (const auto& result : allResults) {
                if (result.rating == "stretch") {
                    filtered.push_back(result);
                }
            }
        }
        if (filtered.empty()) {
            filtered = allResults;
        }

        cout << string(88, '=') << endl;
        cout << "logN=" << logN << " (N=" << (1U << logN) << ")  e1_inf(max)=" << GetNormExpE1Inf(logN) << endl;
        cout << string(88, '-') << endl;
           cout << left << setw(14) << "rating"
               << setw(13) << "ptm"
               << setw(9) << "depth"
               << setw(9) << "Kest"
             << setw(10) << "log2Q"
             << setw(11) << "log2cMin"
             << setw(11) << "log2cMax"
             << setw(11) << "window"
             << setw(11) << "log2cUse"
             << setw(10) << "log2Vn"
             << endl;
        cout << string(88, '-') << endl;

        const size_t count = min(static_cast<size_t>(topK), filtered.size());
        for (size_t i = 0; i < count; ++i) {
            const auto& result = filtered[i];
              cout << left << setw(14) << result.rating
                  << setw(13) << result.ptm
                  << setw(9) << result.setDepth
                  << setw(9) << result.depthEstimate
                 << setw(10) << fixed << setprecision(1) << result.logQ
                 << setw(11) << setprecision(2) << result.log2CMin
                 << setw(11) << setprecision(2) << result.log2CMaxFeasible
                 << setw(11) << setprecision(2) << result.cWindowBits
                 << setw(11) << setprecision(2) << result.log2CSuggest
                 << setw(10) << setprecision(2) << result.log2Vnorm
                 << endl;
            cout << BuildCommand(result, customFactor, retryFactor) << endl;
        }

        if (count == 0) {
            cout << "No candidates found in the requested range." << endl;
        }

        cout << endl;
    }

    cout << "Interpretation:" << endl;
    cout << "  - log2cMin is the minimum c to pass the Q/(2t) lower bound." << endl;
    cout << "  - log2cMax is the maximum c allowed by min(Q/2 wrap bound, t/2 cap)." << endl;
    cout << "  - window = log2cMax - log2cMin. If <= 0, that tuple is infeasible." << endl;
    cout << "  - Kest is the depth estimate derived from the generated logQ and Vnorm." << endl;

    return 0;
}