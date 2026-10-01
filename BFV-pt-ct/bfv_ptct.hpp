#include "openfhe.h"
#include "utils/inttypes.h"
#include <vector>
#include <complex>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <random>
#define _USE_MATH_DEFINES
// #include <bits/stdc++.h>

using namespace lbcrypto;
using namespace std;

#define DEBUG_EXTRACT_A0 1

static inline BigInteger ToBig(const NativeInteger& x) {
    return BigInteger(x.ConvertToInt());
}

vector<int> bitrev_indices(int n) {

    int b = 0;
    while ((1 << b) < n) ++b;

    vector<int> res(n);
    for (int t = 0; t < n; ++t) {
        unsigned x = static_cast<unsigned>(t);
        unsigned y = 0;
        for (int k = 0; k < b; ++k) {
            y = (y << 1) | ((x >> k) & 1u);
        }
        res[t] = static_cast<int>(y);
    }
    return res;
}

vector<complex<double>> apply_T(const vector<complex<double>>& x) {
    size_t n = x.size();
    if (n == 1)
        return x;
    if (n & (n - 1))
        throw invalid_argument("n must be a power of two");

    size_t h = n / 2;
    vector<complex<double>> x0(x.begin(), x.begin() + h);
    vector<complex<double>> x1(x.begin() + h, x.end());
    vector<complex<double>> y0 = apply_T(x0);
    vector<complex<double>> y1 = apply_T(x1);
    vector<complex<double>> out(n);

    uint64_t mod = static_cast<uint64_t>(4 * n);
    uint64_t expv = 1;

    for (size_t k = 0; k < h; ++k) {
        double angle = 2.0 * M_PI * static_cast<double>(expv) / static_cast<double>(mod);
        complex<double> w = polar(1.0, angle);
        complex<double> t = w * y1[k];
        out[k]     = y0[k] + t;
        out[k + h] = y0[k] - t;
        expv = (expv * 5ULL) % mod;
    }
    return out;
}

vector<complex<double>> apply_T_herm(const vector<complex<double>>& x) {
    size_t n = x.size();
    if (n == 1)
        return x;
    if (n & (n - 1))
        throw invalid_argument("n must be a power of two");

    size_t h = n / 2;
    vector<complex<double>> x0(x.begin(), x.begin() + h);
    vector<complex<double>> x1(x.begin() + h, x.end());
    vector<complex<double>> y0(h), y1(h);

    uint64_t mod = static_cast<uint64_t>(4 * n);
    uint64_t expv = 1;

    for (size_t k = 0; k < h; ++k) {
        double angle = -2.0 * M_PI * static_cast<double>(expv) / static_cast<double>(mod);
        complex<double> w = polar(1.0, angle);
        y0[k] = x0[k] + x1[k];
        y1[k] = w * (x0[k] - x1[k]);
        expv = (expv * 5ULL) % mod;
    }

    vector<complex<double>> out0 = apply_T_herm(y0);
    vector<complex<double>> out1 = apply_T_herm(y1);
    vector<complex<double>> out;
    out.reserve(n);
    out.insert(out.end(), out0.begin(), out0.end());
    out.insert(out.end(), out1.begin(), out1.end());
    return out;
}


vector<complex<double>> fft(const vector<double>& m) {
    size_t N = m.size();
    if (N % 2 != 0)
        throw invalid_argument("N must be even");

    size_t n = N / 2;
    if (n & (n - 1))
        throw invalid_argument("n must be a power of two");

    vector<int> rev = bitrev_indices(static_cast<int>(n));
    vector<complex<double>> m_concatenate(n);

    for (size_t i = 0; i < n; ++i) 
        m_concatenate[i] = complex<double>(m[rev[i]], m[rev[i] + n]);

    vector<complex<double>> out = apply_T(m_concatenate);

    return out;
}




vector<complex<double>> normalizedFFT(PrivateKey<DCRTPoly> privateKey, NativePoly& a){

    const auto cryptoParams = static_pointer_cast<CryptoParametersBFVRNS>(privateKey->GetCryptoParameters());
    CryptoContext<DCRTPoly> cryptoContext = privateKey->GetCryptoContext();

    auto paramsQ = cryptoParams->GetElementParams();
    usint N = paramsQ->GetRingDimension();
    // usint n = N / 2;
    BigInteger modulusQ = cryptoParams->GetElementParams()->GetModulus();

    vector<double> coeff(N, 0.0);

    NativeInteger qi = a.GetModulus();
    NativePoly a_tmp = a;
    a_tmp.SetFormat(COEFFICIENT);

    for (size_t k = 0; k < N; ++k) 
        coeff[k] = static_cast<double>(a_tmp.GetValues()[k].ConvertToInt()) / static_cast<double>(qi.ConvertToInt());


    for (auto& x : coeff)
        x -= round(x);

    vector<complex<double>> fft_result = fft(coeff);
    return fft_result;

}

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif
typedef std::complex<double> ComplexLD;

double evalNorm(DCRTPoly& a, size_t k_idx) {
    a.SetFormat(Format::COEFFICIENT);
    auto poly = a.CRTInterpolate(); 
    auto Q = poly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = poly.GetLength();

    double theta = (2.0 * k_idx - 1.0) * M_PI / static_cast<double>(N);
    ComplexLD omega_k(std::cos(theta), std::sin(theta));

    ComplexLD acc(0.0, 0.0);
    ComplexLD omega_pow(1.0, 0.0);
    
    for (size_t j = 0; j < N; ++j) {
        auto v = poly[j];
        double c_j;

        if (v > Qhalf) {
            auto neg_v = Q - v;
            c_j = - static_cast<double>(neg_v.ConvertToDouble());
        } else {
            c_j = static_cast<double>(v.ConvertToDouble());
        }

        acc += c_j * omega_pow;
        omega_pow *= omega_k;
    }

    return static_cast<double>(std::abs(acc));
}

//-- for Debugging --//
double evalNorm(DCRTPoly& a, size_t k_idx, std::string debug_info) {
    a.SetFormat(Format::COEFFICIENT);
    auto poly = a.CRTInterpolate(); 
    auto Q = poly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = poly.GetLength();

    double theta = (2.0 * k_idx - 1.0) * M_PI / static_cast<double>(N);
    ComplexLD omega_k(std::cos(theta), std::sin(theta));

    ComplexLD acc(0.0, 0.0);
    ComplexLD omega_pow(1.0, 0.0);

    // //-- Debug --//
    // std::cout << "\nlog2 Q: " << std::log2(Q.ConvertToDouble()) << std::endl;
    // std::cout << "[[[" << debug_info << "]]]\nCoefficient: ";
    // double count = 0;
    // //-----------//


    for (size_t j = 0; j < N; ++j) {
        auto v = poly[j];
        double c_j;

        //-- Debug --//
        if (v > Q) {
            std::cerr << "Error: Coefficient exceeds modulus." << std::endl;
            break;
        }
        //-----------//
        if (v > Qhalf) {
            auto neg_v = Q - v;
            c_j = - static_cast<double>(neg_v.ConvertToDouble());
        } else {
            c_j = static_cast<double>(v.ConvertToDouble());
        }

        // //-- Debug --//
        // if (std::abs(c_j) < 1000000) {
        //     count += 1;
        // }
        // if (j < 20) {
        //     std::cout << c_j << ", ";
        // }
        // //-----------//
        
        // 2. Horner's method
        acc += c_j * omega_pow;
        omega_pow *= omega_k; // omega^j for Honor's method. 
    }
    
    // //-- Debug --//
    // std::cout << "..." << std::endl;
    // // std::cout << "Max coefficient: " << std::log2(max_c_j) << std::endl;
    // std::cout << "Coefficient close to Qhalf counted: " << count << std::endl;
    // //-----------//

    return static_cast<double>(std::abs(acc));
}
//-----------//


vector<complex<double>> normalizedFFT(PrivateKey<DCRTPoly> privateKey, DCRTPoly& a){
    const auto cryptoParams = static_pointer_cast<CryptoParametersBFVRNS>(privateKey->GetCryptoParameters());
    CryptoContext<DCRTPoly> cryptoContext = privateKey->GetCryptoContext();

    size_t sizeQ = a.GetParams()->GetParams().size();

    auto paramsQ = cryptoParams->GetElementParams();
    usint N = paramsQ->GetRingDimension();
    // usint n = N / 2;
    BigInteger modulusQ = cryptoParams->GetElementParams()->GetModulus();

    vector<double> coeff(N, 0.0);

    for (size_t i = 0; i < sizeQ; ++i) {

        const NativePoly& ai = a.GetElementAtIndex(i);
        NativeInteger qi = ai.GetModulus();

        BigInteger QHati = modulusQ / BigInteger(qi);
        NativeInteger QHati_inv = QHati.ModInverse(qi).ConvertToInt();
        
        NativePoly a_tmp = ai;
        a_tmp *= QHati_inv;
        a_tmp.SetFormat(COEFFICIENT);

        for (size_t k = 0; k < N; ++k) 
            coeff[k] += static_cast<double>(a_tmp.GetValues()[k].ConvertToInt()) / static_cast<double>(qi.ConvertToInt());

    }

    for (auto& x : coeff)
        x -= round(x);

    vector<complex<double>> fft_result = fft(coeff);
    return fft_result;

}

void EvalNoiseBFV(PrivateKey<DCRTPoly> privateKey, ConstCiphertext<DCRTPoly> ciphertext, usint ptm,
                double& noise, double& logQ) {
    const auto cryptoParams = static_pointer_cast<CryptoParametersBFVRNS>(privateKey->GetCryptoParameters());

    const vector<DCRTPoly>& cv = ciphertext->GetElements();
    DCRTPoly s                      = privateKey->GetPrivateElement();
    

    size_t sizeQl = cv[0].GetParams()->GetParams().size();
    size_t sizeQs = s.GetParams()->GetParams().size();

    size_t diffQl = sizeQs - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly sPower(scopy);

    DCRTPoly b = cv[0];
    b.SetFormat(Format::EVALUATION);

    DCRTPoly ci;
    for (size_t i = 1; i < cv.size(); i++) {
        ci = cv[i];
        ci.SetFormat(Format::EVALUATION);

        b += sPower * ci;
        sPower *= scopy;
    }

    DCRTPoly res;
    res = b;

    // Converts back to coefficient representation
    res.SetFormat(Format::COEFFICIENT);
    noise        = (log2(res.Norm()));

    // for (uint i=0;i<100;i++) {
    //     auto x = res.GetElementAtIndex(0).GetValues()[i];
    //     if (2*x < res.GetModulus().ConvertToInt())
    //         cout << x << " ";
    //     else
    //         cout << "-" << res.GetModulus().ConvertToInt()-x << " ";
    // }
    // cout << noise << endl;

    // size_t sizeQ = cryptoParams->GetElementParams()->GetParams().size();
    // logQ = 0;
    // for (usint i = 0; i < sizeQ; i++) {
    //     double logqi = log2(cryptoParams->GetElementParams()->GetParams()[i]->GetModulus().ConvertToInt());
    //     logQ += logqi;
    // }

    // cout << "logQ: " << logQ << endl;
    // cout << "noise: " << noise << endl;

}

void EvalNoiseBFV(PrivateKey<DCRTPoly> privateKey, ConstCiphertext<DCRTPoly> ciphertext, 
                double& noise, double& noise_idx, size_t idx) {
    const auto cryptoParams = static_pointer_cast<CryptoParametersBFVRNS>(privateKey->GetCryptoParameters());

    const vector<DCRTPoly>& cv = ciphertext->GetElements();
    DCRTPoly s                      = privateKey->GetPrivateElement();
    

    size_t sizeQl = cv[0].GetParams()->GetParams().size();
    size_t sizeQs = s.GetParams()->GetParams().size();

    size_t diffQl = sizeQs - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly sPower(scopy);

    DCRTPoly b = cv[0];
    b.SetFormat(Format::EVALUATION);

    DCRTPoly ci;
    for (size_t i = 1; i < cv.size(); i++) {
        ci = cv[i];
        ci.SetFormat(Format::EVALUATION);

        b += sPower * ci;
        sPower *= scopy;
    }

    DCRTPoly res;
    res = b;
 
    // Converts back to coefficient representation
    res.SetFormat(Format::COEFFICIENT);
    noise        = (log2(res.Norm()));
    noise_idx    = (log2(evalNorm(res, idx, "EvalNoiseBFV")));

}








Plaintext LargestPtxt(CryptoContext<DCRTPoly> cryptoContext,
                      size_t n,
                      size_t idx,
                      uint64_t ptm) {

    // Number of coeffs is 2n (= ring dimension)
    const size_t N = 2 * n;
    if (N == 0) {
        OPENFHE_THROW("N must be > 0");
    }

    // k in equation (1 <= k <= N)
    // Using idx as 1-based index
    size_t k = idx;
    if (k < 1 || k > N) {
        OPENFHE_THROW("idx must satisfy 1 <= idx <= N");
    }

    // floor((t-1)/2)
    const int64_t scale = static_cast<int64_t>((ptm - 1) / 4);

    std::vector<int64_t> coeffs(N, 0);

    // pt_j = Sign(cos(k*pi*(j+1/2)/N)) * floor((t-1)/2),  j=0..N-1
    for (size_t j = 0; j < N; ++j) {
        const double angle =
            static_cast<double>(k) * M_PI *
            (static_cast<double>(j) + 0.5) /
            static_cast<double>(N);

        const double c = cos(angle);
        const int64_t sgn = (c >= 0.0) ? 1 : -1;  // Sign(0)=+1
        coeffs[j] = sgn * scale;
    }

    return cryptoContext->MakeCoefPackedPlaintext(coeffs);
}

Plaintext LargePtxtApp(CryptoContext<DCRTPoly> cryptoContext, size_t k_idx, int64_t constant = 1) {
    const uint32_t N = cryptoContext->GetRingDimension();
    
    if (k_idx < 1 || k_idx > N) {
        OPENFHE_THROW("k_idx must satisfy 1 <= k_idx <= N");
    }

    // 1. Coeffs Gen (Sign of Cosine)
    std::vector<int64_t> coeffs(N);
    for (uint32_t j = 0; j < N; ++j) {
        const double term_k = 2.0 * static_cast<double>(k_idx) - 1.0;
        const double angle = (term_k * M_PI * (static_cast<double>(j) + 0.5)) / 
                                   static_cast<double>(N);
        const double cos_val = std::cos(angle);
        coeffs[j] = (cos_val >= 0.0) ? 1 : 0;  // -1 was stored as 65536 = ptm -1 ... why?
        coeffs[j] *= constant;
    }

    // std::cout << "\n======= Print Coeffs =======" << std::endl;
    // for (uint32_t j = 0; j < 20; ++j) {
    //     std::cout << coeffs[j] << ", ";
    // }
    // std::cout << std::endl;

    // 2. Plaintext Gen
    Plaintext plaintext = cryptoContext->MakeCoefPackedPlaintext(coeffs);
    DCRTPoly poly = plaintext->GetElement<DCRTPoly>();

    // FOR DEBUGGING //
    // // 3. using evalNorm
    // double val_from_evalNorm = evalNorm(poly, k_idx);

    // // 4. theory
    // double sin_approx = 1.0 / std::sin(M_PI / (2.0 * (double)N));

    // // 5. compare
    // std::cout << "\n======= [LargePtxtApp comparison] =======" << std::endl;
    // std::cout << "- Target Index (k). : " << k_idx << std::endl;
    // std::cout << "- [theo] log2(1/sin): " << std::log2(sin_approx) << std::endl;
    // std::cout << "- [real] evalNorm   : " << std::log2(val_from_evalNorm) << std::endl;
    
    // // any offset?
    // double diff = std::log2(val_from_evalNorm) - std::log2(sin_approx);
    // std::cout << "- diff (bits): " << diff << std::endl;
    // std::cout << "=============================================\n" << std::endl;

    return plaintext;
}



//////////////////////
///// New method /////
//////////////////////


void invNativePoly(const NativePoly& a_in,
                        const shared_ptr<ILNativeParams>& param_qi,
                        NativeInteger qi,
                        NativePoly& a_inv_out) {

    NativePoly a = a_in;
    a.SetFormat(Format::EVALUATION);
    const auto& evals = a.GetValues();

    NativeVector inv_evals(evals.GetLength(), qi);
    for (size_t j = 0; j < evals.GetLength(); ++j)
        inv_evals[j] = evals[j].ModInverse(qi);

    a_inv_out = NativePoly(param_qi);
    a_inv_out.SetValues(std::move(inv_evals), Format::EVALUATION);
}

// DivideOutLastPrime(const NativeInteger& t, const std::vector<NativeInteger>& tModqPrecon,
//                                       const NativeInteger& negtInvModq, const NativeInteger& negtInvModqPrecon,
//                                       const std::vector<NativeInteger>& qlInvModq,
//                                       const std::vector<NativeInteger>& qlInvModqPrecon) {
//     DCRTPolyImpl::PolyType delta(m_vectors.back());
//     delta.SetFormat(Format::COEFFICIENT);
//     delta *= negtInvModq;
//     this->DropLastElement();
//     size_t size{m_vectors.size()};

// #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(size))
//     for (size_t i = 0; i < size; ++i) {
//         auto tmp{delta};
//         tmp.SwitchModulus(m_vectors[i].GetModulus(), m_vectors[i].GetRootOfUnity(), 0, 0);
//         if (m_format == Format::EVALUATION)
//             tmp.SwitchFormat();
//         m_vectors[i] += (tmp *= t);
//         m_vectors[i] *= qlInvModq[i];
//     }
// }


DCRTPoly DivideOutLastPrime(const DCRTPoly& aIn) {

    DCRTPoly a(aIn);
    a.SetFormat(Format::COEFFICIENT);

    const size_t Lp1 = a.GetNumOfElements(); // = L+1
    if (Lp1 == 0)
        return a; // Empty multi-tower

    // Last tower r (mod q_ell)
    NativePoly r = a.GetElementAtIndex(Lp1 - 1);
    r.SetFormat(Format::COEFFICIENT);
    const NativeInteger ql = r.GetModulus();

    // Cache coefficient vector of r
    const NativeVector& rVals = r.GetValues();

    // Drop last tower: Q' = Q / q_ell
    a.DropLastElement();
    const size_t L = a.GetNumOfElements();

    auto paramsQ = a.GetParams();

    for (size_t i = 0; i < L; ++i) {
        NativePoly ai = a.GetElementAtIndex(i);      // COEFF
        ai.SetFormat(Format::COEFFICIENT);

        // Tower params/modulus
        auto param_qi = paramsQ->GetParams()[i];
        const NativeInteger qi = ai.GetModulus();

        // r_i := r mod q_i (COEFF)
        NativeVector rModQi(rVals.GetLength(), qi);
        for (size_t j = 0; j < rVals.GetLength(); ++j)
            rModQi[j] = rVals[j].Mod(qi);

        NativePoly r_i(param_qi);
        r_i.SetValues(std::move(rModQi), Format::COEFFICIENT);

        // inv := q_ell^{-1} mod q_i
        const NativeInteger ql_inv = ql.ModInverse(qi);

        // a_i <- (a_i - r_i) * ql_inv mod q_i
        ai -= r_i;
        ai *= ql_inv;

        a.SetElementAtIndex(i, std::move(ai));
    }

    a.SetFormat(Format::COEFFICIENT);
    return a;
}



// NativePoly superLargeApart(size_t n, size_t idx, NativeInteger qi,
//                                     const shared_ptr<ILNativeParams>& param_qi) {
    
//     vector<complex<double>> z(n, {0.0, 0.0});
//     z[idx] = complex<double>(static_cast<double>(n) / 2.0, 0.0);

//     vector<double> u = ifft(z); // size = 2n; [Re(y) | Im(y)]

//     // Generate NativePoly (coefficients mod q0)
//     const uint32_t rdim = 2 * n;
//     NativeVector coeffs(rdim, qi);

//     const double qi_ld = static_cast<double>(qi.ConvertToInt());

//     for (size_t i = 0; i < rdim; ++i) {
//         double scaled = static_cast<double>(u[i]) * qi_ld;
//         int64_t rounded = static_cast<int64_t>(llround(scaled));
//         if (rounded < 0) coeffs[i] = static_cast<NativeInteger::Integer>(qi)-static_cast<NativeInteger::Integer>(-rounded);
//         else coeffs[i] = static_cast<NativeInteger::Integer>(rounded);
//     }

//     NativePoly a_qi(param_qi);
//     a_qi.SetValues(std::move(coeffs), Format::COEFFICIENT);
//     a_qi.SetFormat(COEFFICIENT);
    
//     return a_qi;
    
// }


// NativePoly Extract_a0(const DCRTPoly& aIn) {
//     DCRTPoly cur(aIn);
//     // Keep as COEFF to avoid unnecessary format conversion in DivideOutLastPrime
//     cur.SetFormat(Format::COEFFICIENT);

//     // Sequentially remove qL, q{L-1}, ..., q1
//     while (cur.GetNumOfElements() > 1) {
//         cur = DivideOutLastPrime(cur); // (cur - [cur]_{last}) / last
//         cur.SetFormat(Format::COEFFICIENT);
//     }

//     // Result is single tower (q0) — i.e., a_0 \in R_{q0}
//     NativePoly a0 = cur.GetElementAtIndex(0);
//     a0.SetFormat(Format::COEFFICIENT);
//     return a0;
// }


NativePoly Extract_a0(const DCRTPoly& aIn) {
    DCRTPoly cur(aIn);
    cur.SetFormat(Format::COEFFICIENT);

#if DEBUG_EXTRACT_A0
    // --- Debug preparation: Collect original constant terms (for CRT reconstruction) ---
    const size_t Lp1_init = cur.GetNumOfElements(); // L+1
    std::vector<NativeInteger> crtResid;  // [c_0, c_1, ..., c_L]
    std::vector<NativeInteger> crtMod;    // [q_0, q_1, ..., q_L]
    crtResid.reserve(Lp1_init);
    crtMod.reserve(Lp1_init);

    for (size_t i = 0; i < Lp1_init; ++i) {
        NativePoly ai0 = cur.GetElementAtIndex(i);
        ai0.SetFormat(Format::COEFFICIENT);
        crtResid.push_back(ai0.GetValues()[0]);    // Constant term
        crtMod.push_back(ai0.GetModulus());        // q_i
    }

    // Collect Garner digits (constant terms): [a_L, a_{L-1}, ..., a_0]
    std::vector<NativeInteger> digits; // Constants only
    std::vector<NativeInteger> moduli; // q_last for each step
    digits.reserve(Lp1_init);
    moduli.reserve(Lp1_init);
#endif

    // Collect digits while removing qL, q_{L-1}, ..., q1
    while (cur.GetNumOfElements() > 1) {
#if DEBUG_EXTRACT_A0
        size_t idxLast = cur.GetNumOfElements() - 1;
        NativePoly last = cur.GetElementAtIndex(idxLast);
        last.SetFormat(Format::COEFFICIENT);
        const NativeInteger qLast = last.GetModulus();
        const NativeInteger cLast = last.GetValues()[0]; // Constant term

        digits.push_back(cLast);  // constant term of a_last
        moduli.push_back(qLast);  // q_last
#endif
        cur = DivideOutLastPrime(cur); // (cur - [cur]_{last}) / last
        cur.SetFormat(Format::COEFFICIENT);
    }

    // Result is single tower (q0)
    NativePoly a0 = cur.GetElementAtIndex(0);
    a0.SetFormat(Format::COEFFICIENT);

#if DEBUG_EXTRACT_A0
    // push last a0 constant/modulus (digits = [a_L,...,a_0])
    {
        const NativeInteger c0 = a0.GetValues()[0];
        const NativeInteger q0 = a0.GetModulus();
        digits.push_back(c0);
        moduli.push_back(q0);
    }

    // Q = prod q_i
    BigInteger Q(1);
    for (const auto& qi : crtMod) Q *= ToBig(qi);

    // CRT Reconstruction (using original residues)
    BigInteger crt(0);
    for (size_t i = 0; i < crtResid.size(); ++i) {
        BigInteger qi = ToBig(crtMod[i]);
        BigInteger ci = ToBig(crtResid[i]);
        BigInteger Mi = Q / qi;
        BigInteger inv = Mi.ModInverse(qi);
        crt += (ci * Mi * inv) % Q;
        crt %= Q;
    }

    // Garner Base expansion reconstruction (digits = [a_L,...,a_0])
    BigInteger g(0), prod(1);
    for (size_t k = 0; k < digits.size(); ++k) {
        const BigInteger ak = ToBig(digits[k]);
        const BigInteger qk = ToBig(moduli[k]);
        g += ak * prod;    // a_k * ∏_{j>k} q_j
        g %= Q;
        prod *= qk;
        prod %= Q;
    }

    // Output
    std::cout << "[Extract_a0::debug] L+1 = " << Lp1_init << "\n";
    std::cout << "  q_i (i=0..L): ";
    for (size_t i = 0; i < crtMod.size(); ++i) {
        if (i) std::cout << ", ";
        std::cout << crtMod[i];
    }
    std::cout << "\n  c_i = const-term(a mod q_i): ";
    for (size_t i = 0; i < crtResid.size(); ++i) {
        if (i) std::cout << ", ";
        std::cout << crtResid[i];
    }
    std::cout << "\n  digits (a_L..a_0 const-terms): ";
    for (size_t i = 0; i < digits.size(); ++i) {
        if (i) std::cout << ", ";
        std::cout << digits[i];
    }
    std::cout << "\n  CRT(const) = " << crt
              << "\n  Garner(const) = " << g
              << "\n  Match? " << ((crt % Q) == (g % Q) ? "YES" : "NO")
              << "\n\n";
#endif

    return a0;
}



// Plaintext makeApartLarge(PrivateKey<DCRTPoly> privateKey, ConstCiphertext<DCRTPoly> ciphertext, size_t idx) {
//     const auto cryptoParams = static_pointer_cast<CryptoParametersBFVRNS>(privateKey->GetCryptoParameters());
//     CryptoContext<DCRTPoly> cryptoContext = privateKey->GetCryptoContext();

//     const vector<DCRTPoly>& cv = ciphertext->GetElements();
//     DCRTPoly a = cv[1];
//     size_t sizeQ = a.GetParams()->GetParams().size();

//     auto paramsQ = cryptoParams->GetElementParams();
//     usint N = paramsQ->GetRingDimension();
//     usint n = N / 2;
//     BigInteger modulusQ = cryptoParams->GetElementParams()->GetModulus();

//     vector<NativeInteger> moduliQ(sizeQ);
//     for (size_t i = 0; i < sizeQ; i++)
//         moduliQ[i] = a.GetElementAtIndex(i).GetModulus();

//     Plaintext ptxt;

//     NativePoly a0 = Extract_a0(a);
//     NativeInteger q0 = moduliQ[0];
//     auto param_q0 = paramsQ->GetParams()[0];
//     NativePoly a0_inv;

//     NativePoly aLarge_q0 = superLargeApart(n, idx, q0, param_q0);

//     vector<complex<double>> evals = normalizedFFT(privateKey, aLarge_q0); // Verification of normalization

//     cout << evals[idx].real() << " + " << evals[idx].imag() << "i"
//         << "  |a_tilde(w_" << idx << ")| = " << abs(evals[idx]) << endl;

//     invNativePoly(a0, param_q0, q0, a0_inv);

//     a0_inv.SetFormat(EVALUATION);

//     aLarge_q0.SetFormat(EVALUATION);
//     aLarge_q0 *= a0_inv;

//     aLarge_q0.SetFormat(COEFFICIENT);
    

//     NativeVector v = aLarge_q0.GetValues();

//     vector<int64_t> coeff(N);
//     const int64_t q0_int  = static_cast<int64_t>(q0.ConvertToInt());
//     const int64_t half_up = (q0_int + 1) / 2;  // ceil(q/2)
//     for (size_t k = 0; k < v.GetLength(); ++k) {
//         int64_t x = static_cast<int64_t>(v[k].ConvertToInt()); // in [0, q-1]
//         if (x >= half_up) x -= q0_int; // center-lift
//         coeff[k] = x;                  // [-q/2, q/2-1] (even q), [-(q-1)/2,(q-1)/2] (odd q)
//     }

//     ptxt = cryptoContext->MakeCoefPackedPlaintext(coeff);

//     return ptxt;
// }


// Plaintext makeApartLarge(PrivateKey<DCRTPoly> privateKey, ConstCiphertext<DCRTPoly> ciphertext, size_t idx) {
//     const auto cryptoParams = static_pointer_cast<CryptoParametersBFVRNS>(privateKey->GetCryptoParameters());
//     CryptoContext<DCRTPoly> cryptoContext = privateKey->GetCryptoContext();

//     const vector<DCRTPoly>& cv = ciphertext->GetElements();
//     DCRTPoly a = cv[1];
//     size_t sizeQ = a.GetParams()->GetParams().size();

//     auto paramsQ = cryptoParams->GetElementParams();
//     usint N = paramsQ->GetRingDimension();
//     usint n = N / 2;
//     BigInteger modulusQ = cryptoParams->GetElementParams()->GetModulus();

//     vector<NativeInteger> moduliQ(sizeQ);
//     for (size_t i = 0; i < sizeQ; i++)
//         moduliQ[i] = a.GetElementAtIndex(i).GetModulus();

//     bool found = false;

//     Plaintext ptxt;

//     for (size_t i = 0; i < sizeQ; ++i) {
//         const NativePoly& ai = a.GetElementAtIndex(i);
//         NativeInteger qi = moduliQ[i];
//         auto param_qi = paramsQ->GetParams()[i];
//         NativePoly ai_inv;

//         if (!invNativePoly(ai, param_qi, qi, ai_inv))
//             continue;

//         // std::cout << "Found inverse modulo q_" << i << std::endl;
//         found = true;
//         ai_inv.SetFormat(EVALUATION);

//         BigInteger QHati = modulusQ / BigInteger(qi);
//         NativeInteger Qhati = QHati.Mod(qi).ConvertToInt();
//         NativePoly aLarge_qi = superLargeApart(n, idx, qi, param_qi);

//         vector<complex<double>> evals = normalizedFFT(privateKey, aLarge_qi); // verification of normalization

//         cout << evals[idx].real() << " + " << evals[idx].imag() << "i"
//             << "  |a_tilde(w_" << idx << ")| = " << abs(evals[idx]) << endl;

//         aLarge_qi.SetFormat(EVALUATION);
//         aLarge_qi *= ai_inv;
//         aLarge_qi *= Qhati;

//         // NativeInteger Qhati_inv = Qhati.ModInverse(qi).ConvertToInt();
//         // aLarge_qi *= ai;
//         // aLarge_qi *= Qhati_inv;

//         aLarge_qi.SetFormat(COEFFICIENT);

//         // vector<complex<double>> evals1 = normalizedFFT(privateKey, aLarge_qi); // verification of normalization

//         // cout << evals1[idx].real() << " + " << evals1[idx].imag() << "i"
//         //     << "  |pt*A0(w_" << idx << ")| = " << abs(evals1[idx]) << endl;

//         NativeVector v = aLarge_qi.GetValues();

//         vector<int64_t> coeff(N);
//         const int64_t qi_int  = static_cast<int64_t>(qi.ConvertToInt());
//         const int64_t half_up = (qi_int + 1) / 2;  // ceil(q/2)
//         for (size_t k = 0; k < v.GetLength(); ++k) {
//             int64_t x = static_cast<int64_t>(v[k].ConvertToInt()); // in [0, q-1]
//             if (x >= half_up) x -= qi_int; // center-lift
//             coeff[k] = x;                  // [-q/2, q/2-1] (even q), [-(q-1)/2,(q-1)/2] (odd q)
//         }

//         ptxt = cryptoContext->MakeCoefPackedPlaintext(coeff);
//         break;
//     }

//     if (!found)
//         OPENFHE_THROW("No invertible component found");

//     return ptxt;
// }

Ciphertext<DCRTPoly> PCMult(Plaintext ptxt, ConstCiphertext<DCRTPoly> ciphertext, size_t idx) {
    
    const vector<DCRTPoly>& cv = ciphertext->GetElements();

    auto bnew(cv[0]);
    bnew.SetFormat(EVALUATION);
    auto anew(cv[1]);
    anew.SetFormat(EVALUATION);

    DCRTPoly plain = ptxt->GetElement<DCRTPoly>();

    plain.SetFormat(COEFFICIENT);
    double norm1 = evalNorm(plain, idx);

    cout << "log2|pt(w)| = " << log2(norm1) << endl;

    plain.SetFormat(EVALUATION);

    double norm2 = evalNorm(plain, idx);

    cout << "log2|pt(w)| = " << log2(norm2) << endl;

    bnew *= plain;
    anew *= plain;

    auto newCtxt = ciphertext->CloneEmpty();
    newCtxt->SetElements({std::move(bnew), std::move(anew)});
    newCtxt->SetNoiseScaleDeg(ciphertext->GetNoiseScaleDeg());
    return newCtxt;

}