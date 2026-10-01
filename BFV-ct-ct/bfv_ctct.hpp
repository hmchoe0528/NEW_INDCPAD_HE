#include "openfhe.h"
#include "math/dftransform.h"
#include "utils/inttypes.h"
#include <vector>
#include <complex>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
#define _USE_MATH_DEFINES
// #include <bits/stdc++.h>

using namespace lbcrypto;
using namespace std;

#define DEBUG_EXTRACT_A0 1
// <cmath> defines M_PI on the platforms we build for; the fallback is for
// those where it is not exposed under the strict standard modes.
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// ─── Zero-plaintext BFV encryption with full randomness trace ───
// Only for Strong KR-, not used. 
// Full public-key encryption of m = 0, returning the sampled randomness
// (u, e0, e1) alongside the ciphertext:
//
//     c0 = pk_b · u + ns · e0
//     c1 = pk_a · u + ns · e1
//     c0 + c1·s = ns · (-u·e_pk + e0 + e1·s)
//
// u ~ ternary (or DGG if the secret-key distribution is GAUSSIAN),
// e0, e1 ~ DGG. The trace stores u, e0, e1 in COEFFICIENT format for
// downstream evaluation at roots of unity.
//
// Only STANDARD encryption technique (OpenFHE default). Throws on EXTENDED.
struct EncryptTrace {
    Ciphertext<DCRTPoly> ciphertext;
    DCRTPoly             u;    // randomness, COEFFICIENT
    DCRTPoly             e0;   // DGG noise for c0 (a.k.a. e_b), COEFFICIENT
    DCRTPoly             e1;   // DGG noise for c1 (a.k.a. e_a), COEFFICIENT
};

inline EncryptTrace EncryptZeroWithTrace(const PublicKey<DCRTPoly>& publicKey) {
    const auto cp = std::dynamic_pointer_cast<CryptoParametersBFVRNS>(
        publicKey->GetCryptoParameters());
    if (cp->GetEncryptionTechnique() != STANDARD)
        throw std::runtime_error("EncryptZeroWithTrace: only STANDARD technique supported");

    const auto& params = cp->GetElementParams();
    const auto& dgg    = cp->GetDiscreteGaussianGenerator();
    typename DCRTPoly::TugType tug;

    DCRTPoly u = (cp->GetSecretKeyDist() == GAUSSIAN)
                     ? DCRTPoly(dgg, params, Format::EVALUATION)
                     : DCRTPoly(tug, params, Format::EVALUATION);
    DCRTPoly e0(dgg, params, Format::EVALUATION);
    DCRTPoly e1(dgg, params, Format::EVALUATION);

    DCRTPoly b = publicKey->GetPublicElements()[0].Clone();
    DCRTPoly a = publicKey->GetPublicElements()[1].Clone();
    const uint32_t sizeQ  = b.GetParams()->GetParams().size();
    const uint32_t sizeQl = params->GetParams().size();
    if (sizeQl != sizeQ) {
        b.DropLastElements(sizeQ - sizeQl);
        a.DropLastElements(sizeQ - sizeQl);
    }
    b.SetFormat(Format::EVALUATION);
    a.SetFormat(Format::EVALUATION);

    const auto ns = cp->GetNoiseScale();
    DCRTPoly c0 = b * u + ns * e0;
    DCRTPoly c1 = a * u + ns * e1;

    DCRTPoly u_c(u), e0_c(e0), e1_c(e1);
    u_c.SetFormat(Format::COEFFICIENT);
    e0_c.SetFormat(Format::COEFFICIENT);
    e1_c.SetFormat(Format::COEFFICIENT);

    auto ct = std::make_shared<CiphertextImpl<DCRTPoly>>(publicKey);
    ct->SetElements({std::move(c0), std::move(c1)});
    ct->SetNoiseScaleDeg(1);
    ct->SetEncodingType(COEF_PACKED_ENCODING);  // matches MakeCoefPackedPlaintext
    return {std::move(ct), std::move(u_c), std::move(e0_c), std::move(e1_c)};
}

vector<complex<double>> fft(const vector<double>& m) {
    const size_t N = m.size();
    if (N % 2 != 0)
        throw invalid_argument("N must be even");

    const size_t n = N / 2;
    if (n & (n - 1))
        throw invalid_argument("n must be a power of two");

    vector<complex<double>> vals(n);
    for (size_t i = 0; i < n; ++i)
        vals[i] = complex<double>(m[i], m[i + n]);

    const uint32_t cyclOrder = static_cast<uint32_t>(2 * N);  // = 4n
    DiscreteFourierTransform::FFTSpecial(vals, cyclOrder);
    return vals;
}

vector<double> ifft(const vector<complex<double>>& zprime) {
    const size_t n = zprime.size();
    if (n & (n - 1))
        throw invalid_argument("n must be a power of two");

    vector<complex<double>> vals = zprime;
    const uint32_t cyclOrder = static_cast<uint32_t>(4 * n);
    DiscreteFourierTransform::FFTSpecialInv(vals, cyclOrder);

    vector<double> m(2 * n);
    for (size_t i = 0; i < n; ++i) {
        m[i]     = vals[i].real();
        m[i + n] = vals[i].imag();
    }
    return m;
}


// old
// complex<double> evalNormalizedNorm(DCRTPoly& a, size_t k_idx) {
//     a.SetFormat(Format::COEFFICIENT);
//     auto poly = a.CRTInterpolate(); 
//     auto Q = poly.GetModulus();
//     auto Qhalf = Q >> 1;
//     const size_t N = poly.GetLength();

//     long exp = 1;
//     for (long i=0; i < k_idx; ++i) {
//         exp = (exp * 5) % (2 * N); 
//     }

//     double theta = exp * M_PI / static_cast<double>(N);
//     complex<double> omega_k(std::cos(theta), std::sin(theta));

//     complex<double> acc(0.0, 0.0);
//     complex<double> omega_pow(1.0, 0.0);
    
//     for (size_t j = 0; j < N; ++j) {
//         auto v = poly[j];
//         double c_j;

//         if (v > Qhalf) {
//             auto neg_v = Q - v;
//             c_j = - static_cast<double>(neg_v.ConvertToDouble())/static_cast<double>(Q.ConvertToDouble());
//         } else {
//             c_j = static_cast<double>(v.ConvertToDouble())/static_cast<double>(Q.ConvertToDouble());
//         }

//         acc += c_j * omega_pow;
//         omega_pow *= omega_k;
//     }

//     return acc;
// }

// new
complex<double> evalNormalizedNorm(DCRTPoly& a, size_t k_idx) {
    a.SetFormat(Format::COEFFICIENT);
    auto poly = a.CRTInterpolate(); 
    auto Q = poly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = poly.GetLength();

    vector<double> coeff(N, 0.0);

    size_t sizeQ = a.GetParams()->GetParams().size();

    for (size_t i = 0; i < sizeQ; ++i) {

        const NativePoly& ai = a.GetElementAtIndex(i);
        NativeInteger qi = ai.GetModulus();

        BigInteger QHati = Q / BigInteger(qi);
        NativeInteger QHati_inv = QHati.ModInverse(qi).ConvertToInt();
        
        NativePoly a_tmp = ai;
        a_tmp = a_tmp * QHati_inv;
        a_tmp.SetFormat(COEFFICIENT);

        for (size_t k = 0; k < N; ++k) 
            coeff[k] += static_cast<double>(a_tmp.GetValues()[k].ConvertToDouble()) / static_cast<double>(qi.ConvertToDouble());

    }

    for (auto& x : coeff)
        x -= round(x);

    long exp = 1;
    for (long i=0; i < k_idx; ++i) {
        exp = (exp * 5) % (2 * N); 
    }

    double theta = exp * M_PI / static_cast<double>(N);
    complex<double> omega_k(std::cos(theta), std::sin(theta));

    complex<double> acc(0.0, 0.0);
    complex<double> omega_pow(1.0, 0.0);
    
    for (size_t j = 0; j < N; ++j) {
        acc += coeff[j] * omega_pow;
        omega_pow *= omega_k;
    }

    return acc;
}

// old
// vector<complex<double>> normalizedFFT(DCRTPoly& a) {

//     a.SetFormat(Format::COEFFICIENT);
//     auto poly = a.CRTInterpolate(); 
//     auto Q = poly.GetModulus();
//     auto Qhalf = Q >> 1;
//     const size_t N = poly.GetLength();

//     vector<double> coeff(N, 0.0);

//     for (size_t j = 0; j < N; ++j) {
//         auto v = poly[j];

//         if (v > Qhalf) {
//             auto neg_v = Q - v;
//             coeff[j] = - static_cast<double>(neg_v.ConvertToDouble()) / static_cast<double>(Q.ConvertToDouble());
//         } else {
//             coeff[j] = static_cast<double>(v.ConvertToDouble()) / static_cast<double>(Q.ConvertToDouble());
//         }

//     }

//     return fft(coeff);
// }

// new
vector<complex<double>> normalizedFFT(DCRTPoly& a) {

    a.SetFormat(Format::COEFFICIENT);
    auto poly = a.CRTInterpolate(); 
    auto Q = poly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = poly.GetLength();

    vector<double> coeff(N, 0.0);

    size_t sizeQ = a.GetParams()->GetParams().size();

    for (size_t i = 0; i < sizeQ; ++i) {

        const NativePoly& ai = a.GetElementAtIndex(i);
        NativeInteger qi = ai.GetModulus();

        BigInteger QHati = Q / BigInteger(qi);
        NativeInteger QHati_inv = QHati.ModInverse(qi).ConvertToInt();
        
        NativePoly a_tmp = ai;
        a_tmp = a_tmp * QHati_inv;
        a_tmp.SetFormat(COEFFICIENT);

        for (size_t k = 0; k < N; ++k) 
            coeff[k] += static_cast<double>(a_tmp.GetValues()[k].ConvertToDouble()) / static_cast<double>(qi.ConvertToDouble());

    }

    for (auto& x : coeff)
        x -= round(x);

    return fft(coeff);
}

double Log2BigInteger(const BigInteger& x) {
    if (x == BigInteger(0))
        return -std::numeric_limits<double>::infinity();

    usint msb = x.GetMSB();
    if (msb <= 53)
        return std::log2(x.ConvertToDouble());

    usint shift = msb - 53;
    BigInteger top = x >> shift;
    return std::log2(top.ConvertToDouble()) + static_cast<double>(shift);
}

pair<double, double> ScaleBigInteger(const BigInteger& x) {
    if (x == BigInteger(0))
        return {0.0, 0.0};

    usint msb = x.GetMSB();
    usint keepBits = std::min<usint>(msb, 53);
    usint shift = msb - keepBits;
    BigInteger top = x >> shift;

    double mantissa = std::ldexp(top.ConvertToDouble(), -static_cast<int>(keepBits - 1));
    double logScale = static_cast<double>(shift + keepBits - 1);
    return {mantissa, logScale};
}

double log2InfinityNorm(DCRTPoly& a) {
    a.SetFormat(Format::COEFFICIENT);
    auto poly = a.CRTInterpolate();
    auto Q = poly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = poly.GetLength();

    BigInteger maxAbs(0);
    for (size_t j = 0; j < N; ++j) {
        BigInteger centered = (poly[j] > Qhalf) ? (Q - poly[j]) : poly[j];
        if (centered.Compare(maxAbs) > 0)
            maxAbs = centered;
    }

    return Log2BigInteger(maxAbs);
}

double log2EvalNorm(DCRTPoly& a, size_t idx) {
    a.SetFormat(Format::COEFFICIENT);
    auto poly = a.CRTInterpolate();
    auto Q = poly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = poly.GetLength();

    long exp = 1;
    for (long i = 0; i < static_cast<long>(idx); ++i) {
        exp = (exp * 5) % (2 * static_cast<long>(N));
    }

    double theta = exp * M_PI / static_cast<double>(N);
    complex<double> omega_k(std::cos(theta), std::sin(theta));

    complex<double> acc(0.0, 0.0);
    double accLog2Scale = 0.0;
    constexpr double RESCALE_THRESHOLD = 1.0e200;
    constexpr double RESCALE_INV = 1.0 / RESCALE_THRESHOLD;
    constexpr double LOG2_RESCALE = 664.3856189774724;  // log2(1e200)

    for (int j = static_cast<int>(N) - 1; j >= 0; --j) {
        acc *= omega_k;

        BigInteger coeffAbs = (poly[j] > Qhalf) ? (Q - poly[j]) : poly[j];
        if (coeffAbs != BigInteger(0)) {
            auto [coeffMant, coeffScale] = ScaleBigInteger(coeffAbs);
            if (poly[j] > Qhalf)
                coeffMant = -coeffMant;

            double targetScale = std::max(accLog2Scale, coeffScale);
            acc *= std::ldexp(1.0, -static_cast<int>(targetScale - accLog2Scale));
            acc += complex<double>(std::ldexp(coeffMant, -static_cast<int>(targetScale - coeffScale)), 0.0);
            accLog2Scale = targetScale;
        }

        double accAbs = std::abs(acc);
        if (accAbs > RESCALE_THRESHOLD) {
            acc *= RESCALE_INV;
            accLog2Scale += LOG2_RESCALE;
        }
    }

    double accAbs = std::abs(acc);
    if (accAbs == 0.0)
        return -std::numeric_limits<double>::infinity();

    return accLog2Scale + std::log2(accAbs);
}


inline complex<double> evalIntPoly(const vector<int64_t>& coeffs, size_t N, size_t idx) {
    long exp_val = 1;
    for (size_t i = 0; i < idx; ++i)
        exp_val = (exp_val * 5) % (2 * static_cast<long>(N));
    double theta = exp_val * M_PI / static_cast<double>(N);
    complex<double> omega(std::cos(theta), std::sin(theta));
    complex<double> acc(0.0, 0.0), pow_om(1.0, 0.0);
    for (size_t j = 0; j < std::min(coeffs.size(), N); ++j) {
        acc    += static_cast<double>(coeffs[j]) * pow_om;
        pow_om *= omega;
    }
    return acc;
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
    noise        = log2InfinityNorm(res);
    noise_idx    = log2EvalNorm(res, idx);

}



inline double sample_uniform_dist(mt19937_64& gen, double B) {
    uniform_real_distribution<double> dist(-B, B);
    return dist(gen);
}

inline double sample_uniform_pm1(mt19937_64& gen) {
    uniform_int_distribution<int> dist(0, 1);
    return dist(gen) *2 - 1;
}



// ---------------------------------------------------------------------------
// a/q is a genuine (near-)uniform element of R_q, conditioned on its evaluation
// at omega_idx: this is what a brute-force search over sums of Enc(0) yields.
// It is built in the coefficient domain -- a uniform draw plus a single sinusoid
// moving (a/q)(omega_idx) to the target -- so that ||a/q||_inf stays close to
// 1/2 and every canonical coordinate remains populated.
inline NativePoly LargeApartBaseUniform(size_t n, size_t idx, NativeInteger q0,
                                        const shared_ptr<ILNativeParams>& params_q0, double kappa) {
    mt19937_64 gen(random_device{}());
    const size_t N = 2 * n;
    const double var = static_cast<double>(n) / 6.0;             // = N/12

    // target value of (a/q)(omega_idx), same law as before
    double rnd  = sample_uniform_dist(gen, kappa);
    double rest = sqrt(max(0.0, kappa * kappa - rnd * rnd));
    int    sgn  = sample_uniform_pm1(gen);
    complex<double> target(rnd * sqrt(var), sgn * rest * sqrt(var));

    // the correction is a single frequency; measure its sup-norm first
    vector<complex<double>> zc(n, {0.0, 0.0});
    zc[idx] = target;                                            // upper bound on |d|
    vector<double> probe = ifft(zc);
    double corrBound = 0.0;
    for (double x : probe) corrBound = max(corrBound, fabs(x));
    corrBound *= 1.6;                                            // |d| <= |target| + |u_unif(w)|
    const double B = 0.5 - corrBound;
    if (B <= 0.0) OPENFHE_THROW("LargeApartBaseUniform: kappa too large for this N");

    // Rejection sampling: draw uniform coefficients in [-B, B], shift the single
    // slot idx to the target, and keep the result only if it is a legitimate
    // element of R_q, that is ||a/q||_inf < 1/2. The 1.6 above is a heuristic
    // margin on the correction, so a draw overshoots now and then; redrawing
    // leaves the law of the accepted sample untouched.
    uniform_real_distribution<double> ud(-B, B);
    vector<double> u(N);
    constexpr int MAX_DRAWS = 1000;
    int draw = 0;
    double mx = 0.0;
    for (; draw < MAX_DRAWS; ++draw) {
        for (size_t i = 0; i < N; ++i) u[i] = ud(gen);

        // shift the single slot idx to the target, leaving all others untouched
        vector<complex<double>> cur = fft(u);
        zc.assign(n, {0.0, 0.0});
        zc[idx] = target - cur[idx];
        vector<double> corr = ifft(zc);
        for (size_t i = 0; i < N; ++i) u[i] += corr[i];

        mx = 0.0;
        for (double x : u) mx = max(mx, fabs(x));
        if (mx < 0.5) break;
    }
    if (draw == MAX_DRAWS)
        OPENFHE_THROW("LargeApartBaseUniform: ||a/q||_inf >= 0.5 in every draw");

    NativeVector coeffs(static_cast<uint32_t>(N), q0);
    const double q0_ld = static_cast<double>(q0.ConvertToInt());
    for (size_t i = 0; i < N; ++i) {
        int64_t rounded = static_cast<int64_t>(llround(static_cast<double>(u[i]) * q0_ld));
        if (rounded < 0) coeffs[i] = static_cast<NativeInteger::Integer>(q0) - static_cast<NativeInteger::Integer>(-rounded);
        else             coeffs[i] = static_cast<NativeInteger::Integer>(rounded);
    }
    NativePoly a_q0(params_q0);
    a_q0.SetValues(std::move(coeffs), Format::COEFFICIENT);
    a_q0.SetFormat(Format::COEFFICIENT);
    return a_q0;
}



struct ReRandomizeApartTrace {
    Ciphertext<DCRTPoly> ciphertext;
    complex<double> aOverQAtW;
    complex<double> bOverQAtW;
};

ReRandomizeApartTrace reRandomizeApartWithTrace(PrivateKey<DCRTPoly> privateKey, ConstCiphertext<DCRTPoly> ciphertext,
                                                usint ptm, uint32_t idx, double kappa, bool verbose = false){
    
    const auto cryptoParams = static_pointer_cast<CryptoParametersBFVRNS>(privateKey->GetCryptoParameters());

    const vector<DCRTPoly>& cv = ciphertext->GetElements();
    DCRTPoly s                 = privateKey->GetPrivateElement();

    auto paramsQ = cryptoParams->GetElementParams();
    auto paramsQP = cryptoParams->GetParamsQP();

    usint N = paramsQ->GetRingDimension();
    usint n = N / 2;

    BigInteger modulusQ = cryptoParams->GetElementParams()->GetModulus();
    NativeInteger q0 = cv[0].GetElementAtIndex(0).GetModulus();
    
    BigInteger QHat = modulusQ / BigInteger(q0);
    NativeInteger Qhat = QHat.Mod(q0).ConvertToInt();

    DCRTPoly aLarge(paramsQ, Format::COEFFICIENT, true);

    auto param_q0 = paramsQ->GetParams()[0];
    NativePoly a_q0 = LargeApartBaseUniform(n, idx, q0, param_q0, kappa);

    aLarge.SetElementAtIndex(0, Qhat * a_q0); 
    aLarge.SetFormat(Format::EVALUATION);

    // Keep aLarge in EVALUATION format for the multiplications below.
    DCRTPoly aLargeForDebug = aLarge;
    vector<complex<double>> evals0 = normalizedFFT(aLargeForDebug);

    // if (verbose) {
    //     cout << "Choose a-part such that (a/Q)(w_" << idx << ") = " << evals0[idx].real() << " + " << evals0[idx].imag() << "i"
    //          << ",  log|(a/Q)(w_" << idx << ")| = " << log2(abs(evals0[idx])) << endl;
    // }
    
    size_t sizeQl = cv[0].GetParams()->GetParams().size();
    size_t sizeQs = s.GetParams()->GetParams().size();

    size_t diffQl = sizeQs - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly b = cv[0];
    b.SetFormat(Format::EVALUATION);

    b += scopy * cv[1];
    b -= scopy * aLarge;

    // Compute (b/Q)(ω) for the recovery model. normalizedFFT mutates its
    // argument (CRT-interpolates → COEFFICIENT), so we work on a copy and
    // keep b in EVALUATION for the new ciphertext.
    DCRTPoly bForDebug = b;
    vector<complex<double>> evalsB = normalizedFFT(bForDebug);

    auto newCtxt = ciphertext->CloneEmpty();
    newCtxt->SetElements({std::move(b), std::move(aLarge)});
    newCtxt->SetNoiseScaleDeg(ciphertext->GetNoiseScaleDeg());
    return {newCtxt, evals0[idx], evalsB[idx]};

}

Ciphertext<DCRTPoly> reRandomizeApart(PrivateKey<DCRTPoly> privateKey, ConstCiphertext<DCRTPoly> ciphertext, usint ptm,
                                      uint32_t idx, double kappa, bool verbose = true){
    return reRandomizeApartWithTrace(privateKey, ciphertext, ptm, idx, kappa, verbose).ciphertext;
}



//////////////////////////////////////////////////////////////////////////////
///// Phase 2/3 helpers: the signature of the scaling map                 /////
//////////////////////////////////////////////////////////////////////////////
//
// A BFV ciphertext of mu has phase  b + a*s = Delta(mu) + e + q*I  over R,
// where Delta : R_t -> R is the scaling encoding.  Delta is not Z-linear
// unless t | q, and the defect
//
//     nu(mu) := q*mu - t*Delta(mu),        ||nu(mu)||_inf < t,
//
// satisfies nu(mu)_j = r*mu_j mod t with r := q mod t.  Consequently, for
// integers c_i whose sum C is a multiple of t, the ciphertext sum_i c_i*ct_i
// has tracked plaintext C*mu = 0 in R_t -- so it may legitimately be sent to
// the decryption oracle -- while its noise is
//
//     sum_i c_i*e_i  -  (C/t)*nu(mu) ,
//
// i.e. it carries a deterministic, message-dependent term.

// Signature of a single plaintext coefficient: the representative of r*m
// modulo t that carries the sign of m (this is what OpenFHE's per-coefficient
// Delta(mu)_j = sign(mu_j)*floor(q*|mu_j|/t) leaves behind).
inline int64_t SigOfCoeff(uint64_t r, uint64_t ptm, int64_t m) {
    if (m == 0) return 0;
    uint64_t a = static_cast<uint64_t>(m > 0 ? m : -m);
    int64_t  v = static_cast<int64_t>((static_cast<__uint128_t>(r) * a) % ptm);
    return (m > 0) ? v : -v;
}

// Plaintext maximising |nu(pt)(omega_idx)|.  This is the extremal problem of
// Appendix "Choice of Plaintext for Each Primitive Root of Unity", but solved
// in the nu-coordinate and pulled back through Sig (legal since r is invertible
// modulo t): we ask for the +-F sign pattern in nu, not in the plaintext.
// The resulting phase is F/sin(pi/(2N)) ~ 2*N*t/pi, the largest allowed by
// ||nu||_inf < t.  Writes the realised signature to sigOut when non-null.
inline Plaintext LargeSigPtxt(CryptoContext<DCRTPoly> cryptoContext, size_t idx,
                              uint64_t ptm, uint64_t r, vector<int64_t>* sigOut = nullptr) {
    const uint32_t N = cryptoContext->GetRingDimension();
    uint64_t F = 0, mstar = 1;
    for (uint64_t m = 1; m <= ptm / 2; ++m) {
        uint64_t v = (static_cast<__uint128_t>(r) * m) % ptm;
        if (v > F) { F = v; mstar = m; }
    }
    uint32_t g = 1;
    for (size_t i = 0; i < idx; ++i) g = (g * 5) % (2 * N);

    vector<int64_t> coeffs(N);
    if (sigOut) sigOut->assign(N, 0);
    for (uint32_t j = 0; j < N; ++j) {
        const bool pos = cos(M_PI * static_cast<double>(g) * (2.0 * j + 1.0) / (2.0 * N)) > 0;
        coeffs[j] = pos ? static_cast<int64_t>(mstar) : -static_cast<int64_t>(mstar);
        if (sigOut) (*sigOut)[j] = SigOfCoeff(r, ptm, coeffs[j]);
    }
    return cryptoContext->MakeCoefPackedPlaintext(coeffs);
}

// k * ct by double-and-add, using ciphertext additions only.  Note that this is
// *not* a plaintext multiplication: EvalMult by the plaintext t would encode t
// as 0 in R_t and destroy the ciphertext.
inline Ciphertext<DCRTPoly> IntegerMulCT(CryptoContext<DCRTPoly> cryptoContext,
                                         ConstCiphertext<DCRTPoly> ciphertext, uint64_t k) {
    Ciphertext<DCRTPoly> acc = nullptr;
    Ciphertext<DCRTPoly> cur = ciphertext->Clone();
    while (k) {
        if (k & 1) acc = acc ? cryptoContext->EvalAdd(acc, cur) : cur->Clone();
        k >>= 1;
        if (k) cur = cryptoContext->EvalAdd(cur, cur);
    }
    return acc;
}

inline DCRTPoly PolyFromInt64(const vector<int64_t>& v,
                              const shared_ptr<typename DCRTPoly::Params>& params) {
    DCRTPoly out(params, Format::COEFFICIENT, true);
    const size_t n = params->GetRingDimension();
    for (size_t i = 0; i < params->GetParams().size(); ++i) {
        NativeInteger qi = params->GetParams()[i]->GetModulus();
        NativeVector vec(static_cast<uint32_t>(n), qi);
        for (size_t j = 0; j < n; ++j) {
            const int64_t x = v[j];
            vec[j] = (x >= 0) ? NativeInteger(static_cast<uint64_t>(x)).Mod(qi)
                              : qi.ModSub(NativeInteger(static_cast<uint64_t>(-x)).Mod(qi), qi);
        }
        NativePoly p(params->GetParams()[i]);
        p.SetValues(std::move(vec), Format::COEFFICIENT);
        out.SetElementAtIndex(i, std::move(p));
    }
    return out;
}

// sum_{i<M} c_i * Enc(pt) with near-uniform c_i summing to C, a multiple of t.
// The tracked plaintext of the result is C*pt = 0 in R_t.
//
// exact = true  performs the M encryption queries, which is what the adversary
//               really does but costs M public-key encryptions;
// exact = false draws the aggregated randomness directly.  Since the sum is
//               sum_i c_i*u_i * pk + sum_i c_i*e_0i , sum_i c_i*e_1i , with the
//               u_i ternary uniform and the e_i discrete Gaussian, the aggregate
//               is the same law up to the central limit theorem, which is tight
//               for the values of M used here.  This mirrors the way phase 1
//               emulates the |a(omega)| search.
inline Ciphertext<DCRTPoly> EncryptLinComb(CryptoContext<DCRTPoly> cryptoContext,
                                           const KeyPair<DCRTPoly>& keyPair,
                                           Plaintext plaintext,
                                           uint32_t M, uint64_t C, bool exact = false) {
    const auto cp = static_pointer_cast<CryptoParametersBFVRNS>(keyPair.secretKey->GetCryptoParameters());
    const auto params = cp->GetElementParams();
    const size_t N = params->GetRingDimension();
    const double c = static_cast<double>(C) / static_cast<double>(M);

    if (exact) {
        Ciphertext<DCRTPoly> acc = nullptr;
        const uint64_t clo = static_cast<uint64_t>(C / M), rem = C - clo * M;
        for (uint32_t i = 0; i < M; ++i) {
            auto ci  = cryptoContext->Encrypt(keyPair.publicKey, plaintext);
            auto sci = IntegerMulCT(cryptoContext, ci, clo + (i < rem ? 1 : 0));
            acc = acc ? cryptoContext->EvalAdd(acc, sci) : sci;
        }
        return acc;
    }

    DCRTPoly bpk = keyPair.publicKey->GetPublicElements()[0].Clone();
    DCRTPoly apk = keyPair.publicKey->GetPublicElements()[1].Clone();
    const uint32_t szPk = bpk.GetParams()->GetParams().size();
    const uint32_t szQl = params->GetParams().size();
    if (szPk != szQl) { bpk.DropLastElements(szPk - szQl); apk.DropLastElements(szPk - szQl); }
    bpk.SetFormat(Format::EVALUATION); apk.SetFormat(Format::EVALUATION);

    mt19937_64 gen(random_device{}());
    const double sigma = cp->GetDistributionParameter();
    auto draw = [&](double sd) {
        normal_distribution<double> nd(0.0, sd);
        vector<int64_t> v(N);
        for (size_t j = 0; j < N; ++j) v[j] = llround(nd(gen));
        return PolyFromInt64(v, params);
    };
    DCRTPoly U  = draw(c * sqrt(2.0 * M / 3.0));  U.SetFormat(Format::EVALUATION);
    DCRTPoly E0 = draw(c * sigma * sqrt((double)M)); E0.SetFormat(Format::EVALUATION);
    DCRTPoly E1 = draw(c * sigma * sqrt((double)M)); E1.SetFormat(Format::EVALUATION);

    auto zeroSum = make_shared<CiphertextImpl<DCRTPoly>>(keyPair.publicKey);
    zeroSum->SetElements({bpk * U + E0, apk * U + E1});
    zeroSum->SetNoiseScaleDeg(1);
    zeroSum->SetEncodingType(COEF_PACKED_ENCODING);

    // add C * Delta(pt), noiselessly
    auto shift = cryptoContext->EvalSub(cryptoContext->EvalAdd(zeroSum, plaintext), zeroSum);
    return cryptoContext->EvalAdd(zeroSum, IntegerMulCT(cryptoContext, shift, C));
}

//////////////////////
///// New method /////
//////////////////////
