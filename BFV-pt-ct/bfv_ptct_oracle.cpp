// The attack of Algorithm 1, using only what an IND-CPA-D adversary has: for each
// root w_k, K pt-ct multiplications and one decryption query, with K fixed by the
// heuristic bound (d_leg). The recovered key is checked against the public key.
// The secret key is used by the decryption oracle and, after the attack, to report
// the match ratio. bfv_ptct_sk_general is the version we used to tune the attack:
// it also reads the noise with the secret key, to pick and check parameters.
#include "bfv_ptct.hpp"
#include <complex>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <limits>
#include <vector>
#include <map>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

using namespace lbcrypto;
using namespace std;

typedef std::complex<double> ComplexLD;

/* --------------
--- FUNCTIONS ---
-------------- */
ComplexLD evalComplexValue(DCRTPoly& a, size_t k_idx) {
    a.SetFormat(Format::COEFFICIENT);
    auto poly = a.CRTInterpolate(); 
    auto Q = poly.GetModulus();
    auto Qhalf = Q >> 1;
    const size_t N = poly.GetLength();

    double theta = (2.0 * k_idx - 1.0) * M_PI / static_cast<double>(N);
    ComplexLD omega_k(std::cos(theta), std::sin(theta));
    ComplexLD acc(0.0, 0.0);

    for (int j = static_cast<int>(N) - 1; j >= 0; --j) {
        auto v = poly[j];
        double c_j = (v > Qhalf) ? 
            -static_cast<double>((Q - v).ConvertToDouble()) : 
             static_cast<double>(v.ConvertToDouble());
        acc = (acc * omega_k) + c_j;
    }
    return acc;
}

double GetVnorm(CryptoContext<DCRTPoly> cryptoContext) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersBFVRNS>(cryptoContext->GetCryptoParameters());
    if (!cryptoParams) {
        throw std::runtime_error("Invalid crypto parameters: Expected BFVRNS.");
    }

    // Extract parameters from the context
    uint32_t n = cryptoContext->GetRingDimension();
    double sigma = cryptoParams->GetDistributionParameter();
    double alpha = cryptoParams->GetAssuranceMeasure();
    double Berr = sigma * std::sqrt(alpha);
    
    // Determine Bkey based on the secret key distribution
    double Bkey = 0;
    uint32_t thresholdParties = cryptoParams->GetThresholdNumOfParties();
    
    if (cryptoParams->GetSecretKeyDist() == GAUSSIAN) {
        Bkey = std::sqrt(thresholdParties) * Berr;
    } else {
        // Default for UNIFORM_TERNARY and other non-Gaussian distributions
        Bkey = static_cast<double>(thresholdParties);
    }

    // Expansion factors defined in bfvrns-parametergeneration.cpp
    auto delta = [](uint32_t n_val) -> double {
        return 4. * std::sqrt(static_cast<double>(n_val));
    };

    auto deltaMS = [](uint32_t n_val) -> double {
        return 4. * std::sqrt(static_cast<double>(n_val));
    };

    EncryptionTechnique encTech = cryptoParams->GetEncryptionTechnique();

    // Replicate Vnorm logic from the library source
    if (encTech == EXTENDED) {
        return (1. + deltaMS(n) * Bkey) / 2.;
    } else {
        return Berr * (1. + 2. * delta(n) * Bkey);
    }
}

int main(int argc, char* argv[]) {

    // Default values
    int logN = 12;
    int64_t ptm = 786433;
    // "SetDepth" is used for roughly choosing Q.
    // There is no other tool in OpenFHE to set a custom-length Q. */
    int SetDepth = 6;
    double customFactor = 1.0;
    // "ct0Num" is used to rerandomize ct0 (but still with fixed seeds for the same "ct0Num").
    int ct0Num = 1;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--logN" && i + 1 < argc) {
            logN = std::stoi(argv[++i]);
        } else if (arg == "--ptm" && i + 1 < argc) {
            ptm = std::stoll(argv[++i]);
        } else if (arg == "--depth" && i + 1 < argc) {
            SetDepth = std::stoi(argv[++i]);
        } else if (arg == "--cF" && i + 1 < argc) {
            customFactor = std::stod(argv[++i]);
        } else if (arg == "--ct0Num" && i + 1 < argc) {
            ct0Num = std::stoi(argv[++i]);
        } else if (arg == "--help") {
            std::cout << "Usage: [OMP_NUM_THREADS=1] " << argv[0] 
                      << " [--logN N] [--ptm t] [--depth SetDepth] [--cF customFactor (1.0)] [--ct0Num ct0Num (optional)]" << std::endl;
            return 0;
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << std::endl;
            return 1;
        }
    }

    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetRingDim(1 << logN);
    parameters.SetPlaintextModulus(ptm);
    parameters.SetMultiplicativeDepth(SetDepth);
    parameters.SetMultiplicationTechnique(HPSPOVERQLEVELED);
    parameters.SetSecurityLevel(SecurityLevel::HEStd_NotSet);

    // --- WORKING PARAMS --- //
    /** 
    cmake -S . -B build -DWITH_FIXED_SEED=ON
    cmake --build build -j
    OMP_NUM_THREADS=1 
    // 12, 786433, 4, 240
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 12 --ptm 786433 --depth 4 --cF 1.0
    // 12, 8388609, 5, 300
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 12 --ptm 8388609 --depth 5 --cF 1.0
    // ----- //
    // 13, 4194305, 4, 240
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 13 --ptm 4194305 --depth 4 --cF 1.0
    // 13, 8388609, 10, 480
    // 13, 270794753, 6 ,360
    // ----- //
    // 14, 67108865, 4, 240
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 14 --ptm 67108865 --depth 4 --cF 1.0
    // 14, 2097153, 5, 300
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 14 --ptm 2097153 --depth 5 --cF 1.0
    // ----- //
    // 15, 335544305, 4, 240
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 15 --ptm 33423361 --depth 4 --cF 1.0
    // 15, 270794753, 4, 300
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 15 --ptm 270794753 --depth 4 --cF 1.0
    // ----- //
    // 16, 134086657, 4, 300
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 16 --ptm 134086657 --depth 4 --cF 1.0
    // 16, 270794753, 4, 300
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 16 --ptm 270794753 --depth 4 --cF 1.0
    // 16, 1048577, 4, 240
    ./build/bin/examples/pke/bfv_ptct_sk_general --logN 16 --ptm 1048577 --depth 4 --cF 1.0
    */

    // --- WORKING PARAMS, PROBABILISTIC --- //
    // Probabilistically works, but not for our seed = 1. 
    // ./build/bin/examples/pke/bfv_ptct_sk_general --logN 13 --ptm 65537 --depth 4 --cF 1.0
    // ./build/bin/examples/pke/bfv_ptct_sk_general --logN 13 --ptm 1048577 --depth 5 --cF 1.0
    // ./build/bin/examples/pke/bfv_ptct_sk_general --logN 13 --ptm 4194305 --depth 5 --cF 1.0
    // ./build/bin/examples/pke/bfv_ptct_sk_general --logN 16 --ptm 786433 --depth 4 --cF 1.0
    // ./build/bin/examples/pke/bfv_ptct_sk_general --logN 16 --ptm 33423361 --depth 4 --cF 1.0
    // ./build/bin/examples/pke/bfv_ptct_sk_general --logN 16 --ptm 262145 --depth 5 --cF 1.0
    // ./build/bin/examples/pke/bfv_ptct_sk_general --logN 16 --ptm 65537 --depth 5 --cF 1.0
    // ---------------------- //

    // ------------- //
    // Parameter Gen //
    // ------------- //

    double normEXP_e1_inf;
    switch (logN) {
        case 12:
            normEXP_e1_inf = 16.1807;
            break;
        case 13:
            normEXP_e1_inf = 17.3890;
            break;
        case 14:
            normEXP_e1_inf = 18.1876;
            break;
        case 15:
            normEXP_e1_inf = 19.2704;
            break;
        case 16:
            normEXP_e1_inf = 20.2343;
            break;
        default:
            throw std::runtime_error("Invalid logN choice. Please select between 12 and 16.");
    }

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);

    // ------------------ //
    // Context generation //
    // ------------------ //
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);
    
    const auto cryptoParams = std::static_pointer_cast<CryptoParametersBFVRNS>(cryptoContext->GetCryptoParameters());
    double logQ = 0;
    for (auto& param : cryptoParams->GetElementParams()->GetParams()) {
        logQ += log2(param->GetModulus().ConvertToDouble());
    }

    // ----------- //
    // Initial ct0 //
    // ----------- //
    usint N = cryptoContext->GetRingDimension();
    std::vector<int64_t> vectorOfInts1(N, 0);
    Plaintext plaintext0 = cryptoContext->MakeCoefPackedPlaintext(vectorOfInts1);
    auto ct0 = cryptoContext->Encrypt(keyPair.publicKey, plaintext0);
    // For better seed... LOL
    for (int i = 1; i < ct0Num; ++i) {
        ct0 = cryptoContext->Encrypt(keyPair.publicKey, plaintext0);
    }

    // Parameter Summary
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "CONFIGURATION SUMMARY" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Ring Dimension (N)    : " << N << std::endl;
    std::cout << "Log2 Q (Ciphertext)   : " << std::fixed << std::setprecision(2) << logQ << std::endl;
    std::cout << "Plaintext Modulus (p) : " << ptm << " (~" << std::log2((double)ptm) << " bits)" << std::endl;
    std::cout << "SetMultiplicativeDepth: " << SetDepth << std::endl;

    // Depth based on estimate
    // Q/(2t) > (delta_N * t / 2)^K * Vnorm 
    // delta_N = 2 * sqrt(N) & Vnorm as in OpenFHE: bfvrns-parametergeneration.cpp
    double vNorm = GetVnorm(cryptoContext);
    double delta = 2.0 * std::sqrt(static_cast<double>(N));
    uint32_t depth_estimate = static_cast<uint32_t>(std::floor(  (logQ - 1.0 - log2((double)ptm) - log2(vNorm)) / (log2((double)ptm) - 1.0 + log2(delta))  ));
    if (depth_estimate == 0) {
        throw std::runtime_error("depth_estimate must be > 0 for e1-based constant formula.");
    }
    
    // Rough estimate on the constant required for attack (e1-based, c=1 profiled)
    // Q / 2 ~ (c * N / pi)^(K-1) * (c * 2^(normEXP_e1_inf))
    //       = c^K * (N / pi)^(K-1) * 2^(normEXP_e1_inf)
    // => log2(c) = [log2(Q) - 1 - normEXP_e1_inf + (K-1)*(log2(pi)-log2(N))] / K
    double constant_factor = (logQ - 1.5 - normEXP_e1_inf +
                              (static_cast<double>(depth_estimate) - 1.0) *
                                  (log2((double)M_PI) - log2(static_cast<double>(N)))) /
                             static_cast<double>(depth_estimate);
    int64_t constant = static_cast<int64_t>(pow(2.0, constant_factor) * customFactor);
    if (static_cast<double>(constant) > static_cast<double>(ptm) / 2.0) {
        std::cout << "!!!!!! The constant = " << constant << " is too large (>t/2) !!!!!!" << std::endl;
        constant = std::min(constant, static_cast<int64_t>(ptm) >> 1);
    }
    // We now use exactly K depths. 
    // uint32_t depth_estimate_attack = static_cast<uint32_t>(std::floor((logQ - 1 - log2(vNorm))) / (log2(constant) + log2(static_cast<double>(N)) - log2(M_PI))) + 1;
    
    std::cout << "Vnorm = " << std::setprecision(4) << log2(vNorm) << std::endl;
    std::cout << "normEXP_e1_inf = " << normEXP_e1_inf << " (from e1_inf)" << std::endl;
    std::cout << "Estimated Factor = " << (log2((double)ptm) - 1.0 + log2(delta)) << " (delta_N * t / 2), delta_N = 2 * sqrt(N)." << std::endl;
    std::cout << "Theoretical Est. Depth: " << depth_estimate << std::endl;
    std::cout << "Constant used for attack: " << constant << " (2^" << constant_factor << ")" << std::endl;
    std::cout << "Custom Factor: " << customFactor << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    // The decryption oracle: the only place where the attack touches the secret key.
    size_t oracle_calls = 0;
    auto oracleDec = [&](ConstCiphertext<DCRTPoly> ct) {
        Plaintext pt;
        cryptoContext->Decrypt(keyPair.secretKey, ct, &pt);
        oracle_calls++;
        return pt->GetCoefPackedValue();
    };

    // --------------------------------------------- //
    // Attack (Algorithm 1): no secret key used here //
    // --------------------------------------------- //
    // For each root w_k, K = depth_estimate pt-ct multiplications by mu_k, then one
    // decryption query; e0(w_k) is read off the decrypted result divided by mu_k(w_k)^K.
    std::cout << "--- Starting Error Recovery (Indices 1 to " << N << "), K = " << depth_estimate << " ---" << std::endl;
    auto collect_start = std::chrono::high_resolution_clock::now();
    std::vector<ComplexLD> estimated_e0_list(N);
    size_t log_interval = (N >> 3);
    for (size_t idx = 1; idx <= N; ++idx) {
        auto ptxt = LargePtxtApp(cryptoContext, idx, constant);
        DCRTPoly plain = ptxt->GetElement<DCRTPoly>();
        ComplexLD pt_w = evalComplexValue(plain, idx);

        auto ct = ct0->Clone();
        ComplexLD pt_w_pow(1.0, 0.0);
        for (uint32_t d = 0; d < depth_estimate; ++d) {
            ct = cryptoContext->EvalMult(ct, ptxt);
            pt_w_pow *= pt_w;
        }
        std::vector<int64_t> m_vec = oracleDec(ct);

        double phi = (2.0 * idx - 1.0) * M_PI / static_cast<double>(N);
        ComplexLD w(std::cos(phi), std::sin(phi));
        ComplexLD m_dec_w(0.0, 0.0);
        for (int j = static_cast<int>(m_vec.size()) - 1; j >= 0; --j) {
            m_dec_w = (m_dec_w * w) + static_cast<double>(m_vec[j]);
        }
        double Q_curr = 1.0;
        for (auto& p : ct->GetElements()[0].GetParams()->GetParams()) Q_curr *= p->GetModulus().ConvertToDouble();
        estimated_e0_list[idx - 1] = m_dec_w * (Q_curr / (double)ptm) / pt_w_pow;

        if (idx % log_interval == 0) {
            std::cout << "[Progress " << idx / log_interval << "/8] Index " << idx << " | Dec Calls (so far): " << oracle_calls << std::endl;
        }
    }
    auto collect_end = std::chrono::high_resolution_clock::now();

    // e0 by the inverse DFT, then s = (e0 - c0) * c1^{-1} mod q, rounded to {-1, 0, 1}.
    auto recovery_start = std::chrono::high_resolution_clock::now();
    std::vector<int64_t> e0_est_coeffs(N);
    for (size_t j = 0; j < N; ++j) {
        ComplexLD sum_val(0, 0);
        for (size_t k = 1; k <= N; ++k) {
            double theta = (2.0 * k - 1.0) * M_PI * j / static_cast<double>(N);
            sum_val += estimated_e0_list[k - 1] * ComplexLD(std::cos(theta), -std::sin(theta));
        }
        e0_est_coeffs[j] = std::round(sum_val.real() / static_cast<double>(N));
    }

    const vector<DCRTPoly>& cv = ct0->GetElements();
    DCRTPoly c0 = cv[0], c1 = cv[1];
    c0.SetFormat(Format::EVALUATION);
    c1.SetFormat(Format::EVALUATION);
    if (!c1.InverseExists()) {
        std::cout << "[WARN] c1 is not invertible in R_q; no key recovered." << std::endl;
        return 0;
    }
    DCRTPoly e_est(cv[0].GetParams(), Format::COEFFICIENT, true);
    e_est = e0_est_coeffs;
    e_est.SetFormat(Format::EVALUATION);
    DCRTPoly s_raw = (e_est - c0) * c1.MultiplicativeInverse();
    s_raw.SetFormat(Format::COEFFICIENT);
    auto s_raw_int = s_raw.CRTInterpolate();
    auto s_raw_mod = s_raw_int.GetModulus();
    std::vector<int64_t> sk_guess(N, 0);
    for (size_t j = 0; j < N; ++j) {
        auto v = s_raw_int[j];
        if (v == 0) continue;
        sk_guess[j] = (v > (s_raw_mod >> 1)) ? -1 : 1;
    }

    // The attacker's own check, from public data only: pk0 + pk1 * s_guess must be small.
    const vector<DCRTPoly>& pk = keyPair.publicKey->GetPublicElements();
    DCRTPoly s_guess(pk[0].GetParams(), Format::COEFFICIENT, true);
    s_guess = sk_guess;
    s_guess.SetFormat(Format::EVALUATION);
    DCRTPoly pk0 = pk[0], pk1 = pk[1];
    pk0.SetFormat(Format::EVALUATION);
    pk1.SetFormat(Format::EVALUATION);
    DCRTPoly pk_res = pk0 + pk1 * s_guess;
    pk_res.SetFormat(Format::COEFFICIENT);
    double log2_pk_res = log2(pk_res.Norm());
    bool pk_check = log2_pk_res < logQ / 2.0;
    auto recovery_end = std::chrono::high_resolution_clock::now();

    // ------------------------------------------------- //
    // Evaluation only: compare with the true secret key //
    // ------------------------------------------------- //
    DCRTPoly s_true = keyPair.secretKey->GetPrivateElement();
    size_t sizeQl = cv[0].GetParams()->GetParams().size();
    size_t sizeQs = s_true.GetParams()->GetParams().size();
    if (sizeQs > sizeQl) s_true.DropLastElements(sizeQs - sizeQl);
    s_true.SetFormat(Format::COEFFICIENT);
    auto s_true_int = s_true.CRTInterpolate();
    auto s_mod = s_true_int.GetModulus();
    size_t sk_match_cnt = 0;
    for (size_t j = 0; j < N; ++j) {
        auto v = s_true_int[j];
        int64_t centered = (v > (s_mod >> 1)) ? -static_cast<int64_t>((s_mod - v).ConvertToInt())
                                              : static_cast<int64_t>(v.ConvertToInt());
        if (centered == sk_guess[j]) sk_match_cnt++;
    }
    double sk_coeff_match_ratio = static_cast<double>(sk_match_cnt) / static_cast<double>(N);

    std::chrono::duration<double> time_collect = collect_end - collect_start;
    std::chrono::duration<double> time_recovery = recovery_end - recovery_start;
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "SUMMARY (oracle-only attack)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Multiplications (K)     : " << depth_estimate << " for every index" << std::endl;
    std::cout << "O_Dec calls             : " << oracle_calls << std::endl;
    std::cout << "Public-key check        : " << (pk_check ? "passed" : "FAILED")
              << " (log2 ||pk0 + pk1*s_guess||_inf = " << std::setprecision(4) << log2_pk_res << ")" << std::endl;
    std::cout << "SK Coeff Match Ratio    : " << std::fixed << std::setprecision(6) << sk_coeff_match_ratio
              << std::defaultfloat << "   (evaluation only, uses the true key)" << std::endl;
    std::cout << "Collecting Time (e(w))  : " << time_collect.count() << " s" << std::endl;
    std::cout << "Recovery Time (e, s)    : " << time_recovery.count() << " s" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    return 0;
}
