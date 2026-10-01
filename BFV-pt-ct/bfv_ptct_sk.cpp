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

void EvalNoiseBFVComplex(PrivateKey<DCRTPoly> privateKey, ConstCiphertext<DCRTPoly> ciphertext, 
                         double& noise_inf, ComplexLD& noise_w_val, size_t idx) {
    const vector<DCRTPoly>& cv = ciphertext->GetElements();
    DCRTPoly s = privateKey->GetPrivateElement();
    size_t sizeQl = cv[0].GetParams()->GetParams().size();
    size_t sizeQs = s.GetParams()->GetParams().size();
    size_t diffQl = sizeQs - sizeQl;

    auto scopy(s);
    if (diffQl > 0) scopy.DropLastElements(diffQl);

    DCRTPoly sPower(scopy);
    DCRTPoly b = cv[0];
    b.SetFormat(Format::EVALUATION);

    for (size_t i = 1; i < cv.size(); i++) {
        DCRTPoly ci = cv[i];
        ci.SetFormat(Format::EVALUATION);
        b += (sPower * ci);
        if (i < cv.size() - 1) sPower *= scopy;
    }

    b.SetFormat(Format::COEFFICIENT);
    noise_inf = log2(b.Norm()); 
    noise_w_val = evalComplexValue(b, idx); 
}

ComplexLD AnalyzeAndRecoverInitialError(
    size_t idx, 
    Ciphertext<DCRTPoly> ct0, 
    CryptoContext<DCRTPoly> cryptoContext, 
    KeyPair<DCRTPoly> keyPair,
    uint64_t ptm,
    double logQ, 
    int64_t constant,
    uint32_t& last_iter,
    int depthNum, 
    bool verbose = false,
    size_t* evalnoise_call_counter = nullptr) 
{
    auto ptxt = LargePtxtApp(cryptoContext, idx, constant); 
    DCRTPoly plain = ptxt->GetElement<DCRTPoly>();
    ComplexLD pt_w = evalComplexValue(plain, idx);
    double log_pt_w = log2(std::abs(pt_w));

    double noise_inf = 0;
    ComplexLD e_0_w; 
    EvalNoiseBFVComplex(keyPair.secretKey, ct0, noise_inf, e_0_w, idx);
    if (evalnoise_call_counter) {
        (*evalnoise_call_counter)++;
    }
    
    auto current_ct = ct0->Clone();
    ComplexLD e_w_prev = e_0_w;
    ComplexLD pt_w_pow(1.0, 0.0); 
    
    Ciphertext<DCRTPoly> last_valid_ct;
    ComplexLD last_actual_e_w, last_pt_pow_w;
    last_iter = 0;

    if (verbose) {
        std::cout << " ... PT info ... " << std::defaultfloat << std::setprecision(5) << std::endl;
        std::cout << "||pt||_inf: 2^" << log2(static_cast<double>(constant)) << std::endl;
        std::cout << "|pt(w)|   : 2^" << log_pt_w << std::endl;
        std::cout << std::string(60, '-') << std::endl;
    }

    for (uint32_t i = 1; i <= depthNum; i++){ 
        ComplexLD expected_e_curr = e_0_w * pt_w_pow * pt_w;
        pt_w_pow *= pt_w; 

        current_ct = cryptoContext->EvalMult(current_ct, ptxt);

        // O(N) decryption queries. Experimentally, 2N to 5N is enough. 
        if (verbose || (i > depthNum - 4)) {
            // For debugging..
            if (verbose) {
                std::cout << "pt-ct mult done " << i << " times for idx " << idx << std::endl;
                std::cout << "EvalNoiseBFVComplex running.. " << std::endl;
            }

            // Decryption Oracle outputing the scaled and rounded results. 
            EvalNoiseBFVComplex(keyPair.secretKey, current_ct, noise_inf, e_w_prev, idx); 
            if (evalnoise_call_counter) {
                (*evalnoise_call_counter)++;
            }

            double rel_err_curr = std::abs(e_w_prev - expected_e_curr) / std::abs(expected_e_curr);

            // For debugging..
            if (verbose) {
                std::cout << "real e_i (L-inf):  " << std::scientific << noise_inf << std::endl;
                // std::cout << "expected e_i(w) :  " << expected_e_curr << std::endl;
                std::cout << "real e_i(w)     :  " << e_w_prev << std::endl;
                // std::cout << "relative error for e(w): " << rel_err_curr << std::defaultfloat << std::endl;
            }

            if (verbose || (!verbose && rel_err_curr < 0.01)) {
                last_iter = i;
                last_valid_ct = current_ct->Clone();
                last_actual_e_w = e_w_prev;
                last_pt_pow_w = pt_w_pow;
            } else {
                break;
            }
        }
        
        if (noise_inf > logQ - 1) break;
    }

    if (last_valid_ct) {
        Plaintext last_dec_pt;
        cryptoContext->Decrypt(keyPair.secretKey, last_valid_ct, &last_dec_pt);
        std::vector<int64_t> m_vec = last_dec_pt->GetCoefPackedValue();

        const size_t N_dim = cryptoContext->GetRingDimension();
        double phi = (2.0 * idx - 1.0) * M_PI / static_cast<double>(N_dim);
        ComplexLD w(std::cos(phi), std::sin(phi));

        ComplexLD m_dec_w(0.0, 0.0);
        for (int j = static_cast<int>(m_vec.size()) - 1; j >= 0; --j) {
            m_dec_w = (m_dec_w * w) + static_cast<double>(m_vec[j]);
        }

        auto current_params = last_valid_ct->GetElements()[0].GetParams()->GetParams();
        double Q_curr = 1.0;
        for(auto& p : current_params) Q_curr *= (double)p->GetModulus().ConvertToDouble();

        ComplexLD recovered_e_K_w = m_dec_w * (Q_curr / (double)ptm);

        // e_initial_estimated
        return recovered_e_K_w / last_pt_pow_w;
    }

    return ComplexLD(0, 0);
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

/* --------------
--- CANDIDATE ---
-------------- */
// ------ Candidate plaintext moduli ------ //
// uint64_t ptm =      65537; // 16.00 bits (Existing)
// uint64_t ptm =     262145; // 18.00 bits
// uint64_t ptm =     786433; // 19.58 bits (Existing)
// uint64_t ptm =    1048577; // 20.00 bits
// uint64_t ptm =    2097153; // 21.00 bits
// uint64_t ptm =    4194305; // 22.00 bits
// uint64_t ptm =    8388609; // 23.00 bits
// uint64_t ptm =   16777217; // 24.00 bits
// uint64_t ptm =   67108865; // 26.00 bits
// uint64_t ptm =  270794753; // 28.01 bits (Existing)

// (logN, ptm,  logC, d)
// (12, 16777217,  3, 2)
// (12, 16777217,  8, 2)
// (13, 8388609,   3, 2)
// (13, 8388609,   4, 2)
// (14, 786433,    3, 2)
// (15, 4194305,   2, 2)
// (15, 4194305,   3, 4)
// (16, 270794753, 3, 2)

int main(int argc, char* argv[]) {
    // Default values
    int logN = 12;
    uint64_t ptm = 786433;
    int d = 4;
    int logC = 8;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--logN" && i + 1 < argc) {
            logN = std::stoi(argv[++i]);
        } else if (arg == "--ptm" && i + 1 < argc) {
            ptm = static_cast<uint64_t>(std::stoull(argv[++i]));
        } else if (arg == "--logC" && i + 1 < argc) {
            logC = std::stoi(argv[++i]);
        } else if (arg == "--depth" && i + 1 < argc) {
            d = std::stoi(argv[++i]);
        } else if (arg == "--help") {
            std::cout << "Usage: [OMP_NUM_THREADS=1] " << argv[0]
                      << " [--logN N] [--ptm t] [--logC log2||pt||] [--depth SetDepth]" << std::endl;
            return 0;
        } else {
            std::cerr << "Unknown or incomplete argument: " << arg << std::endl;
            return 1;
        }
    }

    if (logN < 12 || logN > 16) {
        std::cerr << "Invalid --logN: expected 12..16" << std::endl;
        return 1;
    }
    if (d <= 0) {
        std::cerr << "Invalid --depth: must be > 0" << std::endl;
        return 1;
    }
    if (logC < 0 || logC > 60) {
        std::cerr << "Invalid --logC: expected 0..60" << std::endl;
        return 1;
    }

    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetPlaintextModulus(ptm);
    parameters.SetMultiplicationTechnique(HPSPOVERQLEVELED);
    parameters.SetMultiplicativeDepth(d);
    parameters.SetSecurityLevel(SecurityLevel::HEStd_NotSet);
    parameters.SetRingDim(1 << logN);

    int64_t constant = static_cast<int64_t>(1) << logC;

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);

    // ------------------ //
    // Context generation //
    // ------------------ //
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);
    
    const auto cryptoParams = std::static_pointer_cast<CryptoParametersBFVRNS>(keyPair.secretKey->GetCryptoParameters());
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

    // Collecting e0(w)
    std::vector<ComplexLD> estimated_e0_list(N);
    double total_abs_error = 0;
    double max_abs_error = 0;
    size_t main_loop_evalnoise_calls = 0;
    uint32_t min_last_iter_main = std::numeric_limits<uint32_t>::max();
    uint32_t max_last_iter_main = 0;
    size_t log_interval = (N >> 3);

    // Parameter Summary
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "CONFIGURATION SUMMARY" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Ring Dimension (N)   : " << N << std::endl;
    std::cout << "Log2 Q (Ciphertext)  : " << std::fixed << std::setprecision(2) << logQ << std::endl;
    std::cout << "Plaintext Modulus (p): " << ptm << " (~" << std::log2((double)ptm) << " bits)" << std::endl;

    // Depth based on estimate
    // Q/(2t) > (delta_N * ||pt||_inf)^K * Vnorm 
    // delta_N = 2 * sqrt(N) & Vnorm as in OpenFHE: bfvrns-parametergeneration.cpp
    double vNorm = GetVnorm(cryptoContext);
    double delta = 2.0 * std::sqrt(static_cast<double>(N));
    uint32_t depth_estimate = static_cast<uint32_t>(std::floor(  (logQ - 1.0 - log2((double)ptm) - log2(vNorm)) / (log2((double)constant) - 1.0 + log2(delta))  ));
    
    std::cout << "Vnorm = " << log2(vNorm) << std::endl;
    std::cout << "Theoretical Est. Depth: " << depth_estimate << std::endl;

    // Rough estimate on the depth for attack
    // Q/(2) < ||e_K|| ~ |pt(w)|^(K-1) * ||e_1|| ~ |pt(w)|^K *Vnorm 
    uint32_t depth_estimate_attack = static_cast<uint32_t>(std::floor((logQ - 1 - log2(vNorm))) / (log2(constant) + log2(static_cast<double>(N)) - log2(M_PI))) + 1;
    std::cout << "Theoretical ATK. Depth: " << depth_estimate_attack << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    // ----------------- //
    // Pre-run for DEBUG //
    // ----------------- //
    // Pre-run a single index for verbose debugging before the main loop.
    const bool run_single_index_debug = true;
    const size_t debug_idx = 7;
    if (run_single_index_debug && debug_idx >= 1 && debug_idx <= N) {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "PRE-RUN SINGLE INDEX DEBUG (idx=" << debug_idx << ")" << std::endl;
        std::cout << std::string(60, '-') << std::endl;

        std::cout << "... Initial ct0 info ... " << std::endl;
        double debug_inf;
        ComplexLD debug_actual_e0_w;
        EvalNoiseBFVComplex(keyPair.secretKey, ct0, debug_inf, debug_actual_e0_w, debug_idx);
        std::cout << "real e_i (L-inf):  " << std::scientific << debug_inf << std::endl;
        std::cout << "real e_i(w)     :  " << debug_actual_e0_w << std::endl;
        std::cout << std::string(60, '-') << std::endl;

        uint32_t debug_iter = 0;
        ComplexLD debug_est_e0 = AnalyzeAndRecoverInitialError(
            debug_idx,
            ct0,
            cryptoContext,
            keyPair,
            ptm,
            logQ,
            constant,
            debug_iter,
            depth_estimate,
            true,
            nullptr);

        double debug_abs_error = std::abs(debug_est_e0 - debug_actual_e0_w);

        std::cout << "[Pre-run Debug] Iter used     : " << debug_iter << std::endl;
        std::cout << "[Pre-run Debug] abs error     : " << std::scientific << debug_abs_error << std::defaultfloat << std::endl;
        std::cout << std::string(60, '=') << std::endl;
    }
    
    // --------------//
    // Main for loop //
    // --------------//
    std::cout << "--- Starting Error Recovery Analysis (Indices 1 to " << N << ") ---" << std::endl;
    auto collect_start = std::chrono::high_resolution_clock::now();
    std::map<int, int> iter_counts;

    for (size_t idx = 1; idx <= N; ++idx) {
        uint32_t current_iter = 0;
        ComplexLD est_e0 = AnalyzeAndRecoverInitialError(idx, ct0, cryptoContext, keyPair, ptm, logQ, constant, current_iter, depth_estimate_attack, false, &main_loop_evalnoise_calls);
        estimated_e0_list[idx - 1] = est_e0;
        iter_counts[current_iter]++;
        min_last_iter_main = std::min(min_last_iter_main, current_iter);
        max_last_iter_main = std::max(max_last_iter_main, current_iter);

        // Real Values
        double dummy_inf;
        ComplexLD actual_e0_w;
        EvalNoiseBFVComplex(keyPair.secretKey, ct0, dummy_inf, actual_e0_w, idx);
        main_loop_evalnoise_calls++;

        double current_abs_error = std::abs(est_e0 - actual_e0_w);
        total_abs_error += current_abs_error;
        if (current_abs_error > max_abs_error) max_abs_error = current_abs_error;

        // Printing for Debugging
        if (idx % log_interval == 0) {
            std::cout << "[Progress " << idx/log_interval << "/8] Index " << std::setw(4) << idx 
                      << " | Iter: " << std::setw(2) << current_iter
                      << " | Abs Error: " << std::scientific << std::setprecision(5) << current_abs_error 
                      << " | Est: (" << est_e0.real() << ", " << est_e0.imag() << "i)" 
                      << std::defaultfloat << std::endl;
        }

        if (current_abs_error > 0.5) {
            std::cout << "!!! [Critical] Error > 0.5 at index " << idx << ". Stopping." << std::endl;
            break;
        }
    }
    auto collect_end = std::chrono::high_resolution_clock::now();
    
    // ---------------//
    // Error recovery //
    // ---------------//
    double l2_error = 0;
    double l2_residual_b_plus_as_minus_e = -1.0;
    double max_abs_residual_b_plus_as_minus_e = -1.0;
    double sk_coeff_match_ratio = -1.0;
    double time_recovery_e_sec = 0.0;
    double time_recovery_s_sec = 0.0;
    auto recovery_start = std::chrono::high_resolution_clock::now();
    if (max_abs_error < 0.5) {
        auto recovery_e_start = std::chrono::high_resolution_clock::now();
        std::vector<int64_t> e0_est_coeffs(N);
        for (size_t j = 0; j < N; ++j) {
            ComplexLD sum_val(0, 0);
            for (size_t k = 1; k <= N; ++k) {
                double theta = (2.0 * k - 1.0) * M_PI * j / static_cast<double>(N);
                ComplexLD inv_omega_pow(std::cos(theta), -std::sin(theta));
                sum_val += estimated_e0_list[k - 1] * inv_omega_pow;
            }
            e0_est_coeffs[j] = std::round(sum_val.real() / static_cast<double>(N));
        }

        // Error Comparison
        const vector<DCRTPoly>& cv = ct0->GetElements();
        DCRTPoly s = keyPair.secretKey->GetPrivateElement();
        DCRTPoly actual_e_poly = cv[0] + cv[1] * s;
        actual_e_poly.SetFormat(Format::COEFFICIENT);
        auto actual_e_poly_int = actual_e_poly.CRTInterpolate();
        auto Q_mod = actual_e_poly_int.GetModulus();
        auto Q_half = Q_mod >> 1;

        double l2_error_sq = 0;
        for (size_t j = 0; j < N; ++j) {
            auto v = actual_e_poly_int[j];
            int64_t actual_coeff = (v > Q_half) ? -static_cast<int64_t>((Q_mod - v).ConvertToInt()) : static_cast<int64_t>(v.ConvertToInt());
            double diff = (double)e0_est_coeffs[j] - (double)actual_coeff;
            l2_error_sq += diff * diff;
        }
        l2_error = std::sqrt(l2_error_sq);

        // Verify ciphertext relation after removing the recovered error estimate.
        // actual_e_poly already represents c0 + s*c1 in the active modulus, so the
        // residual is actual_e_poly - e_est modulo Q.
        auto e_est_poly_int = actual_e_poly_int;
        for (size_t j = 0; j < N; ++j) {
            int64_t c = e0_est_coeffs[j];
            if (c >= 0) {
                e_est_poly_int[j] = static_cast<uint64_t>(c);
            } else {
                e_est_poly_int[j] = Q_mod - static_cast<uint64_t>(-c);
            }
        }

        auto residual_poly_int = actual_e_poly_int - e_est_poly_int;

        double residual_l2_sq = 0.0;
        double residual_max = 0.0;
        for (size_t j = 0; j < N; ++j) {
            auto v = residual_poly_int[j];
            double centered = (v > Q_half)
                                       ? -static_cast<double>((Q_mod - v).ConvertToDouble())
                                       : static_cast<double>(v.ConvertToDouble());
            residual_l2_sq += centered * centered;
            residual_max = std::max(residual_max, std::abs(centered));
        }
        l2_residual_b_plus_as_minus_e = static_cast<double>(std::sqrt(residual_l2_sq));
        max_abs_residual_b_plus_as_minus_e = static_cast<double>(residual_max);
        auto recovery_e_end = std::chrono::high_resolution_clock::now();
        time_recovery_e_sec = std::chrono::duration<double>(recovery_e_end - recovery_e_start).count();

        // ------------//
        // SK recovery //
        // ------------//
        // Recover secret key in R_q directly: s_hat = (e_hat - c0) * c1^{-1} mod q.
        auto recovery_s_start = std::chrono::high_resolution_clock::now();
        DCRTPoly c0_for_sk = cv[0];
        DCRTPoly c1_for_sk = cv[1];
        c0_for_sk.SetFormat(Format::EVALUATION);
        c1_for_sk.SetFormat(Format::EVALUATION);

        DCRTPoly e_est_dcrt(cv[0].GetParams(), Format::COEFFICIENT, true);
        e_est_dcrt = e0_est_coeffs;
        e_est_dcrt.SetFormat(Format::EVALUATION);

        std::vector<int64_t> sk_est_coeffs(N, 0);
        std::vector<int64_t> sk_true_coeffs(N, 0);
        std::vector<int64_t> sk_est_clipped(N, 0);

        if (c1_for_sk.InverseExists()) {
            DCRTPoly c1_inv = c1_for_sk.MultiplicativeInverse();
            DCRTPoly s_est = (e_est_dcrt - c0_for_sk) * c1_inv;
            s_est.SetFormat(Format::COEFFICIENT);
            auto s_est_int = s_est.CRTInterpolate();
            auto s_est_mod = s_est_int.GetModulus();
            auto s_est_half = s_est_mod >> 1;
            for (size_t j = 0; j < N; ++j) {
                auto v = s_est_int[j];
                int64_t centered = (v > s_est_half) ? -static_cast<int64_t>((s_est_mod - v).ConvertToInt())
                                                    : static_cast<int64_t>(v.ConvertToInt());
                sk_est_coeffs[j] = centered;
            }

            DCRTPoly s_true = keyPair.secretKey->GetPrivateElement();
            size_t sizeQl_sk = cv[0].GetParams()->GetParams().size();
            size_t sizeQs_sk = s_true.GetParams()->GetParams().size();
            if (sizeQs_sk > sizeQl_sk) {
                s_true.DropLastElements(sizeQs_sk - sizeQl_sk);
            }
            s_true.SetFormat(Format::COEFFICIENT);
            auto s_true_int = s_true.CRTInterpolate();
            auto s_mod = s_true_int.GetModulus();
            auto s_mod_half = s_mod >> 1;

            size_t sk_match_cnt = 0;
            for (size_t j = 0; j < N; ++j) {
                auto v = s_true_int[j];
                int64_t centered = (v > s_mod_half) ? -static_cast<int64_t>((s_mod - v).ConvertToInt())
                                                    : static_cast<int64_t>(v.ConvertToInt());
                sk_true_coeffs[j] = centered;

                int64_t clip = 0;
                if (sk_est_coeffs[j] > 0) {
                    clip = 1;
                } else if (sk_est_coeffs[j] < 0) {
                    clip = -1;
                }
                sk_est_clipped[j] = clip;
                if (clip == centered) {
                    sk_match_cnt++;
                }
            }
            sk_coeff_match_ratio = static_cast<double>(sk_match_cnt) / static_cast<double>(N);

            const size_t debug_coeffs = std::min<size_t>(50, N);
            std::cout << "\n" << std::string(60, '=') << std::endl;
            std::cout << "DEBUG: FIRST " << debug_coeffs << " SK COEFFS (TRUE VS RECOVERED)" << std::endl;
            std::cout << std::string(60, '-') << std::endl;
            for (size_t j = 0; j < debug_coeffs; ++j) {
                std::cout << "j=" << std::setw(3) << j
                          << " | sk_true=" << std::setw(2) << sk_true_coeffs[j]
                          << " | sk_est_raw=" << std::setw(4) << sk_est_coeffs[j]
                          << " | sk_est_clip=" << std::setw(2) << sk_est_clipped[j]
                          << ((sk_true_coeffs[j] == sk_est_clipped[j]) ? " | OK" : " | XX")
                          << std::endl;
            }
            std::cout << "SK Coeff Match Ratio    : " << std::fixed << std::setprecision(6)
                      << sk_coeff_match_ratio << std::defaultfloat << std::endl;
            std::cout << std::string(60, '=') << std::endl;
        } else {
            std::cout << "[WARN] c1 is not invertible in R_q; skipped direct sk recovery." << std::endl;
        }
        auto recovery_s_end = std::chrono::high_resolution_clock::now();
        time_recovery_s_sec = std::chrono::duration<double>(recovery_s_end - recovery_s_start).count();
    } else {
        l2_error = -1.0;
    }
    auto recovery_end = std::chrono::high_resolution_clock::now();

    // --------------//
    // Final Summary //
    // --------------//
    std::chrono::duration<double> time_collect = collect_end - collect_start;
    std::chrono::duration<double> time_recovery = recovery_end - recovery_start;

    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "SUMMARY OF ANALYSIS" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Average Abs Error (e(w)): " << std::scientific << (total_abs_error / N) << std::endl;
    std::cout << "Maximum Abs Error (e(w)): " << max_abs_error << std::defaultfloat << std::endl;
    std::cout << "Final L2 Error (e)      : " << l2_error << std::endl;
    std::cout << "L2 Residual (c0-e+sc1)  : " << l2_residual_b_plus_as_minus_e << std::endl;
    std::cout << "Max Residual (c0-e+sc1) : " << max_abs_residual_b_plus_as_minus_e << std::endl;
    std::cout << "SK Match Ratio (coeff)  : " << sk_coeff_match_ratio << std::endl;
    std::cout << "EvalNoiseBFVComplex calls (main loop): " << main_loop_evalnoise_calls << " ~ " << static_cast<double>(main_loop_evalnoise_calls) / static_cast<double>(N) << " N" << std::endl;
    if (min_last_iter_main != std::numeric_limits<uint32_t>::max()) {
        std::cout << "last_iter range (main loop): [" << min_last_iter_main << ", " << max_last_iter_main << "]";
    } else {
        std::cout << "last_iter range (main loop): N/A";
    }
    if (min_last_iter_main != max_last_iter_main) {
        std::cout << " --> [";
        for (auto it = iter_counts.begin(); it != iter_counts.end(); ++it) {
            std::cout << it->first << ": " << it->second;
            if (std::next(it) != iter_counts.end()) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    } else {
        std::cout << std::endl;
    }
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Collecting Time (e(w))  : " << time_collect.count() << " s" << std::endl;
    std::cout << "Recovery Time (e)       : " << time_recovery.count() << " s ("
              << time_recovery_e_sec << " s for e + "
              << time_recovery_s_sec << " s for s)" << std::endl;
    std::cout << "Total Elapsed Time      : " << (time_collect.count() + time_recovery.count()) << " s" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    return 0;
}