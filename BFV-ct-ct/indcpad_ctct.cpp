//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=================================================================================
#include "bfv_ctct.hpp"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <complex>
#include <chrono>

using namespace lbcrypto;
using namespace std;

int main(int argc, char** argv) {
    CCParams<CryptoContextBFVRNS> parameters;
    uint64_t ptm   = (argc > 1) ? (uint64_t)atoll(argv[1]) : 1153;
    uint32_t depth = (argc > 2) ? (uint32_t)atoi(argv[2])  : 45;
    double kappa   = (argc > 3) ? atof(argv[3])            : 8.15;
    usint dcrtBits = 52;

    // phase 2/3 parameters: M challenge queries per sample, total scale C = rho*t
    uint32_t M      = (argc > 4) ? (uint32_t)atoi(argv[4])  : 262144;
    uint64_t rho    = (argc > 5) ? (uint64_t)atoll(argv[5]) : 59600370;
    uint32_t trials = (argc > 6) ? (uint32_t)atoi(argv[6])  : 10;
    uint32_t K23    = depth - 1;
    double   ksMin  = (argc > 7) ? atof(argv[7])            : 0.0;   // regenerate the key until k_s >= ksMin
    double   ksMax  = (argc > 8) ? atof(argv[8])            : 1e9;  // ... and k_s <= ksMax, so that e_K stays under q/2
    double   msFac  = (argc > 9) ? atof(argv[9])            : 4.0;  // delta_MS/sqrt(N): 4 since OpenFHE v1.4.2.4, 2 before

    parameters.SetPlaintextModulus(ptm);
    parameters.SetMultiplicationTechnique(HPS);
    parameters.SetMultiplicativeDepth(depth);
    parameters.SetKeySwitchTechnique(HYBRID);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetNumLargeDigits(5);

    std::cout << "depth: " << depth << "  log(t): " << log2(ptm) << "\n";
    std::cout << "Parameter generation start... ";
    parameters.SetSecurityLevel(SecurityLevel::HEStd_128_classic);
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);

    // The figure needs one key on which the attack lands inside the readout window;
    // the distribution of k_s over keys is the business of ctct_multikey, not of this
    // run, so we simply draw again when k_s comes out too small.
    KeyPair<DCRTPoly> keyPair;
    size_t  max_idx = 0;
    double  max_val = -1.0, ks_val = 0.0;
    uint32_t draws  = 0;
    do {
        ++draws;
        keyPair = cryptoContext->KeyGen();
        {
            DCRTPoly s0 = keyPair.secretKey->GetPrivateElement();
            s0.SetFormat(Format::COEFFICIENT);
            auto sb0 = s0.CRTInterpolate();
            auto q0 = sb0.GetModulus(), qh0 = q0 >> 1;
            const size_t n0 = sb0.GetLength();
            std::vector<double> r0(n0);
            for (size_t j = 0; j < n0; ++j)
                r0[j] = (sb0[j] > qh0) ? -(q0 - sb0[j]).ConvertToDouble() : sb0[j].ConvertToDouble();
            auto ev0 = fft(r0);
            max_val = -1.0;
            for (size_t i = 0; i < n0 / 2; ++i)
                if (std::abs(ev0[i]) > max_val) { max_val = std::abs(ev0[i]); max_idx = i; }
            ks_val = max_val / std::sqrt(2.0 * n0 / 3.0);
        }
        if (ksMin > 0.0)
            std::cout << "  keygen draw " << draws << ": k_s = " << ks_val
                      << ((ks_val >= ksMin && ks_val <= ksMax) ? "  [accepted]" : "  [redrawn]")
                      << "\n" << std::flush;
    } while (ks_val < ksMin || ks_val > ksMax);
    if (ksMin > 0.0) {
        // Selecting the key is not something the adversary can do: k_s is fixed by sk and
        // is not observable. The redrawing here only serves to obtain a figure on a key
        // that lands inside the readout window; its price is that the attack succeeds on
        // a fraction P of the keys, and its cost has to be divided by that fraction.
        const double nn = 0.5 * (double)cryptoContext->GetRingDimension();
        auto F = [&](double k) { return 1.0 - std::exp(-nn * std::exp(-k * k)); };
        const double P = F(ksMin) - (ksMax < 1e8 ? F(ksMax) : 0.0);
        std::cout << "\n  The key was redrawn until k_s fell in [" << ksMin << ", " << ksMax
                  << "], which took " << draws << " draws.\n"
                  << "  An adversary cannot do this: k_s is fixed by sk and is not observable.\n"
                  << "  A random key meets the band with probability P = " << P
                  << ", so the cost of the\n  attack is to be multiplied by 1/P = " << 1.0 / P
                  << " (equivalently, its advantage multiplied by P).\n"
                  << "  See ctct_multikey for the outcome over keys that are not selected.\n\n"
                  << std::flush;
    }
    cryptoContext->EvalMultKeysGen(keyPair.secretKey);

    const auto cryptoParams = std::static_pointer_cast<CryptoParametersBFVRNS>(keyPair.secretKey->GetCryptoParameters());
    std::cout << "End!\n";

    usint  N      = cryptoParams->GetElementParams()->GetRingDimension();
    size_t sizeQ  = cryptoParams->GetElementParams()->GetParams().size();
    size_t sizeQP = cryptoParams->GetParamsQP()->GetParams().size();

    double logQ = 0;
    for (usint i = 0; i < sizeQ; i++)
        logQ += log2(cryptoParams->GetElementParams()->GetParams()[i]->GetModulus().ConvertToInt());

    double logP = 0;
    for (size_t j = sizeQ; j < sizeQP; j++)
        logP += log2(cryptoParams->GetParamsQP()->GetParams()[j]->GetModulus().ConvertToInt());

    std::cout << "log(Q): " << logQ << "  (sizeQ=" << sizeQ << ", sizeQP=" << sizeQP << ")\n";
    std::cout << "log(P): " << logP << "  log(N): " << log2(N) << "\n";
    std::cout << "depth: " << depth << "\n";

    DCRTPoly s = keyPair.secretKey->GetPrivateElement();
    s.SetFormat(Format::COEFFICIENT);
    auto s_big = s.CRTInterpolate();
    auto sQ    = s_big.GetModulus();
    auto sQh   = sQ >> 1;
    std::vector<double> s_raw(N);
    for (size_t j = 0; j < N; ++j) {
        auto v   = s_big[j];
        s_raw[j] = (v > sQh) ? -static_cast<double>((sQ - v).ConvertToDouble())
                              :  static_cast<double>(v.ConvertToDouble());
    }
    std::vector<complex<double>> s_evals = fft(s_raw);

    max_val = -1.0;
    for (size_t i = 0; i < N / 2; ++i) {
        if (std::abs(s_evals[i]) > max_val) {
            max_val = std::abs(s_evals[i]);
            max_idx = i;
        }
    }

    uint32_t Mo = 2 * N;
    uint32_t g = 1;
    for (size_t i = 0; i < max_idx; i++) {
        g = (g * 5) % Mo;
    }
    size_t best_k = (g + 1) / 2;

    double kappaDenom = std::sqrt(2.0 * N / 3.0);
    ks_val = max_val / kappaDenom;

    std::cout << "\n[Target Identification]\n";
    std::cout << "Selected max_idx: " << max_idx << " with |s(omega)|: " << max_val << " (ks = " << ks_val << ")\n";
    std::cout << "Inferred k: " << best_k << " (exponent g: " << g << ")\n";

    // ---- the challenge messages of phases 2 and 3 ----------------------------
    // mu_0 = 0 and mu_1 = pt*, the plaintext whose signature nu(pt*) has the
    // largest possible evaluation at omega.  See LargeSigPtxt in bfv_ctct.hpp.
    BigInteger Q = cryptoParams->GetElementParams()->GetModulus(), T(ptm);
    uint64_t r = (Q - (Q / T) * T).ConvertToInt();
    std::vector<int64_t> sig;
    Plaintext pt_star = LargeSigPtxt(cryptoContext, max_idx, ptm, r, &sig);
    Plaintext pt_zero = cryptoContext->MakeCoefPackedPlaintext(std::vector<int64_t>(N, 0));
    const uint64_t C = rho * ptm;

    std::cout << "\n[Challenge Message]\n";
    std::cout << "r = Q mod t: " << r << "\n";
    std::cout << "log2 |nu(pt*)(omega)|          : " << log2(std::abs(evalIntPoly(sig, N, max_idx))) << "\n";
    std::cout << "log2 |rho * nu(pt*)(omega)| : " << log2(std::abs(evalIntPoly(sig, N, max_idx))) + log2((double)rho) << "\n";
    std::cout << "M = " << M << ", rho = " << rho << ", C = rho*t = 2^" << log2((double)C) << "\n";
    std::cout << "----------------------------------------------------------\n";

    const double decBound   = logQ - 1.0 - log2(ptm);
    const double wrapThresh = logQ - 1.1;
    const uint32_t K        = depth;

    const double sigma_val = 3.2, alpha_val = 36.0;
    const double Berr_th   = sigma_val * std::sqrt(alpha_val);
    const double Bkey_th   = 1.0;
    const double delta_th  = 2.0 * std::sqrt(static_cast<double>(N));
    const double deltaMS_th= msFac * std::sqrt(static_cast<double>(N));  // 4 since OpenFHE v1.4.2.4, 2 before
    const double C1_th     = delta_th * deltaMS_th * static_cast<double>(ptm) * Bkey_th;
    const double Vnorm_th  = Berr_th * (1.0 + 2.0 * delta_th * Bkey_th);
    const double logC1_th  = std::log2(C1_th);

    const double tail = max_val / kappaDenom;

    std::cout << "\nStarting iterative multiplications for idx=" << max_idx << " (tail=" << tail << ", global kappa=" << kappa << ")\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  " << std::setw(4) << "d"
              << std::setw(18) << "log2||e||_inf"
              << std::setw(16) << "log2|e(omega)|"
              << std::setw(16) << "log2(B_d)_th"
              << "  status\n"
              << "  " << std::string(64, '-') << "\n";

    double logB_th = std::log2(Vnorm_th);
    auto ciphertext = EncryptZeroWithTrace(keyPair.publicKey).ciphertext;

    string final_status = "SUC";

    double noise_inf, noise_w;
    EvalNoiseBFV(keyPair.secretKey, ciphertext, noise_inf, noise_w, max_idx);
    std::cout << "  " << std::setw(4) << 0
              << std::setw(18) << noise_inf
              << std::setw(16) << noise_w
              << std::setw(16) << logB_th << "\n";

    for (uint32_t d = 1; d <= K; ++d) {
        auto rr_ct = reRandomizeApart(keyPair.secretKey, ciphertext, ptm, max_idx, kappa);
        ciphertext = cryptoContext->EvalMult(rr_ct, rr_ct);

        const double numTowers_d = static_cast<double>(sizeQ - (d - 1));
        const double noiseKS_d   = numTowers_d * (delta_th * Berr_th + deltaMS_th * Bkey_th + 1.0);
        const double C2_d        = delta_th * deltaMS_th * Bkey_th * Bkey_th / 2.0 + noiseKS_d;
        const double logC2_d     = std::log2(C2_d);
        const double a = logC1_th + logB_th, b = logC2_d;
        logB_th = (a >= b) ? a + std::log2(1.0 + std::pow(2.0, b - a))
                            : b + std::log2(1.0 + std::pow(2.0, a - b));

        EvalNoiseBFV(keyPair.secretKey, ciphertext, noise_inf, noise_w, max_idx);
        const char* st = (noise_inf < decBound)   ? "SUC"
                       : (noise_inf < wrapThresh) ? "FAIL" : "WRAP";
        final_status = st;

        std::cout << "  " << std::setw(4) << d
                  << std::setw(18) << noise_inf
                  << std::setw(16) << noise_w
                  << std::setw(16) << logB_th
                  << "  [" << st << "]\n";
    }

    {   // phases 2 and 3 start from a much larger noise than phase 1, so they are run regardless
        // ---- phases 2 and 3 ---------------------------------------------------
        // For gamma in {0,1}: query the encryption oracle M times on mu_gamma,
        // combine with integer coefficients summing to C (a multiple of t), so
        // that the tracked plaintext is 0 and the noise is
        //     sum_i c_i e_i  -  rho * nu(mu_gamma) .
        // Then run the same K-fold amplification as phase 1 and decrypt.  The
        // decryption fails for both gamma, and the distinguisher is the size of
        // the returned coefficients.
        const double advertised = log2((double)C) + log2(Vnorm_th) + K23 * log2(2.0 * msFac * N * (double)ptm);
        std::cout << "\n[Phase 2/3]\n";
        std::cout << "certified bound C*V_N*(" << (uint32_t)(2.0 * msFac) << "Nt)^" << K23 << " = 2^" << advertised
                  << "  vs  q/2t = 2^" << decBound
                  << "  -> " << (advertised < decBound ? "within budget" : "OVER BUDGET")
                  << " (margin " << decBound - advertised << " bits)\n";

        std::vector<uint32_t> target_indices = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900};

        std::cout << "\nAmplified linear-combination Trials (" << trials << " trials, 10 Indices)\n";
        std::cout << "------------------------------------------------------------------------------------------\n";
        std::cout << "Trial | ";
        for (auto idx : target_indices) std::cout << "idx " << std::setw(5) << idx << " (ct0'|ct1') | ";
        std::cout << "\n------------------------------------------------------------------------------------------\n";

        for (uint32_t i = 1; i <= trials; ++i) {
            const auto t_tr = std::chrono::steady_clock::now();
            Plaintext dec[2];
            double e_inf[2], e_w[2], e0_inf[2], e0_w[2];
            for (int gamma = 0; gamma < 2; ++gamma) {
                auto ct = EncryptLinComb(cryptoContext, keyPair, gamma ? pt_star : pt_zero, M, C);
                EvalNoiseBFV(keyPair.secretKey, ct, e0_inf[gamma], e0_w[gamma], max_idx);
                for (uint32_t d = 1; d <= K23; ++d) {
                    auto rr_ct = reRandomizeApart(keyPair.secretKey, ct, ptm, max_idx, kappa);
                    ct = cryptoContext->EvalMult(rr_ct, rr_ct);
                }
                EvalNoiseBFV(keyPair.secretKey, ct, e_inf[gamma], e_w[gamma], max_idx);
                cryptoContext->Decrypt(keyPair.secretKey, ct, &dec[gamma]);
            }
            const auto& v0 = dec[0]->GetCoefPackedValue();
            const auto& v1 = dec[1]->GetCoefPackedValue();

            std::cout << std::setw(5) << i << " | ";
            for (auto idx : target_indices) {
                int64_t c0 = (idx < v0.size()) ? v0[idx] : 0;
                int64_t c1 = (idx < v1.size()) ? v1[idx] : 0;
                std::cout << std::setw(6) << c0 << " | " << std::setw(6) << c1 << " || ";
            }
            std::cout << "   before: log2||e_0||_inf " << e0_inf[0] << " / " << e0_inf[1]
                      << ", log2|e_0(w)| " << e0_w[0] << " / " << e0_w[1]
                      << " | after: log2||e_K||_inf " << e_inf[0] << " / " << e_inf[1]
                      << " | " << std::fixed << std::setprecision(1)
                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_tr).count()
                      << " s\n" << std::defaultfloat << std::setprecision(6);
            std::cout.flush();
        }
    }

    return 0;
}
