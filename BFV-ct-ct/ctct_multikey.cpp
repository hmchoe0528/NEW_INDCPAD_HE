// Robustness check: repeat phases 2 and 3 across independent key pairs, so that
// sk, pk and all encryption randomness differ from run to run.
//
// The gain C has to be matched to the key: the amplification scales as k_s^K, and
// k_s varies from key to key (3.19 to 3.69 observed over ten keys), a swing of
// K*log2(3.69/3.19) = 9.2 bits against a readout window (q/2t, q/2) only log2(t)
// bits wide. A fixed C therefore saturates for some keys and undershoots for
// others. We rescale it from the measured k_s,
//     log2 C = log2 C_ref - K * log2(k_s / k_s_ref),
// which is the analogue, for an experiment that holds sk, of the binary search an
// adversary would run against the decryption oracle.
#include "bfv_ctct.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace lbcrypto; using namespace std;

int main(int argc, char** argv) {
    uint32_t KEYS   = (argc >  1) ? (uint32_t)atoi(argv[1])  : 10;
    uint32_t TRIALS = (argc >  2) ? (uint32_t)atoi(argv[2])  : 3;
    uint64_t ptm    = (argc >  3) ? (uint64_t)atoll(argv[3]) : 786433;
    uint32_t depth  = (argc >  4) ? (uint32_t)atoi(argv[4])  : 35;
    double   kappa  = (argc >  5) ? atof(argv[5])            : 4.5;
    uint64_t rhoRef = (argc >  6) ? (uint64_t)atoll(argv[6]) : 24800;
    double   ksRef  = (argc >  7) ? atof(argv[7])            : 3.36;
    uint32_t ringd  = (argc >  8) ? (uint32_t)atoi(argv[8])  : 0;      // 0 = let OpenFHE choose
    uint32_t M      = (argc >  9) ? (uint32_t)atoi(argv[9])  : 262144;
    uint32_t adapt  = (argc > 10) ? (uint32_t)atoi(argv[10]) : 0;   // 0 = fixed C (the adversary cannot see k_s)
    const uint32_t K23 = depth - 1;

    CCParams<CryptoContextBFVRNS> P;
    P.SetPlaintextModulus(ptm); P.SetMultiplicationTechnique(HPS);
    P.SetMultiplicativeDepth(depth); P.SetKeySwitchTechnique(HYBRID);
    P.SetScalingModSize(52); P.SetNumLargeDigits(depth >= 10 ? 5 : 2);
    if (ringd) { P.SetRingDim(ringd); P.SetSecurityLevel(SecurityLevel::HEStd_NotSet); }
    else       { P.SetSecurityLevel(SecurityLevel::HEStd_128_classic); }
    auto cc = GenCryptoContext(P);
    cc->Enable(PKE); cc->Enable(KEYSWITCH); cc->Enable(LEVELEDSHE);

    const double logCref = log2((double)rhoRef * (double)ptm);
    printf("t = %llu, depth = %u (K = %u), k_a = %.2f, M = 2^%.0f, %u keys x %u runs\n",
           (unsigned long long)ptm, depth, K23, kappa, log2((double)M), KEYS, TRIALS);
    printf("adaptive gain: log2 C = %.2f - %u*log2(k_s/%.3f)\n\n", logCref, K23, ksRef);
    printf("%4s %8s %7s %8s | %-19s | %-19s | %7s %7s %7s | %8s %8s | %8s %8s | %s\n",
           "key", "omega", "k_s", "log2 C", "max|y| D0 (runs)", "max|y| D1 (runs)",
           "maxD0", "minD1", "ratio", "rms D0", "rms D1", "lg|e0|D0", "D1", "note");
    fflush(stdout);

    uint32_t good = 0;
    uint32_t apart = 0;
    for (uint32_t k = 0; k < KEYS; ++k) {
        auto kp = cc->KeyGen(); cc->EvalMultKeysGen(kp.secretKey);
        const auto cp = static_pointer_cast<CryptoParametersBFVRNS>(kp.secretKey->GetCryptoParameters());
        usint N = cp->GetElementParams()->GetRingDimension();
        BigInteger Q = cp->GetElementParams()->GetModulus(), T(ptm);
        uint64_t r = (Q - (Q / T) * T).ConvertToInt();

        DCRTPoly s = kp.secretKey->GetPrivateElement(); s.SetFormat(Format::COEFFICIENT);
        auto sb = s.CRTInterpolate(); auto sQ = sb.GetModulus(), sQh = sQ >> 1;
        vector<double> sr(N);
        for (size_t j = 0; j < N; ++j) sr[j] = (sb[j] > sQh) ? -(sQ - sb[j]).ConvertToDouble() : sb[j].ConvertToDouble();
        auto sev = fft(sr); size_t mi = 0; double mv = -1;
        for (size_t i = 0; i < N / 2; ++i) if (abs(sev[i]) > mv) { mv = abs(sev[i]); mi = i; }
        const double ks = mv / sqrt(2.0 * N / 3.0);
        if (k == 0) printf("N = 2^%.0f, log2 q = %.1f, q/2t = 2^%.1f, (4Nt)^K = 2^%.1f\n\n",
                           log2((double)N), Log2BigInteger(Q),
                           Log2BigInteger(Q) - 1 - log2((double)ptm),
                           K23 * log2(4.0 * N * (double)ptm));

        // rescale the gain for this key, then round to a multiple of t so that t | C
        const double logC = adapt ? logCref - (double)K23 * log2(ks / ksRef) : logCref;
        uint64_t rho = (uint64_t)llround(pow(2.0, logC) / (double)ptm);
        if (rho < 1) rho = 1;
        const uint64_t C = rho * ptm;

        Plaintext pt_star = LargeSigPtxt(cc, mi, ptm, r);
        Plaintext pt_zero = cc->MakeCoefPackedPlaintext(vector<int64_t>(N, 0));

        vector<int64_t> mx[2]; vector<double> rms[2]; double e0w[2] = {0, 0};
        const vector<uint32_t> tgt = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900};
        vector<vector<int64_t>> coef[2];              // per run, the values at tgt
        const auto t_key = chrono::steady_clock::now();
        for (uint32_t t = 0; t < TRIALS; ++t)
            for (int g = 0; g < 2; ++g) {
                auto ct = EncryptLinComb(cc, kp, g ? pt_star : pt_zero, M, C);
                double a, b; EvalNoiseBFV(kp.secretKey, ct, a, b, mi); e0w[g] += b;
                for (uint32_t d = 1; d <= K23; ++d) {
                    Ciphertext<DCRTPoly> rr;
                    for (int n = 0; n < 50; ++n) { try { rr = reRandomizeApart(kp.secretKey, ct, ptm, mi, kappa, false); break; } catch (...) { if (n == 49) throw; } }
                    ct = cc->EvalMult(rr, rr);
                }
                Plaintext dec; cc->Decrypt(kp.secretKey, ct, &dec);
                int64_t m = 0; double ss = 0.0; size_t nc = 0;
                for (auto v : dec->GetCoefPackedValue()) {
                    if (llabs(v) > m) m = llabs(v);
                    ss += (double)v * (double)v; ++nc;
                }
                mx[g].push_back(m);
                rms[g].push_back(sqrt(ss / (double)nc));   // = A/sqrt(2) for a sinusoid of amplitude A
                const auto& dv = dec->GetCoefPackedValue();
                vector<int64_t> row; for (auto idx : tgt) row.push_back(idx < dv.size() ? dv[idx] : 0);
                coef[g].push_back(row);
            }
        // report every run, not just the extreme over them
        auto fmt = [](const vector<int64_t>& v) {
            string s; for (size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + to_string(v[i]); return s;
        };
        const int64_t w0 = *max_element(mx[0].begin(), mx[0].end());
        const int64_t w1 = *min_element(mx[1].begin(), mx[1].end());
        auto avg = [](const vector<double>& v){ double t=0; for(double x:v) t+=x; return t/v.size(); };
        const int64_t half = (int64_t)(ptm / 2);
        // No pass mark is imposed. What matters is whether the worst run under
        // D1 still stands clear of the best under D0, so both bounds and their
        // ratio are printed and the reader judges. The two notes that remain are
        // not thresholds but the failure modes the construction has: the noise
        // never reaching the readout window, and its wrapping past q/2.
        const bool disjoint = (w1 > w0);
        apart += disjoint;
        good += disjoint;
        const char* note = (w1 == 0)         ? "undershoot, below q/2t"
                         : (w0 > half / 2)   ? "saturated, wrapped past q/2"
                         : disjoint          ? ""
                         :                     "ranges overlap";
        printf("%4u %8zu %7.3f %8.2f | %-19s | %-19s | %7lld %7lld %7.1f | %8.1f %8.1f | %8.2f %8.2f | %s\n",
               k, mi, ks, logC, fmt(mx[0]).c_str(), fmt(mx[1]).c_str(),
               (long long)w0, (long long)w1,
               (double)w1 / (double)max<int64_t>(w0, 1),
               avg(rms[0]), avg(rms[1]), e0w[0] / TRIALS, e0w[1] / TRIALS, note);
        // the values themselves, so that per-key distributions can be redrawn later
        const double secs = chrono::duration<double>(chrono::steady_clock::now() - t_key).count();
        printf("  key %u: %u runs at indices 0,100,...,900 (D0|D1), %.1f s\n", k, TRIALS, secs);
        for (uint32_t t = 0; t < TRIALS; ++t) {
            printf("  %4u | ", t + 1);
            for (size_t i2 = 0; i2 < tgt.size(); ++i2)
                printf("%7lld | %7lld || ", (long long)coef[0][t][i2], (long long)coef[1][t][i2]);
            printf("\n");
        }
        fflush(stdout);
    }
    printf("\n%u/%u key pairs separate: the worst run under D1 exceeds the best"
           " under D0.\nThe gain C was %s.\n", apart, KEYS,
           adapt ? "rescaled from the measured k_s" : "fixed in advance, as an adversary must");
    return 0;
}
