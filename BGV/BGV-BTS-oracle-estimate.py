"""
Theoretical estimate of the number of decryption oracle calls for the
KR-D attack on HElib thin bootstrapping (Section 3).

Model (Section 3, "Key recovery from inequality hints"):
  - sigma^2 = (p^{2e'} + p^{2(e+r)}/q~^2 + 1) / 12, q~ = p^e + 1.
  - c_j ~ D^(l): sum of 2^l i.i.d. copies of
        U[-p^{e'}/2, p^{e'}/2) + (p^{e+r}/q~) U[-1/2, 1/2) + (1/q~) U[-q~/2, q~/2).
  - Y_j ~ N(0, (1 + 2^t h phi(m)/m - 2^l) sigma^2) if s_j = +-1.
  - p_1^(l)(x) = p_0^(l)(x) g(x) / Z, g(x) = Pr[|x + Y_j| > p^e/2], Z = E_{p_0}[g].
  - Per-slot failure probability (HS21, Eq. (5)), as in Table 2:
        p_fail = erfc( (p^e/2) / sqrt(2 * 2^t h phi(m)/m * sigma^2) ).

Number of failing slots N to classify all m indices, treating the samples
of c_j from distinct failing slots as independent:
  - N_KL = ln(m) / D_KL(p_1 || p_0): required (Chernoff-Stein, error 1/m).
    This is the estimate given in the paper.
  - N_C = ln(m) / C(p_1, p_0), C = -min_lam ln int p_1^lam p_0^(1-lam):
    sufficient. The likelihood-ratio test has both error probabilities
    <= exp(-N C), so at most one index is misclassified in expectation.
#O_dec = N / (n_slot * p_fail). The paper uses l = 1, since the unique
l = 0 index is always kept as a candidate.

Usage: python3 BGV-BTS-oracle-estimate.py [--dx 0.0078125]
"""
import argparse
import math

import numpy as np

# name: (m, phi(m), n_slot, h, e)
PARAMS = {
    "A":  (21845, 16384, 1024, 120, 12),
    "B":  (28679, 23040,  960, 120, 12),
    "A*": (21845, 16384, 1024,  32,  8),
    "B*": (28679, 23040,  960,  32,  8),
}
P, EP, R, T = 2, 4, 1, 3  # p, e', r, t (number of prime factors of m)


def irwin_hall(x, n, w):
    """Density of the sum of n i.i.d. U[-w/2, w/2)."""
    u = x / w + n / 2
    if n == 1:
        return ((u >= 0) & (u < 1)) / w
    out = np.zeros_like(u)
    for k in range(n + 1):
        out += (-1) ** k * math.comb(n, k) * np.clip(u - k, 0, None) ** (n - 1)
    out /= math.factorial(n - 1) * w
    out[(u < 0) | (u > n)] = 0.0
    return out


def p0_density(n, e, dx):
    """Density of D^(l) for n = 2^l terms, on a grid of step dx."""
    qt = P ** e + 1
    widths = [P ** EP, P ** (e + R) / qt, 1.0]
    dens = None
    for w in widths:
        half = n * w / 2
        x = np.arange(-math.ceil(half / dx), math.ceil(half / dx) + 1) * dx
        d = irwin_hall(x, n, w)
        d /= d.sum() * dx
        dens = d if dens is None else np.convolve(dens, d) * dx
    k = (len(dens) - 1) // 2
    return np.arange(-k, k + 1) * dx, dens


def log_tail(x, thr, s):
    """ln Pr[|x + N(0, s^2)| > thr]."""
    a = np.array([math.erfc((thr - v) / (s * math.sqrt(2))) for v in x])
    b = np.array([math.erfc((thr + v) / (s * math.sqrt(2))) for v in x])
    return np.log(0.5 * (a + b))


def divergences(x, p0, thr, s, dx):
    """D_KL(p_1 || p_0) and Chernoff information C(p_1, p_0)."""
    lg = log_tail(x, thr, s)
    lg -= lg.max()
    lz = math.log(np.sum(p0 * np.exp(lg)) * dx)
    llr = lg - lz                                 # ln(p_1 / p_0)
    p1 = p0 * np.exp(llr)
    kl = float(np.sum(p1 * llr) * dx)
    lams = np.linspace(0.01, 0.99, 99)
    ch = max(-math.log(np.sum(p0 * np.exp(lam * llr)) * dx) for lam in lams)
    return kl, ch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dx", type=float, default=1 / 128)
    dx = ap.parse_args().dx

    for name, (m, phim, nslot, h, e) in PARAMS.items():
        qt = P ** e + 1
        sig2 = (P ** (2 * EP) + P ** (2 * (e + R)) / qt ** 2 + 1) / 12
        thr = P ** e / 2
        sb2 = 2 ** T * h * phim / m * sig2
        pfail = math.erfc(thr / math.sqrt(2 * sb2))
        print(f"{name}: m={m}, n_slot={nslot}, h={h}, e={e}, "
              f"log2(p_fail)={math.log2(pfail):.3f}")
        print(f"  {'l':>2} {'D_KL':>10} {'C':>10} "
              f"{'N_KL':>9} {'log2 #O_KL':>10} {'N_C':>9} {'log2 #O_C':>10}")
        for l in range(T + 1):
            x, p0 = p0_density(2 ** l, e, dx)
            s = math.sqrt((1 + 2 ** T * h * phim / m - 2 ** l) * sig2)
            kl, ch = divergences(x, p0, thr, s, dx)
            n_kl, n_ch = math.log(m) / kl, math.log(m) / ch
            o_kl = n_kl / (nslot * pfail)
            o_ch = math.ceil(n_ch) / (nslot * pfail)
            print(f"  {l:>2} {kl:10.4g} {ch:10.4g} {n_kl:9.4g} "
                  f"{math.log2(o_kl):10.2f} {n_ch:9.4g} {math.log2(o_ch):10.2f}"
                  + (f"   #O_KL={o_kl:.0f}, #O_C={o_ch:.0f}" if e == 8 else ""))


if __name__ == "__main__":
    main()
