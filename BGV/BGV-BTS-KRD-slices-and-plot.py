#!/usr/bin/env python3
"""
Fast analysis of coeff.txt with multiple ciphertext slices.
Tests with first N ciphertexts: 50, 100, 150, 200, ..., 500
"""
import argparse
import math
import numpy as np
from pathlib import Path
from typing import Tuple, List, Dict
from scipy.stats import gaussian_kde
from scipy.signal import fftconvolve
from collections import defaultdict

# ============================================================================
# Parameters
# ============================================================================
### Param A* and B*. Selected with --m on the command line; B* is the default,
### as it was before the switch became an option.
PARAM_SETS = {
    21845: {"phi_m": 16384, "nslots": 1024},   # A*
    28679: {"phi_m": 23040, "nslots": 960},    # B*
}
m = 28679
phi_m = PARAM_SETS[m]["phi_m"]
nslots = PARAM_SETS[m]["nslots"]
###
p = 2
e = 8
hwt = 32
ePrime = 4
r = 1

# Data slicing points (number of ciphertexts). --slices overrides the default,
# which is the range the paper reports; a longer collection can be swept further.
slice_points = list(range(100, 220, 10))

# Selection parameters
n1 = 1
n2 = 32
n4 = 32
n8 = 32

# ============================================================================
# Utility Functions
# ============================================================================

def parse_secret_keys_file(path: Path, section: str = "bts") -> Tuple[List, List]:
    """
    Parse secret_keys.txt
    """
    if not path.exists():
        return [], []

    text = path.read_text()
    if section == "bts":
        marker = "bts secrekey generation start"
        start = text.find(marker)
        if start == -1:
            return [], []
        # HElib may draw the sparse bootstrapping key more than once before one
        # meets its bound, and every attempt is logged. The last is the one that
        # was used; reading the first gives a key that never encrypted anything,
        # and the recovery then matches none of it.
        blocks = [b for b in text[start:].split("Done!") if "sk +1" in b]
        if not blocks:
            return [], []
        target_block = blocks[-1]
    else:
        target_block = text.split("Done!")[0]

    pos, neg = [], []
    for line in target_block.splitlines():
        line = line.strip()
        if line.startswith("sk +1 :"):
            pos = [int(x) for x in line.split(":", 1)[1].replace(",", " ").split()]
        elif line.startswith("sk -1 :"):
            neg = [int(x) for x in line.split(":", 1)[1].replace(",", " ").split()]
    return pos, neg


def build_irwin_hall_dict(matrix_irwin_hall: np.ndarray) -> List[Dict]:
    """
    Build irwin_hall_dictionary from irwin_hall_*.txt matrix
    """
    n_rows, n_cols = matrix_irwin_hall.shape
    result = [None] * n_cols  # Pre-allocate list

    for j in range(n_cols):
        col = matrix_irwin_hall[:, j]
        num = np.max(col) if col.size > 0 else 1
        indices = np.where(col == num)[0].tolist()
        result[j] = {num: indices}

    return result


def count_ciphertexts_until(list_fail_slots: np.ndarray, target_count: int) -> int:
    """
    Find the row index for the target number of ciphertexts
    """
    if len(list_fail_slots) == 0:
        return 0

    # Vectorized approach: find where differences are negative
    diffs = np.diff(list_fail_slots)
    ctxt_boundaries = np.where(diffs < 0)[0] + 1

    if len(ctxt_boundaries) >= target_count:
        return int(ctxt_boundaries[target_count - 1])

    return len(list_fail_slots)  # Return all if target not reached


def compute_variances_optimized(ctxt_lists: np.ndarray,
                               irwin_hall_dict: List[Dict],
                               m: int) -> Tuple[List, List, List, List]:
    """
    Compute variances for all indices, grouped by t
    """
    list_fail_slots = ctxt_lists[:, 0].astype(np.int32)  # Use int32 for memory efficiency
    list_ctxt_a_pwfl = ctxt_lists[:, 2:]

    # Pre-allocate lists
    var_sum1 = []
    var_sum2 = []
    var_sum4 = []
    var_sum8 = []

    slots_array = list_fail_slots

    for idx in range(m):
        t_val, indices = list(irwin_hall_dict[idx].items())[0]

        # Use boolean mask for faster filtering
        mask = np.isin(slots_array, indices)
        if not np.any(mask):
            continue

        var = np.var(list_ctxt_a_pwfl[mask, idx])

        if t_val == 1:
            var_sum1.append((var, idx))
        elif t_val == 2:
            var_sum2.append((var, idx))
        elif t_val == 4:
            var_sum4.append((var, idx))
        elif t_val == 8:
            var_sum8.append((var, idx))

    # Sort each list
    var_sum1.sort(key=lambda x: x[0])
    var_sum2.sort(key=lambda x: x[0])
    var_sum4.sort(key=lambda x: x[0])
    var_sum8.sort(key=lambda x: x[0])

    return var_sum1, var_sum2, var_sum4, var_sum8


def select_candidates(var_sum1: List, var_sum2: List,
                     var_sum4: List, var_sum8: List,
                     n1: int, n2: int, n4: int, n8: int) -> List[int]:
    """
    Select top N candidates from each t group.
    """
    selected1 = var_sum1[-n1:] if n1 > 0 else []
    selected2 = var_sum2[-n2:] if n2 > 0 else []
    selected4 = var_sum4[-n4:] if n4 > 0 else []
    selected8 = var_sum8[-n8:] if n8 > 0 else []

    estimated_nonzero = sorted(set(
        idx
        for var_list in [selected1, selected2, selected4, selected8]
        for _, idx in var_list
    ))

    return estimated_nonzero


def evaluate(estimated_nonzero: List[int],
            sk_pos: List[int], sk_neg: List[int]) -> Dict:
    """
    Evaluate against reference secret.
    """
    reference_available = len(sk_pos) > 0 or len(sk_neg) > 0
    input_sk_nonzero = sorted(set(sk_pos) | set(sk_neg))

    result = {
        'reference_available': reference_available,
        'total_candidates': len(estimated_nonzero),
        'total_true': len(input_sk_nonzero),
        'correct': 0,
        'recall': 0.0,
        'precision': 0.0,
    }

    if reference_available:
        correct = len(set(estimated_nonzero) & set(input_sk_nonzero))
        result['correct'] = correct
        result['recall'] = correct / len(input_sk_nonzero) if input_sk_nonzero else 0
        result['precision'] = correct / len(estimated_nonzero) if estimated_nonzero else 0

    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze coeff.txt slices and optionally plot full ciphertext variance distributions."
    )
    parser.add_argument(
        "--plot-vars",
        action="store_true",
        help="Plot distributions of variances computed on full ciphertext data.",
    )
    parser.add_argument(
        "--plot-coeff-dist",
        action="store_true",
        help=(
            "Plot cj distributions for each (sj, ell) combination. "
            "Reproduces cells 5-6 from BGV-BTS-KRD-coeff-plot.ipynb."
        ),
    )
    parser.add_argument(
        "--slices",
        metavar="START:STOP:STEP",
        help="Ciphertext counts to sweep, as start:stop:step. "
             "Default 100:220:10, the range the paper reports.",
    )
    parser.add_argument(
        "--m",
        type=int,
        choices=sorted(PARAM_SETS),
        default=28679,
        help="Cyclotomic index: 21845 for parameter set A*, 28679 for B* (default).",
    )
    return parser.parse_args()


def plot_variance_distributions(
    var_sum1: List, var_sum2: List, var_sum4: List, var_sum8: List,
    actual_sk_indices: List[int],
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required to plot variance distributions. "
            "Install it with pip install matplotlib"
        ) from exc

    groups = [
        (var_sum1, "$\\ell$=0"),
        (var_sum2, "$\\ell$=1"),
        (var_sum4, "$\\ell$=2"),
        (var_sum8, "$\\ell$=3"),
    ]
    actual_set = set(actual_sk_indices)

    plt.figure(figsize=(14, 10))
    for idx, (var_list, label) in enumerate(groups, start=1):
        plt.subplot(2, 2, idx)
        values = [value for value, _ in var_list]
        if values:
            plt.hist(
                values,
                bins=80,
                color="#4c72b0",
                alpha=0.75,
                edgecolor="black",
            )

            ymax = plt.gca().get_ylim()[1]
            marker_y = ymax * 0.04
            highlight_values = [value for value, j in var_list if j in actual_set]
            if highlight_values:
                plt.scatter(
                    highlight_values,
                    [marker_y] * len(highlight_values),
                    color="red",
                    marker=r"$\downarrow$",
                    s=150,
                    linewidths=0.3,
                    zorder=5,
                    label="Actual SK Index",
                )
                plt.legend(fontsize=9)
        else:
            plt.text(0.5, 0.5, "No values", ha="center", va="center", fontsize=12)

        plt.title(f"\nVariance distribution ({label})")
        if idx == 1 or idx == 3:
            plt.ylabel("Count")
        if idx == 3 or idx == 4:
            plt.xlabel("Variance")
        plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


# ============================================================================
# Coefficient Distribution Functions (cells 5-6 from BGV-BTS-KRD-coeff-plot.ipynb)
# ============================================================================

def _irwin_hall_density(x, t, a):
    ta = t * a
    x_val = (x + ta) / (2 * a)
    if x_val < 0 or x_val > t:
        return 0.0
    k = int(np.floor(x_val))
    s = 0.0
    for j in range(k + 1):
        s += (-1)**j * math.comb(t, j) * (x_val - j)**(t - 1)
    return s / (2 * a)**t / math.factorial(t - 1)


def _symmetric_mod_m_pmf(modulus):
    n = modulus // 2
    if modulus % 2 == 1:
        support = np.arange(-n, n + 1)
        probs = np.ones_like(support, dtype=float) / modulus
    else:
        support = np.arange(-n, n + 1)
        probs = np.full_like(support, 1 / modulus, dtype=float)
        probs[0] = 1 / (2 * modulus)
        probs[-1] = 1 / (2 * modulus)
    return support, probs


def _irwin_hall_discrete_pmf(t, modulus):
    base, prob = _symmetric_mod_m_pmf(modulus)
    hist = prob.copy()
    for _ in range(t - 1):
        hist = np.convolve(hist, prob)
    support = np.arange(t * base[0], t * base[-1] + 1)
    return support, hist


def _expand_pmf_to_modQ_grid(pmf_support, pmf_values, numer):
    N = (len(pmf_support) - 1) // 2
    scaled_support = (pmf_support * numer).astype(int)
    grid_range = np.arange(-numer * N, numer * N + 1)
    padded_pmf = np.zeros(len(grid_range))
    for i, val in enumerate(scaled_support):
        index = val + numer * N
        padded_pmf[index] = pmf_values[i]
    return grid_range, padded_pmf


def _compute_stats(pdf, domain):
    mean = np.trapezoid(domain * pdf, domain)
    var = np.trapezoid((domain - mean) ** 2 * pdf, domain)
    return mean, var


def _compute_kl(p_dist, q_dist, domain, epsilon=1e-20):
    p_dist = np.clip(p_dist, epsilon, None)
    q_dist = np.clip(q_dist, epsilon, None)
    return np.trapezoid(p_dist * np.log(p_dist / q_dist), domain)


def _compute_theoretical_pdfs() -> Tuple[List, List]:
    """Compute theoretical Irwin-Hall convolution PDFs for t in [1, 2, 4, 8]."""
    q = pow(p, e) + 1
    m1, m2, m3 = pow(p, r), pow(p, ePrime), q
    n1_th, n2_th, n3_th = m1 // 2, m2 // 2, m3 // 2
    numer1, numer2, numer3 = q - 1, q, 1
    denom = q
    scale1, scale2, scale3 = numer1 / denom, numer2 / denom, numer3 / denom

    t_candidates = [1, 2, 4, 8]
    c0_list = []
    pdf_list = []
    num_list = [80, 40, 20, 10]

    for t in t_candidates:
        idx_t = int(np.log2(t))
        num = num_list[idx_t]
        dx = 1 / (q * num)

        max_val1 = t * scale1 * n1_th
        n_points1 = 2 * t * num * numer1 * n1_th + 1
        c1 = np.linspace(-max_val1, max_val1, n_points1)

        pdf1 = np.array([_irwin_hall_density(x, t, n1_th) for x in c1])
        pdf1 /= np.trapezoid(pdf1, c1)

        s2, p2 = _irwin_hall_discrete_pmf(t, m2)
        s3, p3 = _irwin_hall_discrete_pmf(t, m3)
        g2, pp2 = _expand_pmf_to_modQ_grid(s2, p2, numer2)
        g3, pp3 = _expand_pmf_to_modQ_grid(s3, p3, numer3)

        pmf23 = np.convolve(pp2, pp3)
        center = len(pmf23) // 2
        supp23 = np.arange(-center, len(pmf23) - center) / denom

        n_points23 = 2 * t * (numer2 * n2_th + numer3 * n3_th) * num + 1
        max_val23 = t * (n2_th * scale2 + n3_th * scale3)
        c23 = np.linspace(-max_val23, max_val23, n_points23)

        pdf23 = np.interp(c23, supp23, pmf23, left=0, right=0)
        pdf23 /= np.trapezoid(pdf23, c23)

        pdf0 = fftconvolve(pdf1, pdf23, mode='full') * dx
        c0 = np.linspace(c1[0] + c23[0], c1[-1] + c23[-1], len(pdf0))

        c0_list.append(c0)
        pdf_list.append(pdf0)

    return c0_list, pdf_list


def plot_coeff_distributions(
    ctxt_lists: np.ndarray,
    irwin_hall_dict: List[Dict],
    sk_pos: List[int],
    sk_neg: List[int],
) -> None:
    """
    Reproduce cells 5 and 6 from BGV-BTS-KRD-coeff-plot.ipynb:
    compute and plot the cj distribution for each (sj, ell) combination.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required. Install with: pip install matplotlib"
        ) from exc

    list_fail_slots = ctxt_lists[:, 0].astype(np.int64)
    list_ctxt_a_pwfl = ctxt_lists[:, 2:]

    sk_pos_set = set(sk_pos)
    sk_neg_set = set(sk_neg)

    def get_si(idx):
        if idx in sk_pos_set:
            return 1
        elif idx in sk_neg_set:
            return -1
        else:
            return 0

    # This stage is the slow one: about 17 minutes per parameter set. Report
    # progress as it goes, so that a long wait does not look like a hang.
    import time as _time
    _t0 = _time.time()
    print("Computing theoretical Irwin-Hall PDFs... "
          "(the whole stage takes about 17 minutes)", flush=True)
    c0_list, pdf_list = _compute_theoretical_pdfs()
    print(f"  done, {_time.time() - _t0:.0f}s elapsed", flush=True)

    t_values = [1, 2, 4, 8]
    si_values = [-1, 0, 1]

    # Select indices per (t, si) combination with balanced counts
    selected_indices = {}
    for t in t_values:
        candidates_by_si = {}
        for si in si_values:
            candidates = [
                idx for idx in range(m)
                if list(irwin_hall_dict[idx].keys())[0] == t and get_si(idx) == si
            ]
            candidates_by_si[si] = candidates

        pos_count = len(candidates_by_si[1])
        neg_count = len(candidates_by_si[-1])
        min_count = min(pos_count, neg_count)
        if min_count == 0:
            min_count = 5

        for si in si_values:
            if candidates_by_si[si]:
                selected_indices[(t, si)] = candidates_by_si[si][:min_count]

    print("Selected indices for plotting (Balanced):")
    for (t, si), idx_list in selected_indices.items():
        print(f"  t={t}, si={si:2}: selected {len(idx_list)} indices")

    # Precompute plot data: aggregate samples, compute KDE and KL divergence
    plot_data = {}
    for (t, si), idx_list in selected_indices.items():
        combined_samples_list = []

        for idx in idx_list:
            target_indices = list(irwin_hall_dict[idx].values())[0]
            for j in range(len(list_fail_slots)):
                if list_fail_slots[j] in target_indices:
                    combined_samples_list.append(list_ctxt_a_pwfl[j][idx])

        if not combined_samples_list:
            print(f"No data for combination (t={t}, si={si})")
            continue

        samples = np.array(combined_samples_list)
        mean_col = np.mean(samples)
        var_col = np.var(samples)

        idx_t = int(np.log2(t))
        c0 = c0_list[idx_t]
        pdf0 = pdf_list[idx_t]
        mean0, var0 = _compute_stats(pdf0, c0)

        try:
            kde = gaussian_kde(samples)
            pdf_samples = kde(c0)
            pdf_samples /= np.trapezoid(pdf_samples, c0)
            kl_div = _compute_kl(pdf_samples, pdf0, c0)
        except Exception as exc:
            print(f"KDE calculation failed for (t={t}, si={si}): {exc}")
            continue

        plot_data[(t, si)] = {
            'idx': f"Aggregated({len(idx_list)})",
            'samples': samples,
            'mean_col': mean_col,
            'var_col': var_col,
            'c0': c0,
            'pdf0': pdf0,
            'pdf_samples': pdf_samples,
            'kl_div': kl_div,
            'mean0': mean0,
            'var0': var0,
            't': t,
            'si': si,
        }

        print(f"Group (t={t}, si={si:2}): Total Indices={len(idx_list)}, "
              f"Total Samples={len(samples)}, KL-div={kl_div:.6f} "
              f"[{_time.time() - _t0:.0f}s elapsed]", flush=True)
        print(f"  Irwin-Hall has mean: {mean0:.3f}, var: {var0:.6f}")
        print(f"  Collected  has mean: {mean_col:.3f}, var: {var_col:.6f}")

    # Plot grouped by t
    plots_by_t = defaultdict(list)
    for key, data in plot_data.items():
        plots_by_t[data['t']].append(data)

    for t, datas in plots_by_t.items():
        if t == 1 or t == 4:
            plt.figure(figsize=(5, 4))
        else:
            plt.figure(figsize=(4.6, 4))

        c0 = datas[0]['c0']
        pdf0 = datas[0]['pdf0']
        plt.plot(c0, pdf0, label=f"Irwin-Hall ($\\ell = {int(np.log2(t))}$)", color='black', linestyle='--', linewidth=2)

        for data in sorted(datas, key=lambda d: (d['si'] != 0, abs(d['si']))):
            si = data['si']
            pdf_samples = data['pdf_samples']
            n_samples = len(data['samples'])

            if si == 0:
                plt.plot(c0, pdf_samples, label=f"Collected ($s_j = {si}$)", color='blue', linestyle='-', linewidth=2.0)
            elif si == 1:
                plt.plot(c0, pdf_samples, label=f"Collected ($s_j = {si}$)", color='red', linestyle=':', linewidth=3.0)
            elif si == -1:
                plt.plot(c0, pdf_samples, label=f"Collected ($s_j = {si}$)", color='green', linestyle=':', linewidth=3.0)

            print(f"For $\\ell = {int(np.log2(t))}$, $s_j = {si}$: {n_samples} samples")

        print(f"└-> Distributions of $(c_j$ | fail$)$ w.r.t. $\\ell= {int(np.log2(t))}$")
        plt.xlabel("$c_j$", fontsize=12)
        if t == 1 or t == 4:
            plt.ylabel("PDF of $c_j$", fontsize=12)
        plt.legend(loc="lower right", fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.show()


# ============================================================================
# Main Analysis
# ============================================================================

def main():
    args = parse_args()

    # Selecting the parameter set rebinds the three module-level constants that
    # the rest of the analysis reads.
    global m, phi_m, nslots, slice_points
    if args.slices:
        a, b, c = (int(x) for x in args.slices.split(":"))
        slice_points = list(range(a, b, c))
    m = args.m
    phi_m = PARAM_SETS[m]["phi_m"]
    nslots = PARAM_SETS[m]["nslots"]

    # Find data directory. The samples ship flat under the parameter set's own
    # directory, a fresh sec3-collect leaves them under Collected_data, and a
    # local copy may keep them one level deeper still, so the subdirectories of
    # each candidate are tried as well.
    stem = f"m{m}_p2_e8_ep4_r1"
    candidate_roots = []
    for parent in (Path.cwd(), Path.cwd() / "BGV", Path.cwd() / "Collected_data"):
        root = parent / stem
        candidate_roots.append(root)
        if root.is_dir():
            candidate_roots.extend(sorted(d for d in root.iterdir() if d.is_dir()))
    # Last the working directory itself, for a run made from inside the data.
    # Its subdirectories are deliberately not searched: they hold the other
    # parameter set, and picking that up silently would be worse than failing.
    candidate_roots.append(Path.cwd())

    # The samples ship as .npy, a quarter the size and read in seconds rather
    # than minutes; a fresh collection leaves the .txt that HElib writes, under
    # a name carrying the time it started, so the newest is the one meant.
    def coeff_in(path):
        for pattern in ("coeff*.npy", "coeff*.txt"):
            found = sorted(path.glob(pattern))
            if found:
                return found[-1]
        return None

    data_root = next((path for path in candidate_roots if coeff_in(path)), None)

    if data_root is None:
        raise FileNotFoundError("Could not find coeff.npy or coeff.txt")

    coeff_path = coeff_in(data_root)

    key_files = sorted(data_root.glob("secret_keys*.txt"))
    if not key_files:
        raise FileNotFoundError(f"Could not find secret_keys*.txt in {data_root}")
    secret_key_file = key_files[-1]

    # The Irwin-Hall matrix is named for m, so asking for it by name is also
    # what keeps a directory holding the other parameter set from being used.
    irwin_candidates = [path for path in (data_root / f"irwin_hall_{m}.npy",
                                          data_root / f"irwin_hall_{m}.txt")
                        if path.exists()]
    if not irwin_candidates:
        raise FileNotFoundError(
            f"Found {coeff_path.name} in {data_root}, but no irwin_hall_{m}: "
            f"the samples there are not those of m={m}")

    irwin_path = irwin_candidates[0]

    print(f"Data root: {data_root.resolve()}")
    print(f"Using coeff: {coeff_path.name}")
    print(f"Using secret keys: {secret_key_file.name}")
    print(f"Using irwin_hall: {irwin_path.name}")
    print()

    # Load reference secret
    sk_pos, sk_neg = parse_secret_keys_file(secret_key_file)
    reference_available = len(sk_pos) > 0 or len(sk_neg) > 0

    if reference_available:
        print(f"Reference secret loaded: hwt={len(sk_pos) + len(sk_neg)}")
        print(f"  sk +1: {sorted(sk_pos)}")
        print(f"  sk -1: {sorted(sk_neg)}")
    else:
        print("No reference secret available")

    print()

    # Load irwin_hall matrix
    print("Loading irwin_hall matrix...")
    matrix_irwin_hall = (np.load(irwin_path) if irwin_path.suffix == ".npy"
                         else np.loadtxt(irwin_path, dtype=np.int16))
    irwin_hall_dict = build_irwin_hall_dict(matrix_irwin_hall)
    print(f"Irwin-Hall matrix shape: {matrix_irwin_hall.shape}")
    print()

    # Load the coefficients. The .npy is memory-mapped: nothing below writes to
    # it, and the pages a slice touches are the only ones that come in.
    print(f"Loading {coeff_path.name}...")
    ctxt_lists_full = (np.load(coeff_path, mmap_mode="r")
                       if coeff_path.suffix == ".npy"
                       else np.loadtxt(coeff_path, dtype=np.float32))
    print(f"Total rows in {coeff_path.name}: {len(ctxt_lists_full)}")
    print()

    # --- Plotting Logic ---
    if args.plot_vars:
        print("Computing full-data variance groups for plotting...", flush=True)
        v1, v2, v4, v8 = compute_variances_optimized(
            ctxt_lists_full, irwin_hall_dict, m
        )
        actual_sk_indices = sorted(set(sk_pos) | set(sk_neg))
        print(f"  Highlighting {len(actual_sk_indices)} actual secret indices in the plot")
        plot_variance_distributions(v1, v2, v4, v8, actual_sk_indices)
        return  # Exit after plotting if requested

    if args.plot_coeff_dist:
        plot_coeff_distributions(ctxt_lists_full, irwin_hall_dict, sk_pos, sk_neg)
        return  # Exit after plotting if requested

    print("=" * 80)
    print(f"Analysis with configurations: n1={n1}, n2={n2}, n4={n4}, n8={n8}")
    print("=" * 80)
    print()

    # Extract slots once for all slices
    list_fail_slots_full = ctxt_lists_full[:, 0].astype(np.int32)

    # Analyze each slice
    results_by_slice = {}

    # Calculate total ciphertexts correctly
    diffs = np.diff(list_fail_slots_full)
    total_ctxts = len(np.where(diffs < 0)[0]) + 1
    total_rows = count_ciphertexts_until(list_fail_slots_full, total_ctxts)
    print(f"\nTotal ciphertexts: {total_ctxts}\n")

    terminate = False
    for num_ctxts in slice_points:
        if terminate:
            print(f"  Exit: All {total_ctxts} ciphertexts used...")
            break

        print(f"Processing first {num_ctxts} ciphertexts...")

        # Find row index for this number of ciphertexts
        row_idx = count_ciphertexts_until(list_fail_slots_full, num_ctxts)

        print(f"Using {row_idx} among {total_rows} failed slots" + (f" (reaching max, only {total_ctxts} ctxts used)" if row_idx == total_rows else ""))
        if row_idx == 0:
            print(f"  Warning: Could not find {num_ctxts} ciphertexts, skipping...")
            continue

        if row_idx == total_rows:
            terminate = True

        # Slice data
        ctxt_lists = ctxt_lists_full[:row_idx]

        # Compute variances
        var_sum1, var_sum2, var_sum4, var_sum8 = compute_variances_optimized(
            ctxt_lists, irwin_hall_dict, m
        )

        print(f"  Used {row_idx} rows ({len(ctxt_lists)} data points)")

        # Select candidates
        estimated_nonzero = select_candidates(
            var_sum1, var_sum2, var_sum4, var_sum8,
            n1, n2, n4, n8
        )

        # Evaluate
        eval_result = evaluate(estimated_nonzero, sk_pos, sk_neg)
        results_by_slice[num_ctxts] = {
            'estimated': estimated_nonzero,
            'eval': eval_result,
            'rows_used': row_idx
        }

        # Print results
        print(f"  Estimated candidates: {len(estimated_nonzero)}")
        if eval_result['reference_available']:
            print(f"    Correct: {eval_result['correct']}/{eval_result['total_true']}")
            print(f"    Recall: {eval_result['recall']:.2%}")
            print(f"    Precision: {eval_result['precision']:.2%}")
        print()

    # Summary table
    print("=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print(f"{'Ctxts':>8} {'Rows':>8} {'Candidates':>12} {'Correct':>10} {'Recall':>10} {'Precision':>10}")
    print("-" * 70)

    for num_ctxts in slice_points:
        if num_ctxts not in results_by_slice:
            continue
        result = results_by_slice[num_ctxts]['eval']
        rows_used = results_by_slice[num_ctxts]['rows_used']
        correct_str = f"{result['correct']}/{result['total_true']}" if result['reference_available'] else "N/A"
        recall_str = f"{result['recall']:.2%}" if result['reference_available'] else "N/A"
        precision_str = f"{result['precision']:.2%}" if result['reference_available'] else "N/A"

        print(f"{num_ctxts:>8} {rows_used:>8} {result['total_candidates']:>12} {correct_str:>10} {recall_str:>10} {precision_str:>10}")


if __name__ == "__main__":
    main()
