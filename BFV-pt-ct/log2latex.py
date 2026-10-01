import re
import sys
import math
import os
import argparse

LOG2_10 = math.log2(10.0)


def log2_abs_from_num_str(num_str):
    s = num_str.strip()
    if not s:
        raise ValueError("Empty numeric string")

    # Parse scientific notation directly to avoid float overflow for huge exponents.
    if "e" in s.lower():
        coeff_str, exp_str = re.split(r"[eE]", s, maxsplit=1)
        coeff = float(coeff_str)
        exp10 = int(exp_str)
        if coeff == 0.0:
            return float("-inf")
        return math.log2(abs(coeff)) + exp10 * LOG2_10

    val = float(s)
    if val == 0.0:
        return float("-inf")
    return math.log2(abs(val))


def log2_hypot_from_components(a_str, b_str):
    la = log2_abs_from_num_str(a_str)
    lb = log2_abs_from_num_str(b_str)

    if la == float("-inf") and lb == float("-inf"):
        return float("-inf")

    m = max(la, lb)
    term_a = 0.0 if la == float("-inf") else 2.0 ** (2.0 * (la - m))
    term_b = 0.0 if lb == float("-inf") else 2.0 ** (2.0 * (lb - m))
    return m + 0.5 * math.log2(term_a + term_b)

def format_log2_value(v):
    if abs(v - round(v)) < 1e-9:
        return f"{int(round(v))}"
    return f"{v:.3f}"


def parse_log_and_generate_latex(file_path, new_mode=False, small_mode=False,
                                 render_pdf=False):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            log_text = f.read().replace('\xa0', ' ')
    except FileNotFoundError:
        return f"Error: File not found: {file_path}"

    def safe_extract(pattern, text, name, flags=re.IGNORECASE | re.MULTILINE):
        match = re.search(pattern, text, flags)
        if not match:
            raise ValueError(f"'{name}' No info.")
        return match

    try:
        # 1. Parse CONFIGURATION SUMMARY
        n_val = safe_extract(r'Ring\s+Dimension\s+\(N\)\s*:\s*(\d+)', log_text, "N").group(1)
        n_int = int(n_val)
        log2q = float(safe_extract(r'Log2\s+Q\s+\(Ciphertext\)\s*:\s*([\d.]+)', log_text, "Log2 Q").group(1))
        p_val = int(safe_extract(r'Plaintext\s+Modulus\s+\(p\)\s*:\s*(\d+)', log_text, "p").group(1))
        log2p = math.log2(p_val)
        
        vnorm = float(safe_extract(r'Vnorm\s*=\s*([\d.]+)', log_text, "Vnorm").group(1))
        factor_match = re.search(r'(?:Estimated\s+)?Factor\s*=\s*([\d.]+)', log_text, re.IGNORECASE | re.MULTILINE)
        factor = float(factor_match.group(1)) if factor_match else None
        
        idx_match = re.search(r'idx\s*=\s*(\d+)', log_text)
        idx_val = idx_match.group(1) if idx_match else "7"

        # 2. Parse PT info (legacy/new)
        if not new_mode:
            pt_inf_log2 = safe_extract(r'\|\|pt\|\|_inf\s*:\s*2\^([\d.]+)', log_text, "PT Inf").group(1)
            pt_w_log2 = safe_extract(r'\|pt\(w\)\|\s*:\s*2\^([\d.]+)', log_text, "PT w").group(1)
        else:
            constants_match = re.search(
                r'Constants\s*\(c1,\s*c2\)\s*:\s*\((\d+)\s*,\s*(\d+)\)\s*,\s*depthChange\s*:\s*(\d+)',
                log_text,
                re.IGNORECASE | re.MULTILINE
            )
            if constants_match:
                c1_val = int(constants_match.group(1))
                c2_val = int(constants_match.group(2))
                depth_change = int(constants_match.group(3))
            else:
                pt_pair_match = safe_extract(
                    r'\|\|pt1\|\|_inf\s*:\s*2\^(\d+(?:\.\d+)?)\s*\|\s*\|\|pt2\|\|_inf\s*:\s*2\^(\d+(?:\.\d+)?)',
                    log_text,
                    "PT1/PT2 inf"
                )
                depth_change_match = safe_extract(
                    r'Switch\s+depth\s*\(.*?\)\s*:\s*(\d+)',
                    log_text,
                    "Switch depth"
                )
                log2_c1 = float(pt_pair_match.group(1))
                log2_c2 = float(pt_pair_match.group(2))
                depth_change = int(depth_change_match.group(1))
                c1_val = None
                c2_val = None

            if c1_val is not None and c2_val is not None:
                log2_c1 = math.log2(c1_val)
                log2_c2 = math.log2(c2_val)

        # Determine file type (general vs custom) - needed early for factor calculation
        is_general_file = os.path.basename(file_path).lower().startswith("general")

        # If Factor is missing in logs, derive it from current settings.
        # custom: log2(||pt||_inf) + 1 + 0.5*log2(N)
        # general: log2(ptm/2) + 1 + 0.5*log2(N)
        if factor is None:
            if is_general_file:
                factor = math.log2(p_val / 2.0) + 1.0 + 0.5 * math.log2(n_int)
            else:
                if not new_mode:
                    factor = float(pt_inf_log2) + 1.0 + 0.5 * math.log2(n_int)
                else:
                    factor = float(log2_c1) + 1.0 + 0.5 * math.log2(n_int)

        # 3. Extract noise data
        d0_linf = safe_extract(r'Initial ct0 info.*?real e_i \(L-inf\)\s*:\s*([\d.e+-]+)', log_text, "Initial L-inf", re.DOTALL).group(1)
        d0_ew_match = safe_extract(r'Initial ct0 info.*?real e_i\(w\)\s*:\s*\(([\d.e+-]+)\s*,\s*([\d.e+-]+)\)', log_text, "Initial e_w", re.DOTALL)
        d0_ew_log2 = log2_hypot_from_components(d0_ew_match.group(1), d0_ew_match.group(2))

        linf_matches = re.findall(r'EvalNoiseBFVComplex running.*?real e_i \(L-inf\)\s*:\s*([\d.e+-]+)', log_text, re.DOTALL)
        ew_matches = re.findall(r'EvalNoiseBFVComplex running.*?real e_i\(w\)\s*:\s*\(([\d.e+-]+)\s*,\s*([\d.e+-]+)\)', log_text, re.DOTALL)
        
        final_linf = [float(d0_linf)] + [float(m) for m in linf_matches]
        final_ew_log2 = [d0_ew_log2] + [log2_hypot_from_components(m[0], m[1]) for m in ew_matches]

        # 4. Prepare threshold and coordinate calculations
        log2_q_over_2 = log2q - 1
        log2_q_over_2t = log2q - 1 - log2p
        ymax_log2 = int(log2q + 30) if small_mode else int(log2q + 60)
        ticks = [30 * i for i in range(1, (ymax_log2 // 30) + 1)]
        ytick_str = ", ".join([f"2^{t}" for t in ticks])
        yticklabel_str = ", ".join([str(t) for t in ticks])
        axis_width = "0.70\\textwidth" if small_mode else "0.95\\textwidth"
        axis_height = "220pt" if small_mode else "300pt"
        legend_pos = "south east" if small_mode else "north west"
        legend_prefix = "\\scriptsize " if small_mode else ""
        node_font = "\\scriptsize" if small_mode else "\\footnotesize"
        
        # 5. Auto-generate coordinates and arrows
        coords_bound, coords_linf, coords_ew = "", "", ""
        arrow_code = "% --- Automatically Generated Arrows ---\n"
        
        # Prepare legend entry for orange line (L-inf bound)
        bound_legend = "$L^\\infty$-bound of $e_d$"
        # if is_general_file:
        #     bound_legend = "$L^\\infty$-bound of $e_d$ w.r.t. any $\\pt \\in \\ring_t$"
        # else:
        #     bound_legend = "$L^\\infty$-bound of $e_d$ w.r.t. $\\infnorm{\\pt}$"

        # 5-2. Real Dec. Fail (when olive points lie between q/2t and q/2)
        real_fail_points = []
        for d in range(len(final_linf)):
            val = final_linf[d]
            coords_linf += f"({d}, {{2^{val:.2f}}}) "
            coords_ew += f"({d}, {{2^{final_ew_log2[d]:.2f}}}) "
            
            # Real Dec. Fail condition: log2(q/2t) < measured L-inf < log2(q/2)
            if log2_q_over_2t < val < log2_q_over_2 - 1:
                real_fail_points.append((d, val))

        # Generate arrows (there may be multiple)
        avg_x = sum(p[0] for p in real_fail_points) / len(real_fail_points) if real_fail_points else 0.0
        if (len(real_fail_points) == 1):
            avg_x -= 0.5  # Shift left if only one point to avoid overlap with the point itself
        if real_fail_points:
            for i, (d, val) in enumerate(real_fail_points):
                label = "\\scriptsize Real\\\\\\scriptsize Dec. Fail." if (i == 0 and small_mode) else ("Real\\\\ Dec. Fail." if i == 0 else "")
                arrow_code += f"    \\draw[<-, >=stealth] (axis cs:{d}, {{2^{(val-2.0):.2f}}}) -- (axis cs:{avg_x}, {{2^{log2_q_over_2t-40}}}) node[below, font={node_font}, align=center] {{{label}}};\n"
        
        # 5-1. Detection logic for Estimated Dec. Fail
        est_fail_d = -1
        for d in range(len(final_linf) + 15):
            bound_val = vnorm + (d * factor)
            coords_bound += f"({d}, {{2^{bound_val:.2f}}}) "
            
            # Condition: capture only the first point where bound_val exceeds log2_q_over_2t (red line)
            if est_fail_d == -1 and bound_val > log2_q_over_2t:
                est_fail_d = d
                # Generate arrow pointing exactly to that point (d)
                est_label = "\\scriptsize Estimated\\\\\\scriptsize Dec. Fail." if small_mode else "Estimated\\\\Dec. Fail."
                arrow_code += f"    \\draw[<-, >=stealth] (axis cs:{d}, {{2^{(bound_val-2.0):.2f}}}) -- (axis cs:{avg_x+1}, {{2^{bound_val-100}}}) node[left, font={node_font}, align=center] {{{est_label}}};\n"
            
            # Stop after enough points to avoid making the graph too long
            if bound_val > log2_q_over_2t: 
                break
        
        last_d_bound = d
        dynamic_xmax = max(len(final_linf) - 0.5, last_d_bound + 0.5)

        # 5-3. The same three series and two thresholds, as a PDF a reviewer can
        # open directly (the pgfplots above needs LaTeX to become a picture).
        # plt.show() is redirected to savefig by indcpad_figsave when FIGDIR is set.
        if render_pdf:
            import matplotlib.pyplot as plt
            xs = list(range(len(final_linf)))
            bound = [vnorm + d * factor for d in xs]
            plt.figure(figsize=(6, 4))
            plt.plot(xs, bound, marker='^', ls='--', color='orange',
                     label=r'estimated $L^\infty$-bound of $e_d$')
            plt.plot(xs, final_ew_log2, marker='s', color='blue',
                     label=r'$\log_2 |e_d(\omega)|$')
            plt.plot(xs, final_linf, marker='D', color='olive',
                     label=r'$\log_2 \|e_d\|_\infty$')
            plt.axhline(log2_q_over_2, color='red', ls='--', lw=0.8)
            plt.axhline(log2_q_over_2t, color='purple', ls='--', lw=0.8)
            plt.text(0, log2_q_over_2, r'$\log_2 q/2$', color='red',
                     va='bottom', fontsize=8)
            plt.text(0, log2_q_over_2t, r'$\log_2 q/(2t)$', color='purple',
                     va='bottom', fontsize=8)
            plt.xlabel(r'PT-CT multiplicative depth $d$')
            plt.ylabel(r'noise magnitude ($\log_2$)')
            plt.title('BFV pt-ct noise growth (Figure 3)')
            plt.legend(loc='upper left', fontsize=8)
            plt.grid(True, ls=':', alpha=0.5)
            plt.show()

        # 6. Return final template
        dynamic_mark_size = 1.5 if len(final_linf) > 15 else 1.3

        if new_mode:
            caption = (
                f"Noise Analysis for Algorithm~\\ref{{alg:BFV-simple}}. "
                f"$N={n_val}$, $q \\approx 2^{{{log2q:.0f}}}$, $t \\approx 2^{{{log2p:.0f}}}$, index $i={idx_val}$. "
                f"In the noise estimation (\\textcolor{{orange}}{{$\\blacktriangle$}}), $\\infnorm{{e_{{\\mathsf{{sqr}}, d}}}} \\le (\\delta_N \\cdot t / 2)^d \\cdot V_N$ is used. "
                f"In the attack, $\\lVert \\mathsf{{pt}} \\rVert_\\infty = 2^{{{format_log2_value(log2_c1)}}}$ "
                f"is used for the first ${depth_change - 1}$ depths, then "
                f"$\\lVert \\mathsf{{pt}} \\rVert_\\infty = 2^{{{format_log2_value(log2_c2)}}}$ is used.~\\label{{fig:pt-ct-error-{n_val}-{log2q:.0f}-new}}"
            )
        else:
            caption = (
                f"Noise Analysis for Algorithm~\\ref{{alg:BFV-simple}}. "
                f"$N={n_val}$, $q \\approx 2^{{{log2q:.0f}}}$, $t \\approx 2^{{{log2p:.0f}}}$, index $i={idx_val}$, "
                f"$\\lVert \\mathsf{{pt}} \\rVert_\\infty = 2^{{{pt_inf_log2}}}$, and "
                f"$|\\mathsf{{pt}}(\\omega)| \\approx 2^{{{pt_w_log2}}}$. "
                f"{'The estimated $L^\\infty$-bound is for any $\\pt \\in \\ring_t$.' if is_general_file else f'The estimated $L^\\infty$-bound is specific to $\\infnorm{{\\pt}} \\le 2^{{{pt_inf_log2}}}$.'} "
                f"\\label{{fig:pt-ct-error-{n_val}-{log2q:.0f}-{'general' if is_general_file else 'custom'}}}"
            )


        return f"""
\\begin{{figure}}[ht]
\\centering
\\begin{{tikzpicture}}
    \\begin{{axis}}[
        width={axis_width}, height={axis_height},
        xlabel={{PT-CT Multiplicative Depths ($d$)}},
        ylabel={{Magnitude of Noise ($e_d$, in $\\log_2$)}},
        ymode=log, log basis y={{2}},
        ytick={{{ytick_str}}},
        yticklabels={{{yticklabel_str}}},
        grid=both, grid style={{dashed, gray!30}},
        xmin=-0.5, xmax={dynamic_xmax},
        ymin=1, ymax={{2^{ymax_log2}}}, 
        legend pos={legend_pos},
    ]
    \\addplot[color=orange, mark=triangle*, mark options={{solid, line width={dynamic_mark_size}pt}}, mark size={dynamic_mark_size}pt, dashed] coordinates {{{coords_bound}}};
    \\addlegendentry{{{legend_prefix}{bound_legend}}}
    
    \\addplot[color=blue, mark=square*, mark options={{solid, line width={dynamic_mark_size}pt}}, mark size={dynamic_mark_size}pt] coordinates {{{coords_ew}}};
    \\addlegendentry{{{legend_prefix}Measured $|e_d(\\omega)|$}}
    
    \\addplot[color=olive, mark=diamond*, mark options={{solid, line width={dynamic_mark_size}pt}}, mark size={dynamic_mark_size}pt] coordinates {{{coords_linf}}};
    \\addlegendentry{{{legend_prefix}Measured $\\lVert e_d \\rVert_\\infty$}}
    
    % --- Threshold Lines ---
    \\draw[red, dashed] (axis cs:-0.5, {{2^{log2_q_over_2}}}) -- (axis cs:{dynamic_xmax}, {{2^{log2_q_over_2}}}) 
        node[pos=0.01, anchor=north west] {{{legend_prefix}$\\log_2 q/2 \\approx {log2_q_over_2:.1f}$}};
    \\draw[purple, dashed] (axis cs:-0.5, {{2^{log2_q_over_2t}}}) -- (axis cs:{dynamic_xmax}, {{2^{log2_q_over_2t}}}) 
        node[pos=0.01, anchor=north west] {{{legend_prefix}$\\log_2 q/2t \\approx {log2_q_over_2t:.1f}$}};

    {arrow_code}
    \\end{{axis}}
\\end{{tikzpicture}}
\\caption{{{caption}}}
%% File: {os.path.basename(file_path)}
\\end{{figure}}
"""
    except Exception as e:
        return f"Error: {e}"

def format_time_str(t):
    if t >= 1000:
        return f"{t:.0f}s"
    elif t >= 10:
        return f"{t:.1f}s"
    elif t >= 1:
        return f"{t:.2f}s"
    else:
        return f"{t:.3f}s"


def format_ptm_latex(p):
    pm1 = p - 1
    if pm1 > 0 and (pm1 & (pm1 - 1)) == 0:
        k = int(round(math.log2(pm1)))
        return f"$2^{{{k}}}+1$"
    log2p = math.log2(p)
    if abs(log2p - round(log2p)) < 0.005:
        return f"$\\approx 2^{{{int(round(log2p))}}}$"
    return f"$\\approx 2^{{{log2p:.1f}}}$"


def format_pow2_latex(exp_str):
    v = float(exp_str)
    if abs(v - round(v)) < 1e-9:
        return f"$2^{{{int(round(v))}}}$"
    return f"$2^{{{v:.1f}}}$"


def parse_file_for_table(file_path, is_general):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            text = f.read().replace('\xa0', ' ')
    except FileNotFoundError:
        return None

    def grep(pattern, flags=re.IGNORECASE | re.MULTILINE):
        return re.search(pattern, text, flags)

    try:
        N = int(grep(r'Ring\s+Dimension\s+\(N\)\s*[:\s]+(\d+)').group(1))
        log2q = float(grep(r'Log2\s+Q\s+\(Ciphertext\)\s*:\s*([\d.]+)').group(1))
        p = int(grep(r'Plaintext\s+Modulus\s+\(p\)[:\s]+(\d+)').group(1))
        d_allow = int(grep(r'Theoretical\s+Est\.\s+Depth\s*:\s*(\d+)').group(1))
        pt_inf_exp = grep(r'\|\|pt\|\|_inf\s*:\s*2\^([\d.]+)').group(1)
        pt_w_exp = grep(r'\|pt\(w\)\|\s*:\s*2\^([\d.]+)').group(1)
        t_col = float(grep(r'Collecting\s+Time\s+\(e\(w\)\)\s*:\s*([\d.]+)').group(1))
        t_sk = float(grep(r'Recovery\s+Time\s+\(e\)\s*:\s*([\d.]+)').group(1))

        if is_general:
            d_attack = d_allow
            o_dec_str = "$N$"
        else:
            d_attack = int(grep(r'Theoretical\s+ATK\.\s+Depth\s*:\s*(\d+)').group(1))
            m = grep(r'last_iter\s+range\s+\(main\s+loop\)\s*:\s*\[(\d+),\s*(\d+)\]')
            if m:
                count = int(m.group(2)) - int(m.group(1)) + 1
                o_dec_str = f"${count}N$" if count > 1 else "$N$"
            else:
                o_dec_str = "$N$"

        log2N = int(round(math.log2(N)))
        return {
            'N': N,
            'N_str': f"$2^{{{log2N}}}$",
            'q_str': f"$\\approx 2^{{{log2q:.0f}}}$",
            'ptm_str': format_ptm_latex(p),
            'pt_inf_str': format_pow2_latex(pt_inf_exp),
            'pt_w_str': format_pow2_latex(pt_w_exp),
            'd_allow': d_allow,
            'd_attack': d_attack,
            'o_dec_str': o_dec_str,
            't_col_str': format_time_str(t_col),
            't_sk_str': format_time_str(t_sk),
            'is_general': is_general,
        }
    except Exception as e:
        print(f"Warning: could not parse {file_path}: {e}", file=sys.stderr)
        return None


def generate_table(folder="."):
    import glob
    custom_files = sorted(glob.glob(os.path.join(folder, "custom_*.txt")))
    general_files = sorted(glob.glob(os.path.join(folder, "general_*.txt")))

    rows = []
    for f in custom_files:
        row = parse_file_for_table(f, is_general=False)
        if row:
            rows.append(row)
    for f in general_files:
        row = parse_file_for_table(f, is_general=True)
        if row:
            rows.append(row)

    rows.sort(key=lambda r: (r['N'], r['is_general']))

    lines = [
        r"\begin{tabular}{@{}ccc|ccccc|ccc@{}}",
        r"\toprule",
        r"$N$ & $q$ & $t$ & $\infnorm{\pt}$ & $|\pt(\omega)|$ & $d_{\sf allow}$ & $d_{\sf attack}$ & $\#\calO_{\dec}$ & $T_{\sf col}$ & $T_{\sf sk}$ \\ \midrule",
    ]
    for row in rows:
        d_allow_str = f"$~\\,{row['d_allow']}^*$" if row['is_general'] else f"${row['d_allow']}$"
        line = (
            f"{row['N_str']} & {row['q_str']} & {row['ptm_str']} & "
            f"{row['pt_inf_str']} & {row['pt_w_str']} & "
            f"{d_allow_str} & ${row['d_attack']}$ & "
            f"{row['o_dec_str']} & {row['t_col_str']} & {row['t_sk_str']} \\\\"
        )
        lines.append(line)
    lines += [r"\bottomrule", r"\end{tabular}"]
    print("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse BFV logs and output LaTeX figure code.")
    parser.add_argument("filename", nargs="?", help="Input log filename (omit when using --table)")
    parser.add_argument("--new", action="store_true", dest="new_mode", help="Use new input format with c1/c2 and depthChange")
    parser.add_argument("--small", action="store_true", dest="small_mode", help="Use compact graph size and scriptsize in-graph texts")
    parser.add_argument("--table", action="store_true", help="Generate LaTeX table from all custom_*.txt and general_*.txt in the current folder")
    parser.add_argument("--plot", action="store_true", help="Render the figure as a PDF (via FIGDIR) instead of printing the pgfplots code")
    args = parser.parse_args()

    if args.table:
        folder = os.path.dirname(os.path.abspath(__file__))
        generate_table(folder)
    else:
        if not args.filename:
            parser.error("filename is required unless --table is specified")
        fname = args.filename
        if not os.path.exists(fname) and not fname.endswith(".txt"):
            fname += ".txt"
        if args.plot:
            parse_log_and_generate_latex(fname, new_mode=args.new_mode,
                                         small_mode=args.small_mode, render_pdf=True)
        else:
            print(parse_log_and_generate_latex(fname, new_mode=args.new_mode, small_mode=args.small_mode))