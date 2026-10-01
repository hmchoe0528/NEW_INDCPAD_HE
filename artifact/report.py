#!/usr/bin/env python3
"""Summarise what a run produced, next to what the paper reports.

The evaluator should not have to read raw logs to find out whether a claim came
out. This walks the results directory, pulls the few numbers each attack is
judged on, and writes results/REPORT.md with the paper's value beside it.

    python3 report.py [results-dir]
"""

import os
import re
import sys

# What the paper reports, keyed by the file a run writes it to.
PAPER = {
    "C1": {                         # Table 2, at the ciphertext counts published
        21845: {50: ("2,225", "7/32"), 150: ("6,948", "25/32"),
                250: ("11,718", "32/32")},
        28679: {50: ("2,441", "8/32"), 200: ("10,265", "27/32"),
                350: ("17,466", "32/32")},
    },
}


# Every experimental figure and table of the paper, under the claim it belongs
# to, with the target that produces it and the file that shows it was produced.
# Anything analytic or purely descriptive is left out, since no run covers it.
COVERAGE = [
    ("C1", "Table 2",      "c1", "c1-keyrecovery-m21845.txt"),
    ("C1", "Figure 2",     "c1", "figures/fig2-01.pdf"),
    ("C2", "Table 3",      "c2", "table3.tex"),
    ("C2", "Figure 3",     "c2", "figures/fig3-01.pdf"),
    ("C3", "Figure 4",     "c3", "c3-ctct.txt"),
    ("C3", "Figure 5",     "c3", "figures/fig5-01.pdf"),
    ("full", "Table 4",    "c3",             "c3-ctct.txt"),
    ("full", "Table 5",    "appendix-ctct",  "appendix-multikey-4nt.txt"),
    ("full", "Table 7",    "appendix-ctct",  "appendix-multikey-8nt.txt"),
    ("full", "Figure 6",   "appendix-bgv",   "figures/fig6-01.pdf"),
    ("full", "Figures 7, 8", "appendix-bgv", "figures/fig7-fig8-m21845-01.pdf"),
    ("full", "Figure 9",   "appendix-figs",  "fig9-logN15.tex"),
    ("full", "Figure 10",  "c3",             "c3-ctct.txt"),
    ("full", "Figure 11",  "appendix-ctct",  "figures/fig11-01.pdf"),
]

# What each claim says, for the heading of its section.
CLAIMS = {
    "C1": ("BGV bootstrapping", "Table 2 and Figure 2",
           "The secret key is recovered from the intermediate values of "
           "failing bootstrappings."),
    "C2": ("BFV pt-ct multiplication", "Table 3 and Figure 3",
           "The whole secret key is recovered in N decryption queries, "
           "within a depth the noise bound certifies."),
    "C3": ("BFV ct-ct multiplication", "Figures 4 and 5",
           "Repeated squarings carry the noise past the decryption "
           "threshold while the certified bound stays inside, which breaks "
           "indistinguishability."),
}


def sources(results):
    """{target: 'run' | 'authors'}, from the marks the Makefile leaves."""
    out = {}
    try:
        with open(os.path.join(results, "sources.tsv"), encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2:
                    out[parts[0]] = parts[1]
    except OSError:
        pass
    return out


def coverage(results, complete):
    """One row per figure and table, under the claim it belongs to."""
    out = ["", "## What this run covers", "",
           "Table 1 is a parameter list, Figure 1 describes oracles, and Table 6",
           "is analytic, so no run can cover them.", "",
           "The rows under *appendix* belong to the full version of the paper",
           "rather than the camera-ready, and are not part of the evaluation.",
           "", "| claim | item | produced by | this run |",
           "| --- | --- | --- | --- |"]
    last = None
    for claim, item, target, evidence in COVERAGE:
        path = os.path.join(results, evidence)
        if evidence in complete:
            state = complete[evidence]
        elif os.path.exists(path):
            state = "yes"
        else:
            state = "\u2014"
        # Point at the file, so that "yes" is something the reader can open.
        if state != "\u2014":
            state = "[%s](%s)" % (state, evidence)
        shown = "*appendix*" if claim == "full" else "**%s**" % claim
        out.append("| %s | %s | `%s` | %s |"
                   % ("" if claim == last else shown, item, target, state))
        last = claim
    return out + [""]


def timings(results):
    """(target, seconds) for each target that ran, most recent run of each."""
    path = os.path.join(results, "timings.tsv")
    starts, out = {}, {}
    try:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 3:
                    continue
                target, stamp, which = parts
                try:
                    stamp = int(stamp)
                except ValueError:
                    continue
                if which == "start":
                    starts[target] = stamp
                elif which == "end" and target in starts:
                    out[target] = stamp - starts.pop(target)
    except OSError:
        return []
    return sorted(out.items(), key=lambda kv: -kv[1])


def human(seconds):
    if seconds < 90:
        return "%d s" % seconds
    if seconds < 5400:
        return "%d min %02d s" % (seconds // 60, seconds % 60)
    return "%d h %02d min" % (seconds // 3600, (seconds % 3600) // 60)


def read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            return handle.read()
    except OSError:
        return None


def recall_rows(text):
    """(ciphertexts, rows, correct) from the summary table."""
    out = {}
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+(\d+)\s+\d+\s+(\d+/\d+)\s", line)
        if m:
            out[int(m.group(1))] = (m.group(2), m.group(3))
    return out


def main():
    results = sys.argv[1] if len(sys.argv) > 1 else "results"
    src = sources(results)
    lines = ["# Results", ""]
    lines += ["One section per claim of the camera-ready. Figures are written to",
              "`figures/`. Each section says whether what it shows was computed by",
              "this run or read from the logs and samples the authors collected;",
              "only the first is evidence that the attack works on your machine,",
              "and only the second meets the paper's own numbers to the digit.",
              ""]

    # -- C1 ------------------------------------------------------------------
    got1 = {}
    for m in (21845, 28679):
        text = read(os.path.join(results, "c1-keyrecovery-m%d.txt" % m))
        if text is not None:
            got1[m] = recall_rows(text)
    if got1:
        name, items, what = CLAIMS["C1"]
        authors = src.get("c1", src.get("quick")) == "authors"
        lines += ["## C1 \u2014 %s (%s)" % (name, items), "", what, "",
                  "Results at `c1-keyrecovery-m21845.txt` and",
                  "`c1-keyrecovery-m28679.txt`.", ""]
        if authors:
            lines += ["Recovered here, from the samples the authors collected, so",
                      "the counts are the ones Table 2 publishes.", "",
                      "| set | ciphertexts | slots | correct | paper |",
                      "| --- | --- | --- | --- | --- |"]
        else:
            lines += ["Gathered and recovered here, so these are this run's own",
                      "failures rather than the authors'. What reproduces is the",
                      "recovery, not Table 2's counts, which belong to their",
                      "sample; `make c1 AUTHORDATA=1` reproduces those.", "",
                      "| set | ciphertexts | slots | correct |",
                      "| --- | --- | --- | --- |"]
        for m in (21845, 28679):
            if m not in got1:
                continue
            label = "A*" if m == 21845 else "B*"
            for n in sorted(PAPER["C1"][m]):
                mine = got1[m].get(n)
                cells = [label, str(n), mine[0] if mine else "-",
                         mine[1] if mine else "not reached"]
                if authors:
                    slots, correct = PAPER["C1"][m][n]
                    cells.append("%s, %s" % (slots, correct))
                row = "| " + " | ".join(cells) + " |"
                lines.append(row)
        lines.append("")
        if os.path.exists(os.path.join(results, "figures", "fig2-01.pdf")):
            lines += ["Figure 2 is the four panels `figures/fig2-0*.pdf`: the",
                      "collected coefficients against the Irwin-Hall density, one",
                      "per level.", ""]

    # -- C2 ------------------------------------------------------------------
    got2, cut2 = [], []
    for name in sorted(os.listdir(results)) if os.path.isdir(results) else []:
        if not name.startswith("general_logN"):
            continue
        text = read(os.path.join(results, name))
        if text is None:
            continue
        ratio = re.search(r"SK Coeff Match Ratio\s*:\s*(\S+)", text)
        err = re.search(r"Final L2 Error \(e\)\s*:\s*(\S+)", text)
        calls = re.search(r"O_Dec calls\s*:\s*(\S+)", text)
        if ratio and err and calls:
            got2.append((name, ratio.group(1), err.group(1), calls.group(1)))
        else:
            cut2.append(name)
    if got2:
        name, items, what = CLAIMS["C2"]
        authors = src.get("c2") == "authors"
        if "c2" in src:
            lines += ["## C2 \u2014 %s (%s)" % (name, items), "", what, "",
                      "Results at `general_logN*.txt`, tabulated in `table3.tex`",
                      "and drawn as `figures/fig3-01.pdf`.", ""]
            lines += ["Read from the authors' own logs." if authors else
                      "Computed by this run.", ""]
            draws = read(os.path.join(results, "c2-draws.txt")) or ""
            extra = [l.split("\t") for l in draws.splitlines() if "\t" in l]
            extra = [(a, b) for a, b in extra if b.strip() not in ("", "1")]
            if extra:
                lines += ["Some keys amplify their noise past the window the",
                          "decrypted result has to land in, which is `log2 t` bits",
                          "wide, so the decryption wraps and the recovery cannot",
                          "start. It shows at the first index within seconds, and",
                          "the target simply draws again:", ""]
                lines += ["- %s took %s draws" % (a, b.strip()) for a, b in extra]
                lines.append("")
        else:
            lines += ["### The smoke test's pt-ct run", "",
                      "At the smallest ring dimension only; Table 3 needs all",
                      "five, which is `make c2`.", ""]
        lines += ["The claim is a full recovery: the ratio is 1, the error is 0.", "",
                  "| run | SK match ratio | final L2 error | oracle calls |",
                  "| --- | --- | --- | --- |"]
        for name_, ratio, err, calls in got2:
            lines.append("| `%s` | %s | %s | %s |" % (name_, ratio, err, calls))
        lines.append("")
    if cut2:
        lines += ["The following logs stop before the attack reports, so those",
                  "runs say nothing either way. Either they were interrupted, or",
                  "the key amplified its noise past the window the decrypted",
                  "result has to land in, which is a property of the key rather",
                  "than a fault: run `make c2` again to draw another.", ""]
        lines += ["- `%s`" % n for n in cut2]
        lines.append("")

    # -- C3 ------------------------------------------------------------------
    text = read(os.path.join(results, "c3-ctct.txt"))
    if text is not None:
        name, items, what = CLAIMS["C3"]
        authors = src.get("c3") == "authors"
        lines += ["## C3 \u2014 %s (%s)" % (name, items), "", what, "",
                  "Result at `c3-ctct.txt`, whose per-squaring table is what",
                  "Figure 4 plots, and `figures/fig5-01.pdf`, drawn from it.", ""]
        lines += ["Read from the authors' own log." if authors else
                  "Computed by this run.", ""]
        lines += ["The two distributions are centred at zero and told apart by",
                  "their spread. The gain is fixed in advance, as an adversary",
                  "must fix it, while the key is drawn afresh: about 91% of keys",
                  "fall in the band that separates, so roughly one run in eleven",
                  "comes out flat and has to be repeated.", ""]

    # -- the smoke test's own §5 line ----------------------------------------
    text = read(os.path.join(results, "quick-ctct.txt"))
    if text is not None:
        verdict = [l.strip() for l in text.splitlines() if "key pairs separate" in l]
        if verdict:
            lines += ["### The smoke test's ct-ct run", "",
                      "At a reduced ring dimension, five key pairs, in",
                      "`quick-ctct.txt`:", "", "> " + verdict[-1], "",
                      "One gain is fixed against keys drawn at random, so a pair",
                      "that does not separate is the attack rather than a fault.",
                      "`make c3` shows the same at the paper's parameters.", ""]

    # -- the appendices ------------------------------------------------------
    apx = []
    for target, fname, what in (
            ("appendix-ctct", "appendix-multikey-4nt.txt", "Table 5"),
            ("appendix-ctct", "appendix-multikey-8nt.txt", "Table 7")):
        text = read(os.path.join(results, fname))
        if text is None:
            continue
        verdict = [l.strip() for l in text.splitlines() if "key pairs separate" in l]
        if verdict:
            apx += ["**%s**, in `%s`:" % (what, fname), "", "> " + verdict[-1], ""]
    if apx:
        lines += ["## The appendices", "",
                  "These belong to the full version rather than the camera-ready,",
                  "and we do not ask for them to be evaluated.", ""] + apx

    # -- what did not run ----------------------------------------------------
    missing = []
    for claim, item, target, evidence in COVERAGE:
        if not os.path.exists(os.path.join(results, evidence)):
            missing.append((target, item))
    if missing:
        lines += ["## Not run", "",
                  "Nothing here covers the following; run the target to add it.", ""]
        seen = {}
        for target, item in missing:
            seen.setdefault(target, []).append(item)
        for target, items in seen.items():
            lines.append("- `make %s` \u2014 %s" % (target, ", ".join(items)))
        lines.append("")

    ran = timings(results)
    if ran:
        total = sum(sec for _, sec in ran)
        lines += ["## How long it took", "",
                  "Measured on this machine, by the targets themselves.", "",
                  "| target | wall time |",
                  "| --- | --- |"]
        lines += ["| `%s` | %s |" % (t, human(sec)) for t, sec in ran]
        lines += ["| **total** | **%s** |" % human(total), ""]

    # The coverage map goes first, before the detail.
    complete = {}
    for _, _, _, ev in COVERAGE:
        path = os.path.join(results, ev)
        if not os.path.exists(path):
            continue
        if ev.endswith(".txt"):
            text = read(path) or ""
            # Each program ends differently: the pt-ct attack and ctct_multikey
            # print a closing line, the analysis a summary table, and
            # indcpad_ctct simply stops after its last trial, so for that one
            # count the trials against the number it announced.
            done = ("SK Coeff Match Ratio" in text or "key pairs separate" in text
                    or "SUMMARY TABLE" in text or "Ctxts" in text)
            if not done:
                want = re.search(r"\((\d+) trials,", text)
                seen_n = len(re.findall(r"^\s*\d+ \|", text, re.M))
                done = bool(want) and seen_n >= int(want.group(1))
            complete[ev] = "yes" if done else "interrupted"
    lines[7:7] = coverage(results, complete)

    out = os.path.join(results, "REPORT.md")
    with open(out, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
