#!/usr/bin/env python3
"""Convert the collected BGV samples from text to .npy, once.

np.loadtxt spends the better part of an hour on the 6.5 GB coeff.txt, which is
no way to hand an artifact to a reviewer. The same array as .npy is about a
third of the size, is memory-mapped rather than read, and loads in seconds. The
analysis script prefers it over the text whenever both are there, so a fresh
a gathering only has to be run through this once.

    python3 to_npy.py /path/to/data-collected-by-authors

Conversion is lossless: the dtypes are the ones the analysis asks for, float32
for the coefficients and int16 for the Irwin-Hall matrix.
"""

import sys
import time
from pathlib import Path

import numpy as np

# (glob, dtype) pairs, matching the two np.loadtxt calls in the analysis script.
TARGETS = [("coeff*.txt", np.float32), ("irwin_hall_*.txt", np.int16)]


def convert(txt: Path, dtype) -> None:
    npy = txt.with_suffix(".npy")
    if npy.exists():
        print("  %s is already there, skipping" % npy.name)
        return
    start = time.time()
    print("  %s -> %s ..." % (txt.name, npy.name), end="", flush=True)
    array = np.loadtxt(txt, dtype=dtype)
    # Write beside the source under a temporary name, so an interrupted run
    # cannot leave a half-written cache that later looks valid.
    tmp = npy.parent / (npy.name + ".partial")
    with open(tmp, "wb") as handle:      # a file object, so np.save does not
        np.save(handle, array)           # append a second .npy to the name
    tmp.rename(npy)
    print(" %s, %.0f MB, %.0f s" % (array.shape, npy.stat().st_size / 2**20,
                                    time.time() - start))


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "data-collected-by-authors")
    if not root.is_dir():
        print("no such directory: %s" % root, file=sys.stderr)
        return 1
    found = False
    for pattern, dtype in TARGETS:
        for txt in sorted(root.rglob(pattern)):
            found = True
            convert(txt, dtype)
    if not found:
        print("found no coeff*.txt or irwin_hall_*.txt under %s" % root,
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
