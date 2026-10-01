"""Two adjustments that let the analysis scripts run unmodified in a container.

The scripts end in plt.show(), which is right at a desk and useless in a
container: with no display the figures are drawn and then discarded. A .pth file
alongside this module imports it at interpreter startup, so setting FIGDIR
redirects show() to savefig() without touching the scripts themselves.

It cannot be called sitecustomize: Debian ships one of its own in
/usr/lib/python3.12, which comes first on sys.path and would shadow this.

The second is np.loadtxt, which takes the better part of an hour on the 6.5 GB
coeff.txt. If to_npy.py has been run, the .npy sits beside the .txt and is used
instead; it is memory-mapped, so the analysis needs far less RAM as well.
"""

import os

_figdir = os.environ.get("FIGDIR")

if _figdir:
    import matplotlib

    matplotlib.use("Agg")
    import sys
    import matplotlib.pyplot as _plt

    # FIGNAME names the figure, and the Makefile sets it to whatever the paper
    # calls the thing being drawn. Without it, fall back to the script name plus
    # its arguments, which at least keeps two invocations from overwriting each
    # other's output.
    _stem = os.environ.get("FIGNAME")
    if not _stem:
        _stem = os.path.splitext(os.path.basename(sys.argv[0] or "figure"))[0]
        _stem += "".join(
            "-" + "".join(c if c.isalnum() else "" for c in a.lstrip("-"))
            for a in sys.argv[1:]
        )
    _count = [0]

    def _show(*args, **kwargs):
        os.makedirs(_figdir, exist_ok=True)
        for num in _plt.get_fignums():
            _count[0] += 1
            path = os.path.join(_figdir, "%s-%02d.pdf" % (_stem, _count[0]))
            _plt.figure(num).savefig(path, bbox_inches="tight")
            print("[saved] %s" % path, flush=True)
        _plt.close("all")

    _plt.show = _show


def _prefer_npy():
    """Read the .npy written by to_npy.py in place of the .txt, when present."""
    import os.path
    import numpy as _np

    _loadtxt = _np.loadtxt

    def loadtxt(fname, *args, **kwargs):
        try:
            cache = os.path.splitext(os.fspath(fname))[0] + ".npy"
        except TypeError:      # a file object or an iterable of lines
            return _loadtxt(fname, *args, **kwargs)
        if os.path.exists(cache):
            print("[npy] %s" % cache, flush=True)
            return _np.load(cache, mmap_mode="r")
        return _loadtxt(fname, *args, **kwargs)

    _np.loadtxt = loadtxt


try:
    _prefer_npy()
except ImportError:            # numpy absent; nothing to do
    pass
