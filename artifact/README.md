# Artifact: *New IND-CPA-D Attacks on BFV and BGV*

This artifact backs the three attacks of the paper,

```
Hyeongmin Choe, Minsik Kang, and Damien Stehlé, "New IND-CPA-D Attacks on BFV and BGV," ACM CCS 2026.
```

The artifact includes reproducible (at least the trend of the attack) programs for each section of the paper. 

| Sec. | Attack | Library | Figures & Tables |
| --- | --- | --- | --- |
| 3 | BGV, bootstrapping | HElib @ `3e337a6` | Fig. 2, Tab. 2 |
| 4 | BFV, plaintext–ciphertext (PT-CT) | OpenFHE @ `df495ba` (v1.5.0) | Fig. 3, Tab. 3 |
| 5 | BFV, ciphertext–ciphertext (CT-CT) | OpenFHE @ `2bf4524` (v1.3) | Fig. 4 and 5 |
||

### Additional (full version of the paper)

In the full version of the paper, further experimental results are provided in Appendix C, D and E. 
Note that these may not be found in the proceeding version. 

| App. | Title | Figures & Tables | produced by |
| --- | --- | --- | --- |
| C | More Bootstrapping Experiments | Figs. 6-8 for BGV bootstrapping attack | `appendix-bgv` |
| D | More Error Growth Experiments | Fig. 9, pt-ct noise growth | `appendix-figs` |
| D | | Fig. 10, ct-ct noise growth, unzoomed from Fig. 4 | `c3` log |
| E.1 | More CT-CT Attack Results (Noise Behavior) | Tab. 4 | `c3` log |
| E.2 | More CT-CT Attack Results (Independent Key Pairs) | Tab. 5 | `appendix-ctct` |
| E.3 | More CT-CT Attack Results (Cost in Revised OpenFHE) | Tab. 6-7, Fig. 11 | `appendix-ctct` |
||

## Getting Started

Docker is the only requirement; the image builds HElib and both OpenFHE versions
from source, which takes about five minutes from cold. Build from the top of the
repository, since the Dockerfile needs the whole tree as its context.

```bash
docker build -f artifact/Dockerfile -t indcpad-artifact .
docker run --rm indcpad-artifact make help
```

Then run the attacks, mounting a directory for the results. Sections 4 and 5
need nothing else:

```bash
mkdir -p results
docker run --rm -v "$PWD/results:/artifact/results" \
  indcpad-artifact make c2 c3
```

Section 3 gathers its own failing ciphertexts (about two hours at `THREADS=8`)
and then ranks the candidate positions of the key and draws Fig. 2 (about 20
minutes):

```bash
docker run --rm -v "$PWD/results:/artifact/results" \
  indcpad-artifact make c1 THREADS=8
```

- Read `results/REPORT.md`, which puts each result beside the value the paper
  gives, one section per claim.
- **`make quick` is a smoke test** of all three attacks, a few minutes long, at reduced
  parameters. Its Section 3 part reads samples, so run it only after a `make c1`
  has gathered them under `data-collected-by-authors`; mount that directory and
  pass `AUTHORDATA=1`.
- A plain `c1` gathers its own and needs only the results mount, as `c2` and
  `c3` do. `c1 AUTHORDATA=1` and `appendix-bgv` read a gathered set instead.
- `--build-arg JOBS=n` sets the compile parallelism; it defaults to 4 because
  OpenFHE wants a couple of GB per job.
- Everything is compiled inside the image, so x86-64 and arm64 both work.
- No GPU. What the targets need, measured:

  | | RAM | disk |
  | --- | --- | --- |
  | `quick`, `c2`, `c3` | 4 GB | image only |
  | `appendix` | 8 GB | 4 GB of samples |
  | `c1 AUTHORDATA=1`, `appendix-bgv` | 8 GB | 4 GB of samples |
  | `c1` | 10 GB | 12 GB written, 6 per parameter set |

  The samples are memory-mapped rather than read, so the Sec. 3 analysis itself
  stays near 3 GB whatever their size. Two things want room: the density
  estimation behind Fig. 2 and Figs. 6 to 8, peaking at 4.7 GB, and a
  gathering, peaking at 6.6. A whole `appendix` peaked at 7.2. Those are the
  highest we measured, and the column above adds headroom to them rather than
  reporting a floor we have tested. On
  Docker Desktop the VM carries its own limit, under Settings, Resources, and
  its default may be below these.
- Figures are written to `results/figures/` as PDF. The scripts open a window
  otherwise; setting `FIGDIR=` empty restores that.
- The container runs as root, so on Linux the files under `results/` come out
  owned by root. Add `--user "$(id -u):$(id -g)"` to the `docker run` to avoid
  that.
- On Linux, a `permission denied` on the Docker socket means your user is not
  in the `docker` group. Prefix the commands with `sudo`, or run
  `sudo usermod -aG docker "$USER"` and start a new login session, since the
  group change does not reach shells that are already open.

## Targets

One per claim of the camera-ready, and nothing else on that path. The times
here and below are what we measured, and may differ on other machines;
`REPORT.md` records what each target actually took.

| target | time | what |
| --- | --- | --- |
| `quick` | a few min | all three attacks; a smoke test, and it validates no claim on its own |
| `c1` | ~2 h 20 | **C1**, Sec. 3 BGV bootstrapping: Tab. 2 and Fig. 2 |
| `c2` | ~3 h | **C2**, Sec. 4 BFV pt-ct: Tab. 3 and Fig. 3 |
| `c3` | 10 min | **C3**, Sec. 5 BFV ct-ct: Fig. 4 and Fig. 5 |
| `report` | 1 s | collect `results/` into `REPORT.md` |
||

`c1` gathers its own failing ciphertexts (about two hours) and then ranks the
candidate positions of the key and draws Fig. 2 (about 20 minutes). Each
of the three also takes `AUTHORDATA=1`, which reads stored data instead of
running the attack: for `c2` and `c3` the logs we ship under `paper_results/`,
for `c1` a set of samples you gathered earlier:

| | time | |
| --- | --- | --- |
| `make c1 AUTHORDATA=1` | ~19 min | Tab. 2's analysis and Fig. 2, from a gathered set |
| `make c2 AUTHORDATA=1` | 1 s | Tab. 3 and Fig. 3, from our logs |
| `make c3 AUTHORDATA=1` | 2 s | Fig. 5, from our log |

Only a fresh run shows the attack works on your machine. For `c2` and `c3`,
`AUTHORDATA=1` meets the paper's own numbers to the digit. For `c1` it does not:
Table 2's counts belong to our sample, which we do not distribute, and any fresh
gathering is a different sample. `REPORT.md` says which it is showing.

### **The appendices.**

These belong to the full version of the paper rather than the camera-ready, and
we do not ask for them to be evaluated.

| target | time | what |
| --- | --- | --- |
| `appendix` | ~1 h 40 | the three below |
| `appendix-bgv` | ~36 min | Figs. 6, 7 and 8 |
| `appendix-ctct` | ~1 h | Tabs. 5 and 7, Fig. 11 |
| `appendix-figs` | 1 s | Figs. 9 and 11, from the stored logs |
||

- Every target draws its own figures and tables from its own output, so nothing
  else has to be run afterwards.
- `TRIALS=n` sets how many ct-ct trials the box plots use. The paper uses 40 and
  the cost is linear in it; the separation shows from about 10.
- `FAILS=n` sets how many failing ciphertexts `c1` gathers per parameter set.
  The setup of each gathering, which does not depend on `FAILS`, takes about
  26 minutes for A\* and 45 for B\*; a failing ciphertext after it costs about
  2.4 seconds for A\* and 4 for B\*.
- `THREADS=n` raises the OpenMP thread count, and the BFV attacks scale with it.
  The Sec. 3 analysis and the figure-drawing are Python and ignore it; the
  Sec. 3 collection threads part of its coefficient extraction, but spends its
  time in bootstrapping and in a setup that are not threaded, so one thread and
  eight come out within seconds of each other. The paper's timings are
  single-threaded, which is the default here.

## Contents

### **Section 3. BGV** 
Recall and precision per slice, from 50 to 500 ciphertexts in steps of 50.
By default `c1` gathers a fresh sample, so the reported counts may differ from
Tab. 2's. On the samples we gathered for the paper, `AUTHORDATA=1` gives the
same counts as Tab. 2.

What `c1` measures is support identification, not key recovery. It ranks the
coordinates by variance and keeps 97 candidates (32 per level for ℓ = 2, 4, 8,
and one for ℓ = 1); "correct" counts how many of the 32 nonzero positions of the
key are among them, graded against the key the gathering stored. The step that
turns the candidates into the key, described in Sec. 3, is not run here.

On a fresh sample, expect 30 to 32 of 32 at Tab. 2's counts, rising with the
number of ciphertexts; A\* can stay at 31/32 up to 500. What a fresh run should show is that trend: recall climbing towards
32/32 as ciphertexts are added, with the precision near 32/97. An attacker who
falls short cannot tell which position is missing, but can tell that the final
step failed, since a key is checked against the public key; it then gathers more
failures or keeps more candidates per level.

#### Samples

Section 3 runs on collected bootstrapping failures. `make c1` gathers them
itself (about two hours for both parameter sets at `THREADS=8`) and then
ranks the candidate positions of the key and draws Fig. 2 (about 20 minutes),
so nothing has to be downloaded. The gathered data lands
under `results/` as

```
data-collected-by-authors/m21845_p2_e8_ep4_r1/{coeff.npy, irwin_hall_21845.npy, secret_keys.txt}
data-collected-by-authors/m28679_p2_e8_ep4_r1/{coeff.npy, irwin_hall_28679.npy, secret_keys.txt}
```

If you already have a set gathered this way, mount its directory at
`/artifact/data-collected-by-authors` and the other targets read it with
`AUTHORDATA=1`.

A gathering writes the samples as text instead, which `np.loadtxt` reads far
more slowly and which takes three times the space. `to_npy.py` converts them
once, losslessly, and the analysis prefers the `.npy` whenever both are there;
`c1` does this for you, but it can be run on its own:

```bash
docker run --rm -v "$PWD/results:/artifact/results" \
  indcpad-artifact python3 to_npy.py results
```

### **Section 4. BFV pt-ct** 
Every run must end with

```
SK Coeff Match Ratio    : 1.000000
Final L2 Error (e)      : 0
```

which means that the key is fully recovered, from `N` decryption queries. 
Each run draws a fresh key, so the individual noise samples may differ from the paper; but the trend will remain the same. 

Figure 3 is drawn from the run's own log as `results/figures/fig3-01.pdf`. It
should look like the paper's: the measured noise (squares and diamonds) climbing
past the estimated `L^inf`-bound (triangles), through the certified thresholds.

Some keys stop the recovery before it starts. The window the decrypted result
has to land in is only `log2 t` bits wide, so the noise a key amplifies can run
past it, the decryption wraps, and nothing can be read back. Over 100 draws at
each ring dimension we saw this 13 times at `logN 12` and 29 at `logN 16`, the
two smallest plaintext moduli, and never at the other three. It shows at the
first index, within seconds, so such a draw costs nothing: `c2` simply draws
again, up to ten times, and `REPORT.md` says how many draws each dimension
took. A key that stops later is left alone, since rerunning it would cost
hours.

### **Section 5. BFV ct-ct** 
The distinguishing attack separates two distributions: the distributions `D_0` and `D_1` corresponding to the challenge bit $0$ and $1$. 
The Figure (analogous to Fig. 5) must show a clear separation, much larger deviation for $D_1$. 

- **A run may fail (with $\approx$ 9% failure probability). If it fails, run it again.** 
The gain `C = ρt` is fixed in advance, but the amplification scales as `k_s^K` for a freshly drawn key, where `k_s = max_ω |s(ω)| / √(2N/3)`. 
The window that separates the decrypted result is only `log2 t` bits wide: below `k_s ≈ 3.08` the noise stays under `q/2t`, and above `k_s ≈ 3.89` both challenges wrap around. 
Only in between, they separate. For our parameter sets, the band carries probability of `0.908`.
We have seen 8/10, 9/10 and 10/10 over three runs of `appendix-ctct`, as one can expect.
- The key behind Fig. 5 had `k_s = 3.305`, an ordinary draw against a median of 3.281.
- `ctct_multikey` does without any such band and reports what one fixed gain achieves over ten random keys (Tab. 5 and 7) 
It prints, per key, the largest `|y|` under `D_0` and the smallest `|y|` under `D_1`, with their ratio.
It marks a key as separating when those two ranges are disjoint.

## Licence

MIT, except for the files that modify HElib or OpenFHE, which stay under the
Apache 2.0 and BSD 2-Clause licences those libraries carry. `LICENSE` at the top
of the tree says which is which.
