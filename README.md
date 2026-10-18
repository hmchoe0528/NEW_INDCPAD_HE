# New IND-CPA-D Attacks on BFV and BGV
This repository is the implementation of the following paper:

> **New IND-CPA-D Attacks on BFV and BGV**  
> Hyeongmin Choe, Minsik Kang, and Damien Stehlé  
> ACM Conference on Computer and Communications Security (CCS) 2026, https://doi.org/10.1145/3830454.3846665

It includes the following CPA-D-style attacks, one per section of the paper: 

  - Section 3. [BGV Bootstrapping Attack using HElib](#bgv-bootstrapping-attack-using-helib) 
    - Related to Tables 1 and 2, and Figures 2, 6, 7, and 8. 
  - Section 4. [BFV PT-CT Attack using OpenFHE](#bfv-pt-ct-attack-using-openfhe)
    - Related to Tables 3, and Figures 3 and 9. 
  - Section 5. [BFV CT-CT Attack using OpenFHE](#bfv-ct-ct-attack-using-openfhe)
    - Related to Tables 4, 5, 6, and 7, and Figures 4, 5, 10, and 11.
    
# Reproducing the Paper Results
The sections below build each attack by hand against upstream HElib and OpenFHE. 

For a one-command, containerized reproduction of every figure and table — no manual build, no library download — use the Docker artifact in [`artifact/`](artifact/README.md): `make quick` for a two-minute smoke test, and `make c1`, `make c2`, `make c3` for the three claims. 

See [`artifact/README.md`](artifact/README.md) for further details.


# BGV Bootstrapping Attack using HElib
See `BGV` directory. 
It includes the following codes, for the BGV bootstrapping attack using [HElib library](https://github.com/homenc/HElib), as described in Section 3. 

  - (a) A package of KR-D attack scripts, mounting bootstrappings, then decrypting, and collecting the intermediate values from the failed slots, which can be publicly back-tracked. It also stores the statistical mapping depending only on $m$, and the secret keys to grade the KR-D attack's output. It works for the parameters A&#42; and B&#42; described in Table 1. The time taken for collecting the intermediate values are reported in Table 2. 

  - (b) A script `BGV-BTS-KRD-slices-and-plot.py` that analyzes the collected samples to generate the distribution plots (such as Figures 2 and 6 for the reconstructed coefficients; and Figures 7 and 8 for the variances with the corresponding non-zero indices of sk.), launching a KR-D attack to recover the secret key from the failed data (by using certain amount of slices among the whole number of slots). The scores of the KR-D attack for each size of the slices are also reported in Table 2. 

## How to use
### (a) Collecting Bootstrapping Failures
First, download HElib from https://github.com/homenc/HElib (tested using the commit `3e337a6` from July 18th, 2023). 
Then copy the modified files from this repository into the corresponding paths in your downloaded HElib's directory, then overwrite existing files. Required replacement files are:
- In `include/helib/`: `Context.h`, `Ctxt.h`, `keys.h`, `sample.h`, `debugging.h`
- In `src/`: `Context.cpp`, `Ctxt.cpp`, `debugging.cpp`, `keys.cpp`, `recryption.cpp`, `sample.cpp`
- In `tests/`: `Algebra.hpp`, `CMakeLists.txt`, `IND_CPAD_attack.cpp`

Then build, for e.g., 

```bash
mkdir build
cd build
cmake \
  -DGMP_DIR=/usr/local \
  -DNTL_DIR=/usr/local \
  -DENABLE_TEST=ON \
  -DPEDANTIC_BUILD=OFF \
  -DENABLE_THREADS=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  ..
make -j indcpad
```

`ENABLE_THREADS` requires an NTL built with `NTL_THREADS=on`, which the Ubuntu and Homebrew packages are; without it HElib ignores the `nthreads` argument below and collects on a single core. You may need to adjust `GMP_DIR`/`NTL_DIR` for your system. The executable will be generated at `build/bin/indcpad`. 
For other build errors, please see HElib's instruction. 

Then, you can now run the attack executable. A quick test will be:

```bash
/path/to/HElib/build/bin/indcpad \
  p=2 r=1 c=3 L=400 N=0 t=32 nthreads=32 skHwt_sp=32 iter=1000 \
  waitForInput=0 progress=1 progressStep=10 stopAfterFailCount=1
```

You can run until 1000 failing ciphertexts are collected, taking several hours to days, depending on `nthreads` (though, it may not scale great):

```bash
/path/to/HElib/build/bin/indcpad \
  p=2 r=1 c=3 L=400 N=0 t=32 nthreads=32 skHwt_sp=32 iter=100000 \
  waitForInput=0 noPrint=1 progress=1 progressStep=10 stopAfterFailCount=1000
```

Here, `iter` is the maximum number of attack iterations, while `stopAfterFailCount` is the number of failing iterations to collect before stopping.
The output files will be saved under `Collected_data` in the current working directory where you run the executable.

The output files run to over 6 GB per parameter set before compression, and a full collection takes around a day (or slightly more). A shorter run of a few hundred failures already suffices to recover the key, as reported in Table 2.

### (b) KR-D Recovery
Now, we are ready to run the KR-D attack with the collected outputs.
Appropriately update the parameters/paths, then run with these possible options. 

Try 
```bash
python BGV-BTS-KRD-slices-and-plot.py [--plot-vars | --plot-coeff-dist]
```

The default option without any flag analyzes the data with its slices (for Table 2). 
Among the data point from the slices, it chooses the top candidate indices by variance per $\ell$. 

The plot option `--plot-vars` plots a 2×2 histogram grid (for Figures 7 and 8) showing the variance distribution of all $m$ indices, one panel per $\ell$. The actual secret key indices are highlighted with red markers for visual inspection.

The other plot option `--plot-coeff-dist` plots 4 figures (for Figures 2 and 6) showing the distributions of $$c_j$$ per $\ell$, distinct graphs for its corresponding secret key value $$s_j \in \{-1, 0, 1\}$$. 

<!------------------------------
--------------------------------
------------------------------->

# BFV PT-CT Attack using OpenFHE
See `BFV-pt-ct` directory. 
It includes the following codes, for the BFV pt-ct mult attack using [OpenFHE library](https://github.com/openfheorg/openfhe-development), as described in Section 4. 

  - (a) A package of KR-D attack scripts, mounting pt-ct multiplications, then decryptions, collecting the decrypted values that failed to correctly decrypted. It then applies the KR-D attack, fully recovering the secret keys. 
  It works for some selected parameters, including those presented in Table 3. The time taken for pt-ct multiplications and recovering the secret keys are reported. It also generates raw data used to plot Figures 3 and 9. 
  
  - (b) A package of additional scripts for finding easy-to-break parameters are also included. Note that this is specific to our script, but not the attack, meaning that by modifying the attack scripts, you can mount the attack to other parameters, as well. 

  - (c) An additional script that draws the plots from the raw output of (a), such as Figures 3 or 9. 

## How to use
### (a) PT-CT Multiplications and KR-D attack
First, download [OpenFHE library](https://github.com/openfheorg/openfhe-development) (tested using the commit `df495ba` from February 26th, 2026). 

- In `openfhe-development/src/pke/examples`: Put `bfv_ptct.hpp`, and `bfv_ptct_sk_general.cpp` beating error estimations for any $pt \in \mathcal{R}_t$, i.e., 

$$\log_2 (Q / 2t) > K \cdot \log_2(\delta_N \cdot t / 2) +\log_2(\textsf{Vnorm}).$$

  - (additional) `bfv_ptct_sk.cpp` beating circuit-specific error estimations w.r.t.

  $$\Vert pt \Vert_\infty:$$

  $$\log_2 (Q / 2t) > K \cdot \log_2(\delta_N \cdot \Vert pt \Vert_\infty) + \log_2(\textsf{Vnorm}).$$

- In `openfhe-development/src/pke/`: Replace `CMakeLists.txt` by `(pke) CMakeLists.txt`

Both scripts recover the whole sk within $O(N)$ decryption queries and in several hours, at most. 

Optionally, for fixed seeds, 
- in the main directory: Replace `CMakeLists.txt` by `(main) CMakeLists.txt`,
- in `src/core/lib/utils/prng/`: Replace `blake2engine.cpp`, 
and configure with `-DWITH_FIXED_SEED=ON`.

Now try,

```bash
cd ~/openfhe-development 
mkdir build 
cmake -S . -B build
cmake --build build -j
./build/bin/examples/pke/bfv_ptct_sk_general --logN 12 --ptm 786433 --depth 4
```
For build errors, please see OpenFHE's instruction, which usually solves most of the problems. 

Note that, the results in Table 3 and Figures 3 and 9 can be fully reproduced with a fixed seed.
Try, 

```bash
cmake -S . -B build -DWITH_FIXED_SEED=ON
cmake --build build -j
OMP_NUM_THREADS=1 ./build/bin/examples/pke/bfv_ptct_sk_general --logN 12 --ptm 786433 --depth 4
```

One can still randomize by running without the "OMP_NUM_THREADS=1" option, possibly with multiple threads. 
The attack results can be found in `BFV-pt-ct/paper_results` directory. 

To re-run all the parameters in Table 3, try running `run_attacks.sh` in the main directory:

```bash
./run_attacks.sh
```
It will take several hours. 

### (b) Finding Parameters
We used `bfv_ptct_param_prescreener.cpp` to find appropriate parameters that supports the full sk recovery (i.e., with absolute error for $$e_0(\omega_i)$$ to be $<0.5$ for all $i = 1, \dots, N$) for `bfv_ptct_sk_general.cpp`. 
It scans candidate tuples $$(\log N, t, \text{SetDepth})$$ and keeps those for which the attack scale $$c := \Vert pt \Vert_\infty$$ lies in a feasible range satisfying $$Q / 2t < \Vert e_K\Vert_\infty \lesssim Q/2$$, or its empirical variant:

$$\frac{Q}{2t} \cdot 2^{m_{\text{low}}} < \Vert e_K\Vert_\infty \left( \approx \Vert e_1\Vert_\infty \cdot (c \cdot N / \pi)^{K-1} \lesssim (c \cdot E_1(\log N)) \cdot (c \cdot N / \pi)^{K-1} \right) \leq \frac{Q}{2} \cdot 2^{-m_{\text{up}}},$$

where $$E_1(\log N)$$ is the profiled value of $$\log_2 \Vert e_1\Vert_\infty$$ for each ring size, experimentally through `bfv_e1_inf_profile.cpp` (see `BFV-pt-ct/paper_results/e1_inf` for its output), and $$m_{\text{low}}, m_{\text{up}}$$ be lower-/upper-margin bits. 
For e.g., one can try, 

```bash
./build/bin/examples/pke/bfv_ptct_param_prescreener --logN 15 --top 10 --lower-margin 0 --upper-margin 1
```

Note that, the outputs are the candidate parameters, but may not always work with the specific hyper-parameters fixed in the attack codes. 

### (c) Plot
Use `log2latex.py` with output files of (a).

```bash
python3 log2latex.py "file_name"
```


<!------------------------------
--------------------------------
------------------------------->

# BFV CT-CT Attack using OpenFHE
See `BFV-ct-ct` directory. 
It includes the following codes, for the BFV ct-ct mult attack using [OpenFHE library](https://github.com/openfheorg/openfhe-development), as described in Section 5. 

  - (a) A package of IND-CPA-D attack scripts, (phase 1) mounting $K$ ct-ct multiplications, then decryptions, to calibrate

$$(\omega, k_a, K),$$

(phases 2 and 3) building reference distributions

$$\mathcal{D}_0,$$

and

$$\mathcal{D}_1,$$

from $$M$$ encryption queries on $$\mu_0 = 0$$ and $$\mu_1 = \mathcal{H}$$, respectively. The $$M$$ ciphertexts are combined with integer coefficients whose sum $$C$$ is a multiple of $ t $, so that the tracked plaintext is $$C\mu_\gamma = 0$$ and the decryption query is legitimate, while the noise becomes

$$\sum_i c_i \epsilon_i - \rho\cdot\nu(\mu_\gamma), \qquad \rho := C/t,$$

where $$\nu(\mu) := q\mu - t\Delta(\mu)$$ is the defect of the scaling encoding $$\Delta$$ (it vanishes only when $ t $ divides $ q $). The plaintext $$\mathcal{H}$$ is chosen to maximise $$\vert\nu(\mathcal{H})(\omega)\vert$$. The same $K$-fold amplification as in phase 1 is then applied and the result decrypted: it fails for both $$\gamma$$, and the two distributions are separated by the size of the returned coefficients, as shown in Figure 5. 
    - We note that it emulates the $$\vert a(\omega) \vert$$ enlarging procedures in phase 1 (which requires $$\ge 2^{20}$$ trials if truely applied), by using sk. However, the choice of $ a $ is still randomized: $ a $ is drawn uniformly from $$\mathcal{R}_q$$ and then shifted at the single coordinate $$\omega$$, so that $$\Vert a/q \Vert_\infty \approx 1/2$$ and all the canonical coordinates stay populated, as they would for a genuine sum of $$\mathsf{Enc}(0)$$'s. 
    - It likewise emulates the $ M $ encryption queries by drawing the aggregated encryption randomness directly, which is the same law up to the central limit theorem. Set the last argument of `EncryptLinComb` to `true` to perform the $ M $ queries instead.
    - Phase 3, which uses the challenge ciphertexts $$\mathsf{Enc}(\mu_\beta)$$ for the challenge bit $$\beta \in \{0, 1\}$$, is omitted: it is identical to one of the two runs of phase 2. One can simply distinguish by checking the magnitude of the decrypted coefficients. 

  - (b) A script `ctct-boxplot.py` that analyzes the collected decryption results to generate the distribution plots (such as Figure 5). 

  - (c) `ctct_multikey.cpp`, repeating phases 2 and 3 over independently generated key pairs, so that $$\mathsf{sk}$$, $$\mathsf{pk}$$ and every encryption randomness differ from run to run. 

## How to use
### (a) CT-CT Multiplications and IND-CPA-D attack
First, download [OpenFHE library](https://github.com/openfheorg/openfhe-development) (tested using the commit `2bf4524` from May 19th, 2025). We note that the depth-based parameter generation policy has been slightly changed from commit `5ebe2b2` (June 27th, 2025), newly introducing `deltaMS` making the attack less efficient.

- In `openfhe-development/src/pke/examples`: Put `bfv_ctct.hpp` and `indcpad_ctct.cpp`.
- In `openfhe-development/src/pke/`: Replace `CMakeLists.txt` by `(pke) CMakeLists.txt`

Build and run as in the previous section, e.g., 

```bash
cmake -S . -B build
cmake --build build -j
./build/bin/examples/pke/indcpad_ctct 786433 35 4.5 262144 24800 40 0 1e9 2.0
```

The arguments are explained in (d); without them the program runs the parameters of Appendix E.3 instead. It will take several hours; phase 1 alone takes minutes, and each of the $$2 \times 40$$ runs of phases 2 and 3 costs another $K$ ct-ct multiplications. 

You can find the results in the `BFV-ct-ct/paper_results` directory. `boxplot-4nt-40.txt` is the full log of the run reported in the paper: its first table is the phase-1 noise growth of Figures 4 and 10, its trial table gives Table 4, and `data-ctct-40x10.txt` is that trial table extracted for Figure 5. It was produced with the earlier library, as described in (d). 

### (b) Multiple key pairs
Build and run `ctct_multikey.cpp` the same way. It takes the number of key pairs and the number of runs per key pair first, then the parameters listed in (d):

```bash
./build/bin/examples/pke/ctct_multikey 10 3 786433 35 4.5 24800 3.36 0 262144 0
```

The result is stored in `paper_results/ct-ct-multikey-10x3`, and is reported as Table 5. It, too, uses the earlier library. Note that the absolute size of the decrypted coefficients varies a lot between key pairs, since the amplification factor is $$\prod_d 2t\,(a^*_d/q)(\omega)\,s(\omega)$$ and therefore depends on $$k_s = \vert s(\omega)\vert / \sqrt{2N/3}$$ raised to the power $K$; what is stable is the ratio between the two distributions, and $$\vert\epsilon_0(\omega)\vert$$ for $$\mu_\gamma = \mathcal{H}$$, which is deterministic and thus independent of the key.

### (c) Plot
Use `ctct-boxplot.py` with the output results of (a). It currently reads `paper_results/data-ctct-40x10.txt`, the trial table of the stored run; adjust the path for the figure of Appendix E.3, whose table is `paper_results/v150/data-ctct-8nt-40x10.txt`. 

```bash
python ctct-boxplot.py
```

Each row of that table is one run of phases 2 and 3, and each column pair is one coefficient index of the decrypted result: the left entry comes from $$\mu_\gamma = 0$$ and the right entry from $$\mu_\gamma = \mathcal{H}$$. The script collects them into the two distributions $$\mathcal{D}_0$$ and $$\mathcal{D}_1$$ and draws them side by side, which is Figure 5. Note that both are centred at $ 0 $, since the decrypted result is a sinusoid whose phase changes from run to run; what separates them is the spread, which is proportional to $$\vert\epsilon_0(\omega)\vert$$ and therefore carries the signature $$\nu(\mu_\gamma)$$. 

### (d) Mounting on OpenFHE v1.5.0
Since commit `5ebe2b2` the library bounds the rounding error $$I$$ by a separate extension factor $$\delta_{\mathsf{ms}} = 4\sqrt{N}$$ instead of $$\delta_N = 2\sqrt{N}$$, so the certified noise bound is $$(8Nt)^K \cdot V_N$$ rather than $$(4Nt)^K \cdot V_N$$, and the condition for the amplification to outpace it becomes $$k_a k_s > 12\sqrt{2}$$. Appendix E.3 reports the attack under that revised bound; Sections 4 and 5 use the earlier one.

`indcpad_ctct.cpp` takes, in order, `t`, `depth`, $$k_a$$, `M`, $$\rho$$, the number of trials, a lower and an upper bound on $$k_s$$ (the key is redrawn until it falls in that band; pass `0` and a large number to accept any key), and $$\delta_{\mathsf{ms}}/\sqrt{N}$$, which is `4` for OpenFHE since v1.4.2.4 and `2` before. That last argument is what lets one source serve both regimes. `ctct_multikey.cpp` takes the number of keys, the number of runs per key, then `t`, `depth`, $$k_a$$, $$\rho$$, a reference $$k_s$$, a ring dimension (`0` to let OpenFHE choose), `M`, and a flag selecting a fixed or a rescaled gain (`0` for fixed, as an adversary is bound to use).

The results of Sections 4 and 5 need the earlier library, which does not configure with current toolchains unless its bundled `google-benchmark` is disabled:

```bash
git -C openfhe-development worktree add ../openfhe-v13 2bf4524
cmake -S openfhe-v13 -B openfhe-v13/build -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_BENCHMARKS=OFF -DBUILD_UNITTESTS=OFF
make -C openfhe-v13/build indcpad_ctct ctct_multikey -j8
./openfhe-v13/build/bin/examples/pke/indcpad_ctct 786433 35 4.5 262144 24800 40 0 1e9 2.0
```

Those of Appendix E.3 need the current one, built as in (a):

```bash
./build/bin/examples/pke/indcpad_ctct  1153 45 8.15 262144 8737706 40 3.40 3.48
./build/bin/examples/pke/indcpad_ctct   257 48 7.93 262144 189000000000 10 0 1e9
./build/bin/examples/pke/ctct_multikey 10 10 1153 45 8.15 8737706 3.374 0 262144 0
```

The outputs are stored as follows.

| file | reported as |
| --- | --- |
| `paper_results/boxplot-4nt-40.txt` | Figures 4, 5 and 10, Table 4 |
| `paper_results/data-ctct-40x10.txt` | Figure 5 (trial table extracted from the above) |
| `paper_results/ct-ct-multikey-10x3` | Table 5 |
| `paper_results/v150/boxplot-8nt-40.txt` | Figure 11, and the $$t = 1153$$ figures of Appendix E.3 |
| `paper_results/v150/data-ctct-8nt-40x10.txt` | Figure 11 (trial table extracted from the above) |
| `paper_results/v150/boxplot-8nt-t257-10.txt` | the $$t = 257$$ figures of Appendix E.3 |
| `paper_results/v150/multikey-8nt-10x10.txt` | Table 7 |

Table 6 is analytic and needs no run. The multikey log was produced as three concurrent processes of 4, 3 and 3 keys, whose outputs were merged and the keys renumbered; running the command above in one process gives the same thing, more slowly.

Two points are worth knowing before re-running. First, the key is drawn afresh each time, and the gain $$C = \rho t$$ has to suit it: the amplification scales as $$k_s^K$$ while the readout window $$(q/2t,\, q/2)$$ is only $$\log_2 t$$ bits wide, so a $$\rho$$ tuned for one key can leave the noise below $$q/2t$$ or past $$q/2$$ for another. The band given to `indcpad_ctct` above keeps the figure reproducible; `ctct_multikey` deliberately does without it, and reports what a fixed $$C$$ achieves over random keys. Second, a run of forty trials takes some hours: the cost is $$2K$$ homomorphic squarings per trial, and is unaffected by `M`, which `EncryptLinComb` simulates rather than performing.


