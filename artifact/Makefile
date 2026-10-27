# Artifact for "New IND-CPA-D Attacks on BFV and BGV".
#
# One target per claim of the camera-ready, and nothing else on that path.
#
#   make quick     a few min All three attacks, at reduced parameters for
#                            Sections 4 and 5. Validates no claim on its own;
#                            it is there to show that the image works.
#
#   make c1        ~2 h 20   Section 3, BGV bootstrapping: Table 2 and
#                            Figure 2. Gathers its own failing ciphertexts.
#   make c2        ~3 h      Section 4, BFV pt-ct: Table 3 and Figure 3.
#   make c2-oracle ~50 min   The same attack, the secret key used only for
#                            key generation, the decryption oracle and grading.
#   make c3        ~10 min   Section 5, BFV ct-ct: Figures 4 and 5.
#
#   make report    instant   Sets what was produced beside what the paper says.
#
# Each of the three takes AUTHORDATA=1, which reads stored data instead of
# running the attack: our shipped logs for c2 and c3, a set gathered earlier by
# `make c1` for c1. For c2 and c3 that is the fast way to see that the paper's
# own numbers are the ones its runs produced:
#
#   make c1 AUTHORDATA=1   ~19 min      make c3 AUTHORDATA=1   2 s
#   make c2 AUTHORDATA=1   1 s
#
# The appendices are not part of the evaluation. They belong to the full version
# of the paper rather than the camera-ready, and are here for completeness:
#
#   make appendix       ~1 h 40   the three below
#   make appendix-bgv   ~36 min   Figures 6, 7 and 8
#   make appendix-ctct  ~1 h      Tables 5 and 7, Figure 11
#   make appendix-figs  ~1 s      Figures 9 and 11, from the stored logs
#
# The times are what we measured, and may differ on other machines;
# `make report` records what each target actually took.
#
# THREADS sets the OpenMP thread count and defaults to 1, because the timings
# the paper quotes are single-threaded. Raising it changes the wall time and not
# the results; 8 is a good choice, and scaling saturates near there. Section 3
# barely moves with it, its time going to bootstrapping and to a setup that are
# not threaded.
#
# TRIALS sets how many ct-ct trials the box plots are built from. The paper uses
# 40; the cost is linear in it and the separation shows from about 10 onward.
#
# The BGV analysis memory-maps the samples and stays near 3 GB whatever
# their size. The density estimation behind the figures peaked at 4.7 GB and a
# gathering at 6.6, so give the machine 10. Docker Desktop's own VM limit,
# under Settings, Resources, may sit below that.

THREADS ?= 1
TRIALS   ?= 40
RESULTS  ?= results
DATA     ?= data-collected-by-authors

# How many failing ciphertexts c1 gathers per parameter set. The setup of each
# collection, which does not depend on this, takes about 26 minutes for A* and
# 45 for B*; after it, a failing ciphertext costs about 2.4 seconds for A* and 4
# for B*, and 6 GB is written per set. The paper stopped
# at 492 for A*. The analysis sweeps 50 to 500, so that the trend shows past
# Table 2's counts on a fresh sample.
FAILS ?= 500

# Hamming weight of the ordinary secret key in c1 (HElib's t); see the README.
HWT ?= 32

CTCT13 := OMP_NUM_THREADS=$(THREADS) indcpad_ctct
CTCT15 := OMP_NUM_THREADS=$(THREADS) /opt/openfhe-v15/bin/indcpad_ctct
MK13   := OMP_NUM_THREADS=$(THREADS) /opt/openfhe-v13/bin/ctct_multikey
MK15   := OMP_NUM_THREADS=$(THREADS) /opt/openfhe-v15/bin/ctct_multikey
PTCT   := OMP_NUM_THREADS=$(THREADS) bfv_ptct_sk_general
PTCTO  := OMP_NUM_THREADS=$(THREADS) bfv_ptct_oracle
BGV    := indcpad
PLOT   := python3 /artifact/BGV/BGV-BTS-KRD-slices-and-plot.py

# Where Section 3 reads from: what c1 has just gathered, or what we shipped.
ifeq ($(AUTHORDATA),)
  BGVSRC := $(RESULTS)
else
  BGVSRC := $(DATA)
endif

# Every target stamps its own start and end, so that the report can say how long
# each step took without the evaluator having to wrap anything in `time`.
tick = @printf '%s\t%s\t%s\n' '$@' "$$(date +%s)" '$(1)' >> $(RESULTS)/timings.tsv

# And records whether what it produced came from a run here or from stored logs
# and samples, so that the report can say which it is showing.
mark = @printf '%s\t%s\n' '$@' '$(1)' >> $(RESULTS)/sources.tsv

.PHONY: help quick c1 c2 c2-oracle c3 report clean \
        appendix appendix-bgv appendix-ctct appendix-figs

help:
	@sed -n '1,/^$$/p' Makefile | sed 's/^# \{0,1\}//'
	@echo "Current: THREADS=$(THREADS) TRIALS=$(TRIALS) RESULTS=$(RESULTS) DATA=$(DATA) HWT=$(HWT)"

## ------------------------------------------------------------------- quick --
# The smoke test. Section 3 runs at the paper's parameters, on a set gathered
# earlier by `make c1`; Sections 4 and 5 run fresh at reduced parameters, which
# is enough to show the attacks work.
quick: | $(RESULTS)
	$(call tick,start)
	@test -d $(DATA)/m21845_p2_e8_ep4_r1 -o -d $(DATA)/Collected_data/m21845_p2_e8_ep4_r1 \
	  || { echo "$(DATA) holds no samples: mount them, or run 'make c1'"; exit 1; }
	$(call mark,authors)
	@echo "== Section 3, from the samples under $(DATA)"
	for m in 21845 28679; do \
	  cd /artifact/$(DATA) && $(PLOT) --m $$m --slices 50:510:50 \
	    > /artifact/$(RESULTS)/c1-keyrecovery-m$$m.txt 2>&1 ; \
	done
	@echo "-- A* (m=21845):" && sed -n '/Ctxts *Rows/,$$p' $(RESULTS)/c1-keyrecovery-m21845.txt
	@echo "-- B* (m=28679):" && sed -n '/Ctxts *Rows/,$$p' $(RESULTS)/c1-keyrecovery-m28679.txt
	@echo "== Section 4, at the smallest ring dimension"
	log=$(RESULTS)/general_logN12_ptm786433_depth4.txt; \
	for draw in 1 2 3 4 5 6 7 8 9 10; do \
	  $(PTCT) --logN 12 --ptm 786433 --depth 4 --cF 1.0 > $$log 2>&1 ; \
	  if grep -q "SK Coeff Match Ratio" $$log; then break; fi; \
	  if ! grep -qE "idx [1-5]: Cannot increase" $$log; then break; fi; \
	  echo "   this key was one the parameters do not cover; drawing another"; \
	done; \
	grep -E "SK Coeff Match Ratio|Final L2 Error|O_Dec calls" $$log
	@echo "== Section 5, five key pairs at a reduced ring dimension"
	$(MK13) 5 2 786433 35 6.0 3800000000000 2.83 4096 262144 0 \
	  > $(RESULTS)/quick-ctct.txt 2>&1
	@tail -2 $(RESULTS)/quick-ctct.txt
	$(call tick,end)
	@$(MAKE) --no-print-directory report
	@echo "== all three attacks ran; see $(RESULTS)/REPORT.md"

## ---------------------------------------------------------------------- C1 --
# BGV bootstrapping. The failures are gathered here rather than read from ours,
# which takes about two hours and writes 12 GB; ranking the candidate positions
# of the key and drawing Figure 2 take about 20 minutes more. AUTHORDATA=1 reads
# a set gathered earlier instead. Table 2's published counts belong to our
# sample, which is not distributed, so a fresh gathering reproduces the trend,
# not the counts. The
# slices run from 50 to 500 in steps of 50, which covers the ones Table 2 prints:
# 50, 150 and 250 ciphertexts for A*, 50, 200 and 350 for B*.
#
# What the analysis measures is support identification: how many of the 32
# nonzero positions of the key fall among its 97 candidates (32 per level for
# l = 2, 4, 8, one for l = 1). The step that turns the candidates into the key
# is not run here.
#
# Of the plots the analysis can draw, only the A* coefficient distributions are
# Figure 2. The rest belong to the full version; see appendix-bgv.
c1: | $(RESULTS)
	$(call tick,start)
ifeq ($(AUTHORDATA),)
	$(call mark,run)
	@echo "== gathering the failing ciphertexts, both parameter sets"
	@echo "   about 46 minutes for A* and 72 for B*, at THREADS=8"
	for phim in 16384 23040; do \
	  cd /artifact/$(RESULTS) && $(BGV) p=2 r=1 c=3 L=400 N=$$phim t=$(HWT) \
	    nthreads=$(THREADS) skHwt_sp=32 iter=100000 waitForInput=0 noPrint=1 \
	    progress=1 progressStep=10 stopAfterFailCount=$(FAILS) \
	    > c1-collect-N$$phim.txt 2>&1 ; \
	done
	@echo "== converting them, which the analysis then memory-maps"
	python3 to_npy.py $(RESULTS)
else
	@test -d $(DATA)/m21845_p2_e8_ep4_r1 -o -d $(DATA)/Collected_data/m21845_p2_e8_ep4_r1 \
	  || { echo "$(DATA) holds no samples: mount them, or drop AUTHORDATA"; exit 1; }
	$(call mark,authors)
	@echo "== reading the samples under $(DATA)"
endif
	@echo "== Table 2, ranking the candidate positions of the key"
	for m in 21845 28679; do \
	  cd /artifact/$(BGVSRC) && $(PLOT) --m $$m --slices 50:510:50 \
	    > /artifact/$(RESULTS)/c1-keyrecovery-m$$m.txt 2>&1 ; \
	done
	@echo "-- A* (m=21845):" && sed -n '/Ctxts *Rows/,$$p' $(RESULTS)/c1-keyrecovery-m21845.txt
	@echo "-- B* (m=28679):" && sed -n '/Ctxts *Rows/,$$p' $(RESULTS)/c1-keyrecovery-m28679.txt
	@echo "== Figure 2, the A* coefficient distributions; about 18 minutes"
	cd /artifact/$(BGVSRC) && FIGNAME=fig2 $(PLOT) --m 21845 --plot-coeff-dist
	$(call tick,end)

## ---------------------------------------------------------------------- C2 --
# BFV pt-ct, at the five ring dimensions of Table 3. Look for
# "SK Coeff Match Ratio : 1.000000" in each: that is the claim. Two of the three
# hours are the largest ring dimension.
#
# Some keys stop the recovery before it starts: the window the decrypted result
# has to land in is log2 t bits wide, and the noise a key amplifies can run past
# it, so the decryption wraps and nothing can be read back. Over 100 draws each
# we saw it 13 times at logN 12 and 29 at logN 16, the two smallest plaintext
# moduli, and never at 13, 14 or 15. It shows at the first index within seconds,
# so such a draw costs nothing and this simply draws again, up to ten times. A
# key that stops later is left alone, since rerunning it would cost hours.
c2: | $(RESULTS)
	$(call tick,start)
ifeq ($(AUTHORDATA),)
	$(call mark,run)
	@: > $(RESULTS)/c2-draws.txt
	for case in "12 786433" "13 4194305" "14 67108865" "15 33423361" "16 1048577"; do \
	  set -- $$case; \
	  log=$(RESULTS)/general_logN$$1_ptm$$2_depth4.txt; \
	  for draw in 1 2 3 4 5 6 7 8 9 10; do \
	    echo "-- logN $$1, ptm $$2 (draw $$draw)"; \
	    $(PTCT) --logN $$1 --ptm $$2 --depth 4 --cF 1.0 > $$log 2>&1 ; \
	    if grep -q "SK Coeff Match Ratio" $$log; then break; fi; \
	    if ! grep -qE "idx [1-5]: Cannot increase" $$log; then break; fi; \
	    echo "   this key was one the parameters do not cover; drawing another"; \
	  done; \
	  echo "logN $$1	$$draw" >> $(RESULTS)/c2-draws.txt; \
	  grep -E "SK Coeff Match Ratio|Final L2 Error" $$log \
	    || echo "   !! logN $$1 did not recover; see the log"; \
	done
else
	$(call mark,authors)
	@echo "== reading the logs behind the paper's own Table 3"
	cp BFV-pt-ct/paper_results/general_logN*_depth4.txt $(RESULTS)/
	@grep -H -E "SK Coeff Match Ratio" $(RESULTS)/general_logN*_depth4.txt
endif
	@echo "-- Table 3"
	cd BFV-pt-ct && python3 -c "import log2latex; log2latex.generate_table('/artifact/$(RESULTS)')" \
	  | tee /artifact/$(RESULTS)/table3.tex | tail -8
	@echo "-- Figure 3, drawn as a PDF from the run's own log"
	cd BFV-pt-ct && FIGNAME=fig3 python3 log2latex.py \
	  /artifact/$(RESULTS)/general_logN16_ptm1048577_depth4.txt --plot
	$(call tick,end)

# The same attack as Algorithm 1 states it: K = d_leg multiplications and one
# decryption query per root, then a check of the recovered key against the
# public key. The secret key is used only for key generation, the decryption
# oracle and grading the match ratio. c2 runs the version we tuned the attack
# with, which also reads the noise with the secret key; on Table 3's parameters
# both make the same queries.
c2-oracle: | $(RESULTS)
	$(call tick,start)
	for case in "12 786433" "13 4194305" "14 67108865" "15 33423361" "16 1048577"; do \
	  set -- $$case; \
	  log=$(RESULTS)/oracle_logN$$1_ptm$$2_depth4.txt; \
	  echo "-- logN $$1, ptm $$2"; \
	  $(PTCTO) --logN $$1 --ptm $$2 --depth 4 --cF 1.0 > $$log 2>&1 ; \
	  grep -E "Multiplications|O_Dec calls|Public-key check|SK Coeff Match Ratio" $$log; \
	done
	$(call tick,end)

## ---------------------------------------------------------------------- C3 --
# BFV ct-ct at the paper's parameters, on the OpenFHE that predates the revised
# bound. This attack stops where the two challenge distributions separate, which
# is what breaks indistinguishability, so it recovers no key and costs minutes
# rather than hours.
#
# A run comes out flat about one time in eleven, on a key the fixed gain does not
# suit. That is the attack rather than a fault: the gain is fixed before the key
# is drawn, as an adversary must fix it. Run it again.
c3: | $(RESULTS)
	$(call tick,start)
ifeq ($(AUTHORDATA),)
	$(call mark,run)
	$(CTCT13) 786433 35 4.5 262144 24800 $(TRIALS) 0 1e9 2.0 \
	  > $(RESULTS)/c3-ctct.txt 2>&1
else
	$(call mark,authors)
	@echo "== reading the log behind the paper's own Figures 4 and 5"
	cp BFV-ct-ct/paper_results/boxplot-4nt-40.txt $(RESULTS)/c3-ctct.txt
endif
	@tail -4 $(RESULTS)/c3-ctct.txt
	@echo "-- Figure 4 is the per-squaring table in $(RESULTS)/c3-ctct.txt"
	@echo "-- Figure 5, drawn from it"
	cd BFV-ct-ct && FIGNAME=fig5 python3 ctct-boxplot.py /artifact/$(RESULTS)/c3-ctct.txt
	$(call tick,end)

## ---------------------------------------------------------------- appendix --
# NONE OF THIS IS IN THE CAMERA-READY, which is capped at 15 pages; it belongs to
# the full version of the paper, and we do not ask for it to be evaluated.
appendix: appendix-bgv appendix-ctct appendix-figs
	@echo "== appendix experiments done; see $(RESULTS)"

# Figure 6 is the B* counterpart of Figure 2, and Figures 7 and 8 the variances.
appendix-bgv: | $(RESULTS)
	$(call tick,start)
	@test -d $(BGVSRC)/m21845_p2_e8_ep4_r1 -o -d $(BGVSRC)/Collected_data/m21845_p2_e8_ep4_r1 \
	  || { echo "   skipped: $(BGVSRC) holds no samples"; exit 0; }
	@echo "   36 minutes; elapsed time is printed as it goes"
	cd /artifact/$(BGVSRC) && FIGNAME=fig6 $(PLOT) --m 28679 --plot-coeff-dist
	for m in 21845 28679; do \
	  cd /artifact/$(BGVSRC) && FIGNAME=fig7-fig8-m$$m $(PLOT) --m $$m --plot-vars ; \
	done
	$(call tick,end)

# Table 5, ten key pairs under the bound the body assumes; then Table 7 and
# Figure 11, the same attack against the revised bound of the current OpenFHE.
appendix-ctct: | $(RESULTS)
	$(call tick,start)
	$(MK13) 10 3 786433 35 4.5 24800 3.36 0 262144 0 \
	  > $(RESULTS)/appendix-multikey-4nt.txt 2>&1
	@tail -2 $(RESULTS)/appendix-multikey-4nt.txt
	$(CTCT15) 1153 45 8.15 262144 8737706 $(TRIALS) 3.40 3.48 \
	  > $(RESULTS)/appendix-boxplot-8nt.txt 2>&1
	$(CTCT15) 257 48 7.93 262144 189000000000 10 0 1e9 \
	  > $(RESULTS)/appendix-boxplot-8nt-t257.txt 2>&1
	$(MK15) 10 10 1153 45 8.15 8737706 3.374 0 262144 0 \
	  > $(RESULTS)/appendix-multikey-8nt.txt 2>&1
	@tail -2 $(RESULTS)/appendix-multikey-8nt.txt
	@echo "-- Figure 11, from the run just made"
	cd BFV-ct-ct && FIGNAME=fig11 python3 ctct-boxplot.py /artifact/$(RESULTS)/appendix-boxplot-8nt.txt
	$(call tick,end)

# Redrawn from the logs we ship, so that the paper's own figures can be checked
# against the runs that produced them.
appendix-figs: | $(RESULTS)
	$(call tick,start)
	cd BFV-ct-ct && FIGNAME=fig11-from-stored-log \
	  python3 ctct-boxplot.py paper_results/v150/boxplot-8nt-40.txt
	cd BFV-pt-ct && python3 log2latex.py \
	  paper_results/custom_logN14_ptm786433_logC3_d2.txt \
	  > /artifact/$(RESULTS)/fig9-logN14.tex
	cd BFV-pt-ct && python3 log2latex.py \
	  paper_results/general_logN15_ptm33423361_depth4.txt \
	  > /artifact/$(RESULTS)/fig9-logN15.tex
	@echo "-- wrote $(RESULTS)/fig9-logN{14,15}.tex"
	$(call tick,end)

## ------------------------------------------------------------------ report --
# Collects whatever is in $(RESULTS) into one table, next to the paper's values.
report: | $(RESULTS)
	@python3 /artifact/report.py $(RESULTS)

$(RESULTS):
	mkdir -p $(RESULTS)

clean:
	rm -rf $(RESULTS)
