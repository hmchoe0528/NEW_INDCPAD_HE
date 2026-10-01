#!/usr/bin/env zsh
set -e

BIN_GENERAL="./build/bin/examples/pke/bfv_ptct_sk_general"
BIN_CUSTOM="./build/bin/examples/pke/bfv_ptct_sk"
OUTDIR="./paper_results"
mkdir -p "$OUTDIR"

run_case() {
    local label="$1"
    local bin="$2"
    local args="$3"
    local outfile="$OUTDIR/${label}.txt"
    echo "=== Running: $label ===" | tee -a "$outfile"
    echo "Command: OMP_NUM_THREADS=1 $bin $args" | tee -a "$outfile"
    echo "Started: $(date)" | tee -a "$outfile"
    echo "" | tee -a "$outfile"
    OMP_NUM_THREADS=1 eval "$bin $args" 2>&1 | tee -a "$outfile"
    echo "" | tee -a "$outfile"
    echo "Finished: $(date)" | tee -a "$outfile"
    echo "=== Done: $label ===" 
    echo ""
}

# General runs (bfv_ptct_sk_gen)
run_case "general_logN12_ptm786433_depth4"    "$BIN_GENERAL" "--logN 12 --ptm 786433    --depth 4 --cF 1.0"
run_case "general_logN13_ptm4194305_depth4"   "$BIN_GENERAL" "--logN 13 --ptm 4194305   --depth 4 --cF 1.0"
run_case "general_logN14_ptm67108865_depth4"  "$BIN_GENERAL" "--logN 14 --ptm 67108865  --depth 4 --cF 1.0"
run_case "general_logN15_ptm33423361_depth4"  "$BIN_GENERAL" "--logN 15 --ptm 33423361  --depth 4 --cF 1.0"
run_case "general_logN16_ptm1048577_depth4"   "$BIN_GENERAL" "--logN 16 --ptm 1048577   --depth 4 --cF 1.0"

# Customized runs (bfv_ptct_sk)
# run_case "custom_logN12_ptm16777217_logC3_d2"  "$BIN_CUSTOM" "--logN 12 --ptm 16777217  --logC 3 --depth 2"
# run_case "custom_logN12_ptm16777217_logC8_d2"  "$BIN_CUSTOM" "--logN 12 --ptm 16777217  --logC 8 --depth 2"
# run_case "custom_logN13_ptm8388609_logC3_d2"   "$BIN_CUSTOM" "--logN 13 --ptm 8388609   --logC 3 --depth 2"
# run_case "custom_logN13_ptm8388609_logC4_d2"   "$BIN_CUSTOM" "--logN 13 --ptm 8388609   --logC 4 --depth 2"
# run_case "custom_logN14_ptm786433_logC3_d2"    "$BIN_CUSTOM" "--logN 14 --ptm 786433    --logC 3 --depth 2"
# run_case "custom_logN15_ptm4194305_logC2_d2"   "$BIN_CUSTOM" "--logN 15 --ptm 4194305   --logC 2 --depth 2"
# run_case "custom_logN16_ptm270794753_logC3_d2" "$BIN_CUSTOM" "--logN 16 --ptm 270794753 --logC 3 --depth 2"

echo "All runs complete. Results saved to $OUTDIR/"