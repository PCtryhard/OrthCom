#!/bin/bash
# Monomer validation of the designed Lock sequences (chain A) with ColabFold.
#
# WHY --num-recycle 3 (was 0):
#   --num-recycle 0 was a smoke test, not validation. AF2 refines its confidence
#   (pLDDT/pTM) iteratively across recycles; a single pass (0 recycles) gives a
#   systematically depressed, noisy signal. 3 recycles is the standard floor for a
#   usable confidence read.
#   NOTE on VRAM: recycling is applied in-place at inference (it is NOT unrolled),
#   so peak VRAM is set by sequence length x MSA depth x model -- NOT by recycle
#   count. Going 0 -> 3 does not materially raise peak memory for these ~172-res
#   single-sequence monomers; the old 0-recycle was a weak workaround.
#
# ROBUSTNESS (8 GB RTX 4060):
#   We fold ONE sequence per colabfold_batch process instead of pointing it at the
#   whole directory. A fresh process per sequence => clean GPU context, so an OOM on
#   one sequence releases all VRAM and the loop just continues (failures logged, not
#   fatal). Already-folded sequences are skipped, so re-running resumes in place.

set -u

# --- memory environment (PROPOSED -- see TF_FORCE_UNIFIED_MEMORY note) ----------
export XLA_PYTHON_CLIENT_PREALLOCATE=false   # don't grab ~90% VRAM up front
# Let XLA spill allocations to host RAM instead of a hard OOM. No cost when the
# fold fits in VRAM; only slower if it would otherwise have crashed. Remove this
# line for maximum speed if you never hit OOM.
export TF_FORCE_UNIFIED_MEMORY=1
# --------------------------------------------------------------------------------

SEQ_DIR="data/processed/mpnn/seqs"
OUT_DIR="data/processed/alphafold"
FAIL_LOG="${OUT_DIR}/failed_monomer.log"
mkdir -p "$OUT_DIR"
: > "$FAIL_LOG"

shopt -s nullglob
fastas=( "$SEQ_DIR"/*.fa )
if (( ${#fastas[@]} == 0 )); then
    echo "No .fa files in $SEQ_DIR (run src/run_mpnn.sh first)."
    exit 1
fi

for fa in "${fastas[@]}"; do
    name=$(basename "$fa" .fa)
    if compgen -G "${OUT_DIR}/${name}_unrelaxed_*.pdb" > /dev/null; then
        echo "[skip] $name already folded"
        continue
    fi
    echo "[fold] $name"
    localcolabfold/colabfold-conda/bin/colabfold_batch \
        "$fa" \
        "$OUT_DIR" \
        --num-recycle 3 \
        --num-models 1 \
        --model-type alphafold2_ptm \
        --msa-mode single_sequence \
    || { echo "$name" >> "$FAIL_LOG"; echo "[FAIL] $name (logged, continuing)"; }
done

echo "Done. Failures: $(wc -l < "$FAIL_LOG") (see $FAIL_LOG)"
