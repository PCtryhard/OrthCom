#!/bin/bash
# AF2-MULTIMER co-fold of Lock+Key complexes -- the real interface test.
#
# WHY a separate script from run_alphafold.sh:
#   The monomer validation uses alphafold2_ptm, which has no interface notion and
#   never emits iPTM. To score the Lock<->Key interface we MUST use a multimer
#   model: alphafold2_multimer_v3 (verified valid in colabfold/batch.py; it is also
#   colabfold's "auto" choice for complexes). Only multimer runs populate the
#   'iptm' key and a cross-chain PAE in *_scores_*.json.
#
# INPUT: data/processed/multimer_input/*.fa  (from src/prep_multimer.py). Each
#   record is  LockSeq:KeySeq  -- ':' is colabfold's chain-break delimiter.
#
# ROBUSTNESS (8 GB RTX 4060):
#   Multimer is heavier than ptm and the complex is ~187 res, so this is the run
#   most likely to hit the VRAM ceiling. We fold ONE complex per process: an OOM
#   frees all VRAM and the loop continues (failures logged, not fatal); completed
#   jobs are skipped on re-run. TF_FORCE_UNIFIED_MEMORY lets XLA spill to host RAM
#   instead of a hard OOM (proposed safety net -- remove for max speed).

set -u

export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=1

IN_DIR="data/processed/multimer_input"
OUT_DIR="data/processed/multimer"
FAIL_LOG="${OUT_DIR}/failed_multimer.log"
mkdir -p "$OUT_DIR"
: > "$FAIL_LOG"

shopt -s nullglob
fastas=( "$IN_DIR"/*.fa )
if (( ${#fastas[@]} == 0 )); then
    echo "No FASTAs in $IN_DIR."
    echo "Run: python src/prep_multimer.py --lock-seq-dir <your_top_N_dir>"
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
        --model-type alphafold2_multimer_v3 \
        --msa-mode single_sequence \
    || { echo "$name" >> "$FAIL_LOG"; echo "[FAIL] $name (logged, continuing)"; }
done

echo "Done. Failures: $(wc -l < "$FAIL_LOG") (see $FAIL_LOG)"
echo "Rank with: python src/parse_interface.py"
