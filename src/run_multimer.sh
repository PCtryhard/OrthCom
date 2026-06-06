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
# PERFORMANCE / STABILITY (8 GB RTX 4060 under WSL2) -- learned the hard way:
#   The real bottleneck was NOT VRAM size but JAX's default BFC pool allocator,
#   which on this WSL2/driver combo fragmented catastrophically: ~5 min PER
#   recycle, then a hard segfault at recycle=2. Switching to the platform
#   allocator (raw cudaMalloc/cudaFree, XLA_PYTHON_CLIENT_ALLOCATOR below) fixed
#   both -- stable runs at ~1 s/recycle after the first compiled pass.
#   TF_FORCE_UNIFIED_MEMORY (a TensorFlow flag) did nothing here: colabfold infers
#   with JAX/XLA, not the TF GPU runtime, so it is dropped.
#   We still fold ONE complex per process: a failure frees all VRAM and the loop
#   continues (failures logged, not fatal); completed jobs are skipped on re-run.
#
# RECYCLES:
#   Use alphafold2_multimer_v3's intended schedule -- up to 20 recycles with early
#   stop at 0.5 A. Fast-converging complexes stop early; with recycles ~1 s each
#   the cap is nearly free, and test metrics were still climbing past recycle 3.

set -u

export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

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
        --num-recycle 20 \
        --recycle-early-stop-tolerance 0.5 \
        --num-models 1 \
        --model-type alphafold2_multimer_v3 \
        --msa-mode single_sequence \
    || { echo "$name" >> "$FAIL_LOG"; echo "[FAIL] $name (logged, continuing)"; }
done

echo "Done. Failures: $(wc -l < "$FAIL_LOG") (see $FAIL_LOG)"
echo "Rank with: python src/parse_interface.py"
