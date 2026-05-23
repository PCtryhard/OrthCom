#!/bin/bash

# Force JAX to allocate memory dynamically to prevent WSL crashes
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=1

# Run ColabFold using the GPU
localcolabfold/colabfold-conda/bin/colabfold_batch \
    data/processed/mpnn/seqs \
    data/processed/alphafold \
    --num-recycle 3 \
    --model-type alphafold2_ptm
