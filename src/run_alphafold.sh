#!/bin/bash

# Prevent JAX from swallowing all GPU VRAM, but remove the unified memory flag
# to prevent WSL from segfaulting during system RAM paging
export XLA_PYTHON_CLIENT_PREALLOCATE=false

# Run ColabFold using the GPU, restricted to 1 model for rapid screening
localcolabfold/colabfold-conda/bin/colabfold_batch \
    data/processed/mpnn/seqs \
    data/processed/alphafold \
    --num-recycle 3 \
    --num-models 1 \
    --model-type alphafold2_ptm
