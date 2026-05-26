#!/bin/bash

export XLA_PYTHON_CLIENT_PREALLOCATE=false

# Run ColabFold with 0 recycles to bypass the WSL memory crash
localcolabfold/colabfold-conda/bin/colabfold_batch \
    data/processed/mpnn/seqs \
    data/processed/alphafold \
    --num-recycle 0 \
    --num-models 1 \
    --model-type alphafold2_ptm \
    --msa-mode single_sequence
