#!/bin/bash

# Run ColabFold using the GPU for lightning-fast prediction
localcolabfold/colabfold-conda/bin/colabfold_batch \
    data/processed/mpnn \
    data/processed/alphafold \
    --num-recycle 3 \
    --model-type alphafold2_ptm
