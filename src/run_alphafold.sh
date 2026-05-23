#!/bin/bash

# Update the input path to wherever the actual files are!
localcolabfold/colabfold-conda/bin/colabfold_batch \
    data/processed/mpnn/seqs \
    data/processed/alphafold \
    --num-recycle 3 \
    --model-type alphafold2_ptm
