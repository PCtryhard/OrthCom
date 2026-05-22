#!/bin/bash

mkdir -p data/processed/mpnn

# Loop through every PDB file RFdiffusion generated
for file in data/processed/rfdiffusion/*.pdb; do
    echo "Running ProteinMPNN on $file..."
    python ProteinMPNN/protein_mpnn_run.py \
        --pdb_path "$file" \
        --pdb_path_chains "B" \
        --out_folder data/processed/mpnn \
        --num_seq_per_target 2 \
        --batch_size 1
done
