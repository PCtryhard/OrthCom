#!/bin/bash

mkdir -p data/processed/mpnn

python ProteinMPNN/protein_mpnn_run.py \
    --pdb_dir data/processed/rfdiffusion \
    --out_folder data/processed/mpnn \
    --fixed_chain_label A \
    --num_seq_per_target 2 \
    --batch_size 1
