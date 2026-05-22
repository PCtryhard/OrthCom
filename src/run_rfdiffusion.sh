#!/bin/bash

# Create output directory
mkdir -p data/processed/rfdiffusion

# Fire RFdiffusion to hallucinate the Lock
python RFdiffusion/scripts/run_inference.py \
    inference.output_prefix=data/processed/rfdiffusion/lock \
    inference.input_pdb=data/raw/synthetic_key.pdb \
    inference.model_directory_path=models/weights \
    'contigmap.contigs=[A1-15/100-120]' \
    'ppi.hotspot_res=[A8]' \
    inference.num_designs=10
