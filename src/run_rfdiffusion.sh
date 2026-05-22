#!/bin/bash

# Create output directory if Kimi missed it
mkdir -p data/processed/rfdiffusion

# Fire RFdiffusion to hallucinate the Lock
python /path/to/RFdiffusion/scripts/run_inference.py \
    inference.output_prefix=data/processed/rfdiffusion/lock \
    inference.input_pdb=data/raw/synthetic_key.pdb \
    'contigmap.contigs=[A1-15/100-120]' \
    'ppi.hotspot_res=[A8]' \
    inference.num_designs=10
