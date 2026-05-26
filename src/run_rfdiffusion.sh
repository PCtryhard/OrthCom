#!/bin/bash

mkdir -p data/processed/rfdiffusion

# Fire RFdiffusion with the /0 chain break to force a two-chain complex
python RFdiffusion/scripts/run_inference.py inference.output_prefix=data/processed/rfdiffusion/lock inference.input_pdb=data/raw/canonical_key.pdb inference.model_directory_path=models/weights 'contigmap.contigs=[A1-15/0 100-120]' 'ppi.hotspot_res=[A3,A8,A13]' inference.num_designs=100