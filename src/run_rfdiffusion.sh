#!/bin/bash

mkdir -p data/processed/rfdiffusion

# Fire RFdiffusion with the Iron Maiden constraints and a massive 150-200 length limit
python RFdiffusion/scripts/run_inference.py inference.output_prefix=data/processed/rfdiffusion/lock inference.input_pdb=data/raw/CANON.pdb inference.model_directory_path=models/weights 'contigmap.contigs=[A1-15/0 150-200]' 'ppi.hotspot_res=[A1,A4,A8,A12,A15]' inference.num_designs=1000