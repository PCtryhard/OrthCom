#!/bin/bash
# Phase 2, step 8 (NEW strategy) -- RFdiffusion partial-diffusion soft-melt.
#
# Freezes the beta-strand framework AND the docked Key as fixed motifs and
# gently noises ONLY the contacting loops (diffuser.partial_T=10, ~20% of the
# 50-step trajectory). No ppi.hotspot_res: in partial diffusion the fixed Key
# motif already defines the interface (hotspots are a de-novo binder lever).
#
# Run src/prep_partial_diffusion.py FIRST -- it writes data/raw/scaffold_clean.pdb
# and the exact-length contig (data/raw/contig.txt) that this script consumes.
# Designs are written incrementally (loopmelt_0.pdb, loopmelt_1.pdb, ...), so an
# interrupted run still leaves usable backbones.
set -eu

OUTDIR=data/processed/rfdiffusion
mkdir -p "$OUTDIR"

INPUT=data/raw/scaffold_clean.pdb
CONTIG=$(cat data/raw/contig.txt)

echo "input : $INPUT"
echo "contig: $CONTIG"
echo "partial_T=10  num_designs=100  (no hotspots)"

# Double-quote the contig override: the contig contains a space (the '/0 ' chain
# break) and must reach Hydra as a single argument.
python RFdiffusion/scripts/run_inference.py \
  inference.output_prefix="$OUTDIR/loopmelt" \
  inference.input_pdb="$INPUT" \
  inference.model_directory_path=models/weights \
  "contigmap.contigs=$CONTIG" \
  diffuser.partial_T=10 \
  inference.num_designs=100
