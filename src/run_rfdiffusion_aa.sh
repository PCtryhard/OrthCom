#!/bin/bash
# Phase 2 (ligand-aware rebuild) -- RFdiffusionAA: diffuse a protein pocket
# AROUND the corrected fluoroaromatic anchor, so the BACKBONE SHAPE can
# complement the para-F. This is the step the LigandMPNN dry-run pointed to:
# F-vs-PHE discrimination has to come from pocket shape/physics, not from MPNN
# sequence preference (a single para-F barely moved LigandMPNN; see PROGRESS_LOG).
#
# CLOUD-ONLY (your decision 2026-06-07): the 8 GB local box is NOT used for
# RFdiffusionAA. Run this on a >=16 GB GPU host with Apptainer. The command is
# host-agnostic; the portable copy-paste package is cloud/rfdiffusion_aa/README.md.
#
# INPUT: cloud/rfdiffusion_aa/input_L4F_pocket.pdb (tracked; baked-in F-ligand).
#   = scaffold chain A (CARRIER protein) + the corrected L4F anchor as HETATM
#   resname PFF (8 ring/CB/CA carbons + the para-F, with CONECT). Verified vs the
#   RFdiffusionAA source: make_indep->filter_het matches HETATM by resname ==
#   inference.ligand, parse_mol reads it with OpenBabel (atom types from element,
#   bonds from CONECT/geometry, NO CCD lookup). A carrier protein is required
#   because process_target() builds the contig frame from protein CA atoms; with a
#   DE-NOVO contig the carrier is parsed then DISCARDED (only the ligand is kept as
#   motif) -- exactly the upstream 7v11/OQO pattern. Regenerate with
#   src/prep_rfdiffusion_aa_input.py.
#
# =============================================================================
# ONE-TIME SETUP on the cloud host (heavy; ~13 GB + a system install)
# =============================================================================
#   RFdiffusionAA is Apptainer-only upstream (no conda/pip env is provided).
#   1) git clone --recurse-submodules \
#        https://github.com/baker-laboratory/rf_diffusion_all_atom.git && cd $_
#   2) Container (~11.8 GB):
#        wget http://files.ipd.uw.edu/pub/RF-All-Atom/containers/rf_se3_diffusion.sif
#   3) Weights (~1.27 GB, aa.yaml default name, in the repo dir):
#        wget http://files.ipd.uw.edu/pub/RF-All-Atom/weights/RFDiffusionAA_paper_weights.pt
#   4) Install Apptainer: https://apptainer.org/docs/admin/main/installation.html
# =============================================================================
set -eu

REPO="$(cd "$(dirname "$0")/.." && pwd)"
AA="$REPO/rf_diffusion_all_atom"
SIF="$AA/rf_se3_diffusion.sif"
WEIGHTS="$AA/RFDiffusionAA_paper_weights.pt"          # aa.yaml ckpt_path default
INPUT_SRC="$REPO/cloud/rfdiffusion_aa/input_L4F_pocket.pdb"

# ---- design parameters (edit here) ------------------------------------------
NRES=120                      # de-novo pocket length: contig ['NRES-NRES'] (80-150 sane)
T=100                         # denoising steps (README uses 100-200)
NUM=1                         # number of designs (raise once the first succeeds)
LIGCODE=PFF                   # MUST equal the HETATM resname in the input PDB
GPU="--nv"                    # "--nv"=GPU; set GPU="" for CPU (very slow)
OUTPRE="output/l4f_pocket/des"
# -----------------------------------------------------------------------------

# stage the input INSIDE the container bind mount (cwd=$AA at run time)
mkdir -p "$AA/input" "$AA/$(dirname "$OUTPRE")"
cp "$INPUT_SRC" "$AA/input/input_L4F_pocket.pdb"

# preflight: fail clearly if the heavy one-time setup is not done yet
command -v apptainer >/dev/null 2>&1 || { echo "ERROR: apptainer not installed (see ONE-TIME SETUP)"; exit 1; }
[ -s "$SIF" ]     || { echo "ERROR: missing container $SIF (~11.8 GB; see ONE-TIME SETUP)"; exit 1; }
[ -s "$WEIGHTS" ] || { echo "ERROR: missing weights $WEIGHTS (~1.27 GB; see ONE-TIME SETUP)"; exit 1; }

echo "input  : input/input_L4F_pocket.pdb   ligand=$LIGCODE (carrier protein discarded by de-novo contig)"
echo "contig : ['${NRES}-${NRES}']    T=$T   num_designs=$NUM   GPU='${GPU}'"
echo "outputs: $AA/${OUTPRE}_*.pdb"

# Run from $AA so run_inference.py and ./input ./output resolve inside the
# default ($PWD) bind mount. Matches the upstream README example exactly.
cd "$AA"
apptainer run $GPU "$SIF" -u run_inference.py \
  inference.deterministic=True \
  diffuser.T="$T" \
  inference.output_prefix="$OUTPRE" \
  inference.input_pdb=input/input_L4F_pocket.pdb \
  inference.ligand="$LIGCODE" \
  "contigmap.contigs=['${NRES}-${NRES}']" \
  inference.num_designs="$NUM" \
  inference.design_startnum=0
