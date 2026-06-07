#!/bin/bash
# Phase 2 (ligand-aware rebuild) -- RFdiffusionAA: diffuse a protein pocket
# AROUND the corrected fluoroaromatic anchor, so the BACKBONE SHAPE can
# complement the para-F. This is the step the LigandMPNN dry-run pointed to:
# F-vs-PHE discrimination has to come from pocket shape/physics, not from MPNN
# sequence preference (a single para-F barely moved LigandMPNN; see PROGRESS_LOG).
#
# LIGAND INPUT: data/raw/ligand_L4F.pdb  (resname PFF; 8 ring/CB/CA carbons +
#   the para-F; CONECT records included). Verified against the RFdiffusionAA
#   source: make_indep() -> filter_het() matches HETATM by resname == inference.
#   ligand, then parse_mol() reads it with OpenBabel -- atom types come from the
#   ELEMENT column and bonds from CONECT/geometry. There is NO CCD/template
#   lookup, so "PFF" is just a label and our side-chain-only anchor (no free
#   termini) is exactly what we want a pocket built around. Nothing to change.
#
# =============================================================================
# ONE-TIME SETUP  (heavy; ~13 GB + a system install. NOT run by this script.)
# =============================================================================
#   RFdiffusionAA is Apptainer-only upstream (no conda/pip env is provided).
#   1) Install Apptainer in WSL2 (needs sudo). Ubuntu .deb route:
#        https://apptainer.org/docs/admin/main/installation.html
#   2) Container (~11.8 GB) -> into rf_diffusion_all_atom/:
#        wget -P rf_diffusion_all_atom/ \
#          http://files.ipd.uw.edu/pub/RF-All-Atom/containers/rf_se3_diffusion.sif
#   3) Weights (~1.27 GB) -> into rf_diffusion_all_atom/ (aa.yaml default name):
#        wget -P rf_diffusion_all_atom/ \
#          http://files.ipd.uw.edu/pub/RF-All-Atom/weights/RFDiffusionAA_paper_weights.pt
#   4) Submodule lib/rf2aa is already initialized (cloned with --recurse-submodules).
#
# VRAM: this box is an 8 GB RTX 4060 (~5.7 GB free). NVIDIA's RFdiffusion floor
#   is ~12 GB, but that is for large designs; a SMALL de-novo pocket (NRES~80 +
#   9 ligand atoms) is tiny and is expected to fit 8 GB. Start small, scale up.
# =============================================================================
set -eu

REPO="$(cd "$(dirname "$0")/.." && pwd)"
AA="$REPO/rf_diffusion_all_atom"
SIF="$AA/rf_se3_diffusion.sif"
WEIGHTS="$AA/RFDiffusionAA_paper_weights.pt"      # aa.yaml ckpt_path default
LIGSRC="$REPO/data/raw/ligand_L4F.pdb"

# ---- design parameters (edit here) ------------------------------------------
NRES=80                       # pocket length: contig ['NRES-NRES'] (small=fits 8GB)
T=100                         # denoising steps (README uses 100-200)
NUM=1                         # number of designs (raise once the first succeeds)
LIGCODE=PFF                   # MUST equal the HETATM resname in the ligand PDB
GPU="--nv"                    # "--nv"=GPU; set GPU="" to run on CPU (much slower)
OUTPRE="output/l4f_pocket/des"
# -----------------------------------------------------------------------------

# stage the ligand INSIDE the container bind mount (cwd=$AA at run time)
mkdir -p "$AA/input" "$AA/$(dirname "$OUTPRE")"
cp "$LIGSRC" "$AA/input/anchor_L4F.pdb"

# preflight: fail clearly if the heavy one-time setup is not done yet
command -v apptainer >/dev/null 2>&1 || { echo "ERROR: apptainer not installed (see ONE-TIME SETUP)"; exit 1; }
[ -s "$SIF" ]     || { echo "ERROR: missing container $SIF (~11.8 GB; see ONE-TIME SETUP)"; exit 1; }
[ -s "$WEIGHTS" ] || { echo "ERROR: missing weights $WEIGHTS (~1.27 GB; see ONE-TIME SETUP)"; exit 1; }

echo "input  : input/anchor_L4F.pdb   ligand=$LIGCODE"
echo "contig : ['${NRES}-${NRES}']    T=$T   num_designs=$NUM   GPU='${GPU}'"
echo "outputs: $AA/${OUTPRE}_*.pdb"

# Run from $AA so run_inference.py and ./input ./output resolve inside the
# default ($PWD) bind mount. Matches the upstream README example exactly.
cd "$AA"
apptainer run $GPU "$SIF" -u run_inference.py \
  inference.deterministic=True \
  diffuser.T="$T" \
  inference.output_prefix="$OUTPRE" \
  inference.input_pdb=input/anchor_L4F.pdb \
  inference.ligand="$LIGCODE" \
  "contigmap.contigs=['${NRES}-${NRES}']" \
  inference.num_designs="$NUM" \
  inference.design_startnum=0
