# OrthCom — The Orthogonal Communication Axis

**Mission:** build a perfectly isolated synthetic signalling channel for human biology — a "Lock and Key" interface that lets a designed peptide address a designed receptor pocket with high thermodynamic stability and zero cross-reactivity against the native proteome.

- **The Lock (hardware):** a receptor pocket engineered to engulf the Key. *Currently realised computationally as an engineered lipocalin (Anticalin) beta-barrel whose binding loops are re-sculpted around the Key.*
- **The Key (software):** a short synthetic peptide carrying a non-canonical amino acid (**L-4-fluorophenylalanine, L4F**) as its central recognition "tooth" — designed to be orthogonal to native human proteins.

This repository contains **Phase 1 (Key physical chemistry)** and **Phase 2 (Lock structural design)** — the *in silico* stack only.

---

## Strategy (read this before running anything)

We **pivoted** from the original approach. Both are described here so the history is clear, but only the current one is wired up.

- **Current — Scaffold partial diffusion (rigid start, low freedom).** Start from a **pre-existing rigid hydrophobic binding site** (the 4QAE Anticalin barrel), dock the 15-mer Key into it, then use RFdiffusion **partial diffusion** (`diffuser.partial_T`) to gently *soft-melt only the binding loops* around the Key while the beta-strand framework stays frozen. This constrains the model instead of letting it hallucinate freely.
- **Deprecated — De-novo hallucination (high freedom).** Keep the Key as a fixed motif and hallucinate 100–200 brand-new residues around it (`contigmap.contigs=[A1-15/0 150-200]`, `ppi.hotspot_res=[...]`). This produced long alpha-helix "pool noodles" instead of pockets and needed ~1000 designs to find a usable concavity. Preserved only in the design history below.

**Roadmap:** (1) partial-diffusion soft-melt of the loops around the 15-mer Key while preserving contact with the central PHE/L4F "Trojan" anchor → (2) ProteinMPNN re-sequences the scaffold (Key frozen) → (3) AlphaFold validates → (4) **compress the Key from 15 → 3–5 residues by inpainting** → (5) **OpenMM** swaps the canonical PHE anchor back to the real L4F and scores binding thermodynamics. We need a *really good fit on the non-canonical amino acid*; the helical-rigidity rationale for the Key is expected to break at 3–5 residues and will be revisited then.

---

## Pipeline

### Phase 1 — The Key (physical chemistry)

| Step | Tool | Script | Output | Status |
|---|---|---|---|---|
| 1. 3D conformer of L4F | RDKit | [src/generate_conformer.py](src/generate_conformer.py) | `data/raw/L_4_fluorophenylalanine.mol` | OK |
| 2. Quantum electrostatics | Psi4 (B3LYP/6-31G\*) | [src/run_psi4_dft.py](src/run_psi4_dft.py) | SCF energy | **Partial — ESP grid not yet exported** |
| 3. Force-field translation | AmberTools (`antechamber`/RESP) | *(none)* | `*.mol2`, `*.frcmod` | **Unimplemented — current files are external** |
| 4. Helix scaffold | PeptideBuilder | [src/build_scaffold.py](src/build_scaffold.py) | `data/raw/helix_scaffold.pdb` | OK |
| 5. Graft anchor at residue 8 | BioPython + RDKit | [src/graft_residue.py](src/graft_residue.py) (L4F), [src/Dummy/graft_dummy_residue.py](src/Dummy/graft_dummy_residue.py) (PHE) | `data/raw/synthetic_key.pdb` (L4F8), `data/raw/CANON.pdb` (PHE8) | OK — see graft caveat |

Steps 2–3 feed **only** the final OpenMM thermodynamics stage and are not on the critical path for Lock design. They must be completed before that stage (see *Known limitations*).

### Phase 2 — The Lock (scaffold partial diffusion)

| Step | Tool | Script | Notes |
|---|---|---|---|
| 6. Assemble complex | (manual) | — | 4QAE Anticalin (chain A) + Key (chain B) → `data/raw/scaffold_complex.pdb` |
| 7. Prep + contig | BioPython | [src/prep_partial_diffusion.py](src/prep_partial_diffusion.py) | strips waters/ANISOU; emits `scaffold_clean.pdb` + length-checked contig |
| 8. Loop soft-melt | RFdiffusion (partial_T) | [src/run_partial_diffusion.sh](src/run_partial_diffusion.sh) | strands+Key frozen, contacting loops diffused. **No `ppi.hotspot_res`** |
| 9. Sequence design | ProteinMPNN | [src/run_mpnn.sh](src/run_mpnn.sh) | designs scaffold (chain **A**), freezes Key (chain **B**) |
| 10. Validation | AlphaFold2 / ColabFold | [src/run_alphafold.sh](src/run_alphafold.sh) | single-sequence mode; see recycle caveat |
| 11. Complex check | ColabFold multimer | [src/prep_multimer.py](src/prep_multimer.py) | Lock + Key co-fold (pLDDT/iPTM) |

#### The "Trojan Horse" anchor
RFdiffusion/AlphaFold are trained on the 20 canonical residues and crash on L4F (`KeyError: 'L4F'`). For all deep-learning steps the Key presents a plain **Phenylalanine (PHE)** at position 8 (the fluorine-free stand-in). PHE — not Tyr — is deliberate: Tyr's hydroxyl would make ProteinMPNN build a polar H-bond network that later repels the hydrophobic fluorine. The real L4F is reintroduced only in the OpenMM stage.

---

## Repository layout

```
src/                       Pipeline scripts
  generate_conformer.py    Phase 1, step 1 (RDKit)
  run_psi4_dft.py          Phase 1, step 2 (Psi4)  — see status above
  build_scaffold.py        Phase 1, step 4 (PeptideBuilder)
  graft_residue.py         Phase 1, step 5 — L4F key  -> synthetic_key.pdb
  Dummy/graft_dummy_residue.py   Phase 1, step 5 — PHE key -> CANON.pdb
  prep_partial_diffusion.py      Phase 2, step 7 — clean + contig generator
  run_partial_diffusion.sh       Phase 2, step 8 — RFdiffusion partial_T
  run_mpnn.sh              Phase 2, step 9 — ProteinMPNN (designs chain A)
  run_alphafold.sh         Phase 2, step 10 — ColabFold
  prep_multimer.py         Phase 2, step 11 — Lock+Key multimer FASTAs
  parse_scores.py          Rank AlphaFold outputs by pLDDT
  pocket_geometry.py       Pocket I/O helpers (volume calc is a stub)
scripts/download_weights.py    Fetch RFdiffusion checkpoints
data/raw/                  Inputs (git-ignored except scaffold_complex.pdb)
data/processed/            Stage outputs (git-ignored)
models/weights/            RFdiffusion checkpoints (git-ignored)
tests/                     pytest
```

Vendored tools (`RFdiffusion/`, `ProteinMPNN/`, `localcolabfold/`) and `outputs/` are **git-ignored** — install them locally per Setup.

## Setup

```bash
# 1. Python deps for the helper scripts
pip install -r requirements.txt        # numpy scipy biopython pytest (+ rdkit, psi4 via conda for Phase 1)

# 2. Vendored design tools (each has its own conda env / install)
git clone https://github.com/RosettaCommons/RFdiffusion.git
git clone https://github.com/dauparas/ProteinMPNN.git
bash src/setup_colabfold.sh             # installs localcolabfold/ (WSL/Linux)

# 3. RFdiffusion model weights -> models/weights/
python scripts/download_weights.py      # Base_ckpt.pt, Complex_base_ckpt.pt
```

## Running the first partial diffusion

`data/raw/scaffold_complex.pdb` (4QAE Anticalin chain A, residues 6–177 + 15-mer Key chain B) is committed and ready.

```bash
# 1. Clean the complex and auto-generate a length-checked contig
python src/prep_partial_diffusion.py
#   -> data/raw/scaffold_clean.pdb   (ATOM only; waters/ANISOU stripped)
#   -> data/raw/contig.txt           (strands+Key frozen, contacting loops melted)

# 2. Soft-melt the binding loops around the Key
bash src/run_partial_diffusion.sh       # RFdiffusion, diffuser.partial_T=10, no hotspots
```
`prep_partial_diffusion.py --whole` instead emits a whole-scaffold gentle-melt contig (`[172-172/0 B1-15]`) for a plumbing sanity check.

The default loop-melt flags scaffold residues within `--contact-threshold` (6 Å) of the Key, then tidies the raw contact mask into coherent loops: `--close-gap` (merge fragments split by short frozen gaps, default 2), `--min-run` (drop lone strand-side-chain contacts, default 3), and `--pad` (grow each loop for strand-junction freedom, default 1). On the committed complex this yields 6 loops (44/172 residues melted):
```
[A6-34/9-9/A44-47/6-6/A54-66/5-5/A72-77/5-5/A83-92/6-6/A99-122/13-13/A136-177/0 B1-15]
```
Use `--pad 0 --min-run 1 --close-gap 0` to recover the raw per-residue contact mask.

## Known limitations / TODO

- **ESP + AmberTools (Phase 1, steps 2–3) unimplemented.** Needed before the OpenMM swap-back, not before Lock design.
- **Graft matcher is heuristic.** Validate PHE8/L4F8 backbone + χ dihedrals before reintroducing fluorine.
- **AlphaFold runs at `--num-recycle 0`** (WSL memory workaround): treat pLDDT as a coarse filter, not proof.
- **Partial diffusion cannot change loop length or clear a hard clash.** If native 4QAE loops are too short to engulf a 15-mer, switch from `partial_T` to loop *rebuild* (motif scaffolding, full T, variable loop length).
- **Key compression (15 → 3–5mer).** A 3–5mer cannot hold a helix; expect to need a constrained/stapled or macrocyclic micro-key.
- **Hardware ceiling:** validated on an 8 GB RTX 4060 Laptop GPU.

## Scaffold provenance
The Lock scaffold is **PDB 4QAE** — "Crystal structure of an engineered lipocalin (Anticalin) in complex with human hepcidin." We reuse the rigid 8-stranded barrel and re-sculpt its four hepcidin-binding loops around our Key.

## Design history

1. **L4F dictionary crash** → convert the Key's anchor to PHE for all deep-learning steps.
2. **"Pool noodle" helices** (de-novo) → the single-hotspot de-novo route collapsed to long helices; multi-point "claw" hotspots helped but never reliably produced pockets — motivating the scaffold pivot.
3. **Contig fusion** → de-novo runs fused Lock onto Key without the explicit `/0` chain break.
4. **MSA timeout** → AlphaFold runs in `--msa-mode single_sequence` (synthetic proteins have no evolutionary history).
5. **Pivot to scaffold partial diffusion** (current) — start rigid, melt only the loops.

> Deprecated de-novo command, kept for reference:
> `contigmap.contigs=[A1-15/0 150-200] ppi.hotspot_res=[A1,A4,A8,A12,A15]` on `CANON.pdb`.
