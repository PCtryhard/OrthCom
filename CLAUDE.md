# CLAUDE.md — OrthCom Pipeline

Guidance for Claude when working in this repo. Read this before touching the protein-design pipeline.

## What OrthCom is

A computational pipeline to engineer a synthetic **"Lock and Key"** protein interface:
- **Key** = a 15-residue peptide (poly-Ala backbone) carrying a hydrophobic aromatic anchor at the central residue (position 8). The eventual physical key uses the non-canonical amino acid **L-4-fluorophenylalanine (L4F)**; for deep-learning steps the fluorine is removed so the anchor presents as plain **PHE** ("Trojan Horse").
- **Lock** = a receptor pocket that engulfs the Key with high thermodynamic stability and zero native cross-reactivity.

The toolchain is RDKit/Psi4/AmberTools (Key physics) → RFdiffusion (backbone) → ProteinMPNN (sequence) → AlphaFold2/ColabFold (validation). See `README.md` for the original step-by-step narrative.

## IMPORTANT: two competing strategies live in this repo

There has been a **strategy pivot**. `README.md` and the partial-diffusion scripts now reflect the NEW strategy; the OLD de-novo script (`run_rfdiffusion.sh`) is retained only as historical record:

1. **OLD — de-novo hallucination** (what every committed script + every recorded run in `outputs/` actually does): keep the 15-mer Key as a fixed motif and hallucinate 100–200 brand-new residues around it (`contigmap.contigs=[A1-15/0 150-200]`, `ppi.hotspot_res=[...]`). This produced long alpha-helix "pool noodles" instead of pockets and needs ~1000 designs to fish for a usable concavity.

2. **NEW — scaffold partial diffusion** (the current intent; scripted as `src/prep_partial_diffusion.py` -> `src/run_partial_diffusion.sh`): start from a **rigid pre-existing beta-barrel scaffold** with the Key already docked into its hydrophobic site, then use **partial diffusion (`diffuser.partial_T`)** to gently "soft-melt" only the variable binding loops around the Key while the beta-strand framework stays frozen. Endgame: compress the Key from 15→3-5 residues via inpainting once loops are stabilized.

When asked to advance "the pipeline," assume the **NEW** strategy unless told otherwise. Do not silently re-run the old de-novo scripts.

## Canonical artifacts (and landmines)

Key variants in `data/raw/` — **be precise about which one you use**:
- `CANON.pdb` — poly-Ala 15-mer with **PHE8**. This is the authoritative canonical Key (Trojan-Horse target). All recent de-novo runs used this.
- `synthetic_key.pdb` — has **L4F8**. Feeding this to RFdiffusion crashes it (`KeyError: 'L4F'`). Physics steps only, never deep learning.
- `synthetic_key_rfd.pdb` — has **TYR8**. STALE/WRONG: Tyr was explicitly rejected in the design (hydroxyl "polar penalty" repels the fluorine on swap-back). Do not use; candidate for deletion.
- `data/raw/scaffold_complex.pdb` (committed) — the assembled NEW-strategy complex: **chain A** = beta-barrel scaffold (protein residues **6–177**, 172 residues, contiguous) + **chain B** = the Key (poly-Ala + PHE8). The raw file carried crystallographic waters (HOH), ANISOU records, and a P31 `CRYST1`; `src/prep_partial_diffusion.py` strips these into `data/raw/scaffold_clean.pdb` (the actual RFdiffusion input).
- `data/raw/starters/1mel.pdb1` — a camelid VHH nanobody (chain A) bound to lysozyme (chain L); a *separate/alternative* starter, ~132-residue VHH. RESOLVED: the committed `scaffold_complex.pdb` is the **172-residue Anticalin (4QAE)**, not this ~132-res VHH — confirmed by chain-A residue count.

Script/path landmines:
- `src/run_rfdiffusion.sh` is the **DEPRECATED** de-novo path (kept for history). Its input was corrected from the nonexistent `canonical_key.pdb` to `CANON.pdb`; do not run it for the NEW strategy.
- `src/run_mpnn.sh` now designs chain **A** (`--pdb_path_chains "A"`) and extracts `split("/")[0]` — correct for the NEW topology (Scaffold/Lock=**A**, Key=**B**: design the Lock, freeze the Key). ProteinMPNN's `--pdb_path_chains` lists chains to **design**; all others are frozen context. (The OLD topology had Lock=B, hence the previous `"B"` / `[-1]`.)

## RFdiffusion partial-diffusion mechanics (authoritative, from `RFdiffusion/README.md`)

- Partial diffusion **requires the contig to be exactly the same total length as the input PDB** — you cannot use length ranges that change residue count.
- In the contig, a **motif reference** (e.g. `A6-29`, `B1-15`) is **kept fixed**; a **bare length** (e.g. `8-8`) is **diffused/noised**. So to freeze strands+Key and melt loops: write strands and the Key as motif refs and the loops as bare lengths whose sizes exactly equal the native loop lengths.
- `diffuser.partial_T` sets how far to noise (T=50 default). Low values (~5–15) = gentle "soft-melt."
- `contigmap.provide_seq=[start-end]` (zero-indexed) pins the *sequence* of a *diffused* chain; needed only if you also diffuse the Key but want its sequence fixed. Requires a different checkpoint (auto-selected) — verify it exists in `models/weights/`.
- **`ppi.hotspot_res` is a de-novo binder-design lever (Complex model), not a partial-diffusion control.** For partial diffusion on a pre-docked complex, the fixed Key motif defines the interface; drop hotspots.

Two limits to respect for the NEW strategy:
- partial_T **cannot change loop length** and a low partial_T **will not relieve a hard steric clash**. The Key must be docked clash-free first. If native loops are too short to engulf a 15-mer helix, you need loop *rebuild* (motif scaffolding / inpainting, full T, variable loop length), not partial diffusion.
- A 3–5mer cannot hold alpha-helical geometry (needs ~5+ residues/turn). The "rigid helix Key" rationale breaks at the compression endgame — expect to need a constrained/stapled or non-helical micro-key.

## Environment / hardware (observed working)

- GPU: **NVIDIA RTX 4060 Laptop, 8 GB VRAM** (from run logs). Fine for ~190-residue complexes; the 8 GB ceiling is the real bottleneck, not missing deps.
- RFdiffusion, ProteinMPNN, and ColabFold all have execution evidence — the conda envs (incl. SE3Transformer, Hydra) are installed and functional. Weights present: `models/weights/Base_ckpt.pt`, `Complex_base_ckpt.pt`.
- ColabFold (`src/run_alphafold.sh`) runs under WSL with `--num-recycle 0` to dodge a memory crash, `--msa-mode single_sequence` (correct for synthetic proteins). **0 recycles materially weakens AF2 confidence** — treat pLDDT from this config as a coarse filter, not proof.
- PyMOL (README Step 6, fluorine deletion) appears **superseded** by the RDKit graft path (`src/Dummy/graft_dummmy_residue.py`, which writes `CANON.pdb` with PHE8). PyMOL is not invoked by any committed script.

## Repo layout (project files; ignore the huge `localcolabfold/` conda tree)

- `src/` — pipeline scripts. Physics: `generate_conformer.py`, `run_psi4_dft.py`, `build_scaffold.py`, `graft_residue.py` (L4F), `Dummy/graft_dummmy_residue.py` (PHE→`CANON.pdb`). DL: `prep_partial_diffusion.py` (clean + length-checked contig), `run_partial_diffusion.sh` (NEW partial-diffusion soft-melt), `run_rfdiffusion.sh` (OLD de-novo, deprecated), `run_mpnn.sh`, `run_alphafold.sh`, `parse_scores.py`, `prep_multimer.py` (has a `KEY_SEQUENCE = "REPLACE_ME"` placeholder). `pocket_geometry.py` has stub `calculate_volume` (raises NotImplementedError).
- `data/raw/` — inputs/keys (above). `data/processed/{rfdiffusion,mpnn,alphafold}` — stage outputs.
- `outputs/<date>/<time>/` — Hydra run dirs; `.hydra/overrides.yaml` is the cleanest record of what was actually run. `run_inference.log` can be huge (cautious-mode skip spam).
- `RFdiffusion/`, `ProteinMPNN/`, `localcolabfold/` — vendored tools.

## Conventions

- Don't commit `data/`, `outputs/`, `models/weights/`, or the vendored tool trees.
- When changing a run, prefer editing the `src/*.sh` script over ad-hoc CLI so `outputs/.../overrides.yaml` stays an honest record.
- Clean any scratch/helper scripts you create from the repo root when done.
