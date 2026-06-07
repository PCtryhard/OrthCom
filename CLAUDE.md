# CLAUDE.md — OrthCom Pipeline

Guidance for Claude working in this repo. Read this first. **`PROGRESS_LOG.md` is the
authoritative source of truth** for what was actually tried and why; `README.md` is the
human-facing narrative. If this file ever disagrees with `PROGRESS_LOG.md`, the log wins —
and fix this file.

> This file was rewritten on 2026-06-07 for the **ligand-aware rebuild + 2-Nal pivot**. The
> earlier version described the **partial-diffusion** era as current — that era is now
> deprecated (see "Strategy"). If you find partial-diffusion-as-current guidance anywhere,
> it is stale.

## What OrthCom is

A **methods project** with one question: can a computationally designed receptor pocket
**discriminate a synthetic non-canonical amino acid (ncAA) from the canonical residues** —
bind it with a measurable **affinity gap** over its closest natural look-alike, not merely
"bind" it? The win condition is the *gap*. Everything is *in silico*.

Terminology (retained): **Lock** = the designed receptor pocket; **Key** = the short peptide
whose central recognition "tooth" is the ncAA. The earlier grand framing ("orthogonal
communication axis / signalling platform") has been retired as marketing.

## Current target (the 2-Nal pivot)

- **PRIMARY = L-2-naphthylalanine (2-Nal):** a fused **bicyclic** aromatic side chain, much
  larger / differently shaped than any canonical residue → real geometry for a pocket to
  grip, and no canonical look-alike to cross-react with. Commercially available.
- **NEGATIVE CONTROL / baseline = L-4-fluorophenylalanine (L4F):** kept *precisely* as the
  documented ncAA that proved **too subtle** to discriminate (a single para-F moved
  LigandMPNN sequence choice ~0.003 bits ≈ noise). **L4F is no longer the target.**
- **Canonical comparator = PHE** (the F-absent / unsubstituted aromatic). Trp is the tougher
  bicyclic look-alike worth testing against 2-Nal.

## Strategy: ligand-aware design is CURRENT; two earlier eras are DEPRECATED

When asked to advance "the pipeline," assume the **ligand-aware** strategy. Do **not**
silently re-run the partial-diffusion or de-novo scripts.

- **CURRENT — ligand-aware design.** Condition every learnable step on the *actual ncAA
  ligand*: **RFdiffusionAA** diffuses a backbone *around* the ligand (pocket shape) →
  **LigandMPNN** designs the pocket sequence conditioned on the ligand atoms → the
  **AF2-multimer co-fold instrument** scores the interface → **physics** scores the
  ncAA-vs-canonical **affinity gap** (the discrimination test).
- **DEPRECATED — era 2, scaffold partial diffusion.** Dock the helical Key into the rigid
  4QAE Anticalin barrel and soft-melt the loops (`diffuser.partial_T`). Abandoned: the
  helix-vs-barrel-loop interface was frustrated (AF2 deformed the helix, hard clashes, Key
  pLDDT ~30), and a malformed anchor was feeding it. Scripts: `prep_partial_diffusion.py`,
  `run_partial_diffusion.sh`.
- **DEPRECATED — era 1, de-novo hallucination.** Hallucinate 100–200 residues around the
  fixed 15-mer Key. Abandoned: "pool-noodle" helices, not pockets. Script: `run_rfdiffusion.sh`.

Why the pivots (L4F→2-Nal, and off partial diffusion): see "What didn't work" in `README.md`
and the decision log in `PROGRESS_LOG.md`.

## Canonical artifacts (and landmines)

`data/raw/` is git-ignored **except** `scaffold_complex.pdb`. Be precise about which
Key/ligand you use.

**Design ligands** (Stage 1; built by `src/prep_ligand.py` — side-chain-centric, **0 free
oxygens**, no H):
- `ligand_2NAL.pdb` — **PRIMARY** 2-naphthylalanine side chain (to be built in this pivot — Task B).
- `ligand_L4F.pdb` — 4-fluoro-Phe side chain, resname **PFF**, 9 heavy atoms incl. the
  para-F. Now the **negative control**, not the target.
- `ligand_PHE.pdb` — PHE control, 8 atoms (= L4F minus the F); differs by exactly one atom.

**Key structures:**
- `CANON.pdb` — poly-Ala 15-mer with **PHE8**; the authoritative canonical (Trojan-Horse) Key.
- `synthetic_key.pdb` — **L4F8**; physics steps only. Feeding to a DL tool crashes (`KeyError: 'L4F'`).
- `synthetic_key_rfd.pdb` — **TYR8**; STALE/WRONG (Tyr was rejected — hydroxyl polar penalty). Do not use.
- `scaffold_complex.pdb` *(tracked)* — chain A = 172-res **4QAE Anticalin** barrel (res 6–177),
  chain B = 15-mer Key (poly-Ala + PHE8). Confirmed 4QAE, **not** the 1mel VHH.
- `scaffold_clean.pdb` — cleaned `scaffold_complex` (waters/ANISOU/`CRYST1` stripped). The
  pocket/instrument scripts read the **frozen Key from its chain B**. NOTE its baked
  anchor is the *malformed* one (from `graft_residue.py`); the corrected anchor lives in the
  `ligand_*.pdb` files, not here.
- `data/raw/starters/1mel.pdb1` — a camelid VHH + lysozyme; an *alternative* starter, NOT the
  committed scaffold.

**Cloud package** *(tracked)*: `cloud/rfdiffusion_aa/` — `README.md` + a baked input PDB
(`input_L4F_pocket.pdb`; a 2-Nal sibling to be added in this pivot — Task C).

**Script landmines:**
- Deep-learning tools can't represent the ncAA by residue name (e.g. `KeyError: 'L4F'`).
  Peptide-based DL steps therefore present a plain **PHE** "Trojan-Horse" at position 8 (the
  AF2 co-fold reads the Key as PHE8). The ligand-aware steps instead see the **real ncAA** as
  a HETATM **ligand** (OpenBabel, no name lookup); the full physics returns only at the
  affinity-gap stage.
- `run_mpnn.sh` (ProteinMPNN) designs chain **A** (Lock), freezes chain **B** (Key). The
  current sequence stage uses **LigandMPNN** (ligand-conditioned), not this — `run_mpnn.sh` is
  an earlier-era tool.
- `run_alphafold.sh` is **monomer** ColabFold at `--num-recycle 0` (an old memory workaround;
  weak confidence — coarse filter only). The **interface** test is the multimer instrument
  below, a different script with different settings.

## Mechanics that are authoritative now

**SUPERSEDED:** the detailed RFdiffusion *partial-diffusion* contig mechanics (motif-ref vs
bare-length, `partial_T`, `provide_seq`, `ppi.hotspot_res`) that this file used to document.
Partial diffusion is no longer the strategy. If ever needed, they live in `RFdiffusion/README.md`
and the git history.

Current ligand-aware mechanics:
- **RFdiffusionAA needs a carrier protein, not a bare ligand.** Its `process_target()` builds
  the contig frame from protein CA atoms; with a **de-novo** contig the carrier is parsed then
  **discarded** and only the ligand is kept as the motif the fresh pocket is built around (the
  upstream 7v11/OQO pattern). `src/prep_rfdiffusion_aa_input.py` bakes carrier (scaffold chain
  A) + Kabsch-placed ligand into the cloud input.
- **RFdiffusionAA reads the ligand via OpenBabel** — atom identity from the **element column**,
  bonds from **CONECT/geometry**, **no CCD/template lookup**. So `inference.ligand=PFF` is just
  a label matching the HETATM resname; the para-F (and 2-Nal's fused rings) are read directly.
  Nothing to edit for atom-name conventions.
- **Side-chain-centric anchor ("Option A"):** the design ligand is the side chain + a CA cap,
  with **zero free-terminus oxygens** and no hydrogens. The ncAA is mid-peptide in the real
  Key (no free carboxylate); a buried charged carboxylate would swamp the subtle discrimination
  signal. `src/prep_ligand.py` builds from verified CCD ideal coords and geometry-self-checks.
- **RFdiffusionAA is CLOUD-ONLY (decision).** It wants ≥12–16 GB VRAM + Apptainer + ~13 GB of
  container/weights; the local box is 8 GB. `cloud/rfdiffusion_aa/` is the self-contained
  copy-paste package; `src/run_rfdiffusion_aa.sh` is the repo-integrated, preflight-gated
  runner. **No cloud run is launched without an explicit go.**
- **LigandMPNN finding (shapes the whole plan):** with ligand atom-context ON, a single para-F
  moved sequence choice ~0.003 bits (≈ noise) vs up to ~0.60 bits for whole-ligand context.
  → discrimination must come from pocket **shape** (RFdiffusionAA) and/or **physics**, NOT from
  MPNN AA preference. `--ligand_mpnn_use_atom_context 0` is a clean built-in ablation control.
- **AF2-multimer instrument (this WORKS — preserve it):** `prep_multimer.py` →
  `run_multimer.sh` → `parse_interface.py`. Co-fold `LockSeq:KeySeq` with
  `alphafold2_multimer_v3`, `--msa-mode single_sequence`, up to 20 recycles (0.5 Å early-stop).
  On the 8 GB box it is stable ONLY with JAX's `platform` allocator (the default BFC allocator
  fragments and segfaults at recycle 2). Rank by **interface-PAE** (both cross-chain blocks),
  tie-break **iPTM**. (The old `prep_multimer.py` `KEY_SEQUENCE="REPLACE_ME"` placeholder is
  gone — the Key sequence is now read programmatically from `scaffold_clean.pdb` chain B.)

## Environment / hardware (still true)

- GPU: **NVIDIA RTX 4060 Laptop, 8 GB** (~5.7 GB free). Fine for the AF2-multimer instrument
  and LigandMPNN (which runs CPU-only here anyway); **not** for RFdiffusionAA → cloud.
- conda + GPU live in **WSL2**; the agent's Bash is **git-bash/MINGW** (conda is not on its
  PATH — existing scripts call tools by absolute path). The `ligandmpnn` conda env is installed
  (CPU torch 2.2.1).
- Weights present: `models/weights/Base_ckpt.pt`, `Complex_base_ckpt.pt` (RFdiffusion —
  deprecated path). RFdiffusionAA `.sif`/weights are **not** downloaded locally (cloud-only).

## Repo layout (project files; ignore the vendored `RFdiffusion/`, `ProteinMPNN/`, `LigandMPNN/`, `rf_diffusion_all_atom/`, `localcolabfold/` trees)

- `src/` — full annotated map in `README.md`. Current ligand-aware: `prep_ligand.py`,
  `prep_rfdiffusion_aa_input.py`, `run_rfdiffusion_aa.sh`, `prep_ligandmpnn_input.py`,
  `analyze_ligandmpnn_dryrun.py`. **Instrument (preserve):** `prep_multimer.py`,
  `run_multimer.sh`, `parse_interface.py`. **Physics (preserve; partial):**
  `generate_conformer.py`, `run_psi4_dft.py`, `build_scaffold.py`, `graft_residue.py`,
  `Dummy/graft_dummmy_residue.py` (note: three `m`s). **Deprecated eras:**
  `prep_partial_diffusion.py`, `run_partial_diffusion.sh`, `run_rfdiffusion.sh`, `run_mpnn.sh`,
  `run_alphafold.sh`, `parse_scores.py`, `pocket_geometry.py` (volume calc is a stub).
- `cloud/rfdiffusion_aa/` — tracked cloud package (README + baked input PDB).
- `data/raw/` — inputs/keys (above). `data/processed/{rfdiffusion,mpnn,multimer,ligandmpnn}` —
  stage outputs (git-ignored).
- `outputs/<date>/<time>/` — Hydra run dirs (deprecated RFdiffusion path); `.hydra/overrides.yaml`
  is the cleanest record of what a run actually used.
- `PROGRESS_LOG.md` — the authoritative decision log. `README.md` — the human narrative.

## Conventions

- Don't commit `data/`, `outputs/`, `models/weights/`, or the vendored tool trees.
- **Append every decision to `PROGRESS_LOG.md`** (what / why / effect on the codebase; newest
  at the bottom). Keep the "What didn't work" institutional memory intact.
- When changing a run, prefer editing the `src/*.sh` script over ad-hoc CLI so
  `outputs/.../overrides.yaml` stays an honest record.
- **Deprecate by renaming, never by deleting.** Preserve the instrument and the physics path.
- **No heavy spend** (cloud runs, multi-GB downloads, GPU jobs) without showing the exact
  command and getting an explicit yes.
- Clean any scratch/helper scripts you create from the repo root when done.
