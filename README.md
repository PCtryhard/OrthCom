# OrthCom — discriminating a non-canonical amino acid with a designed pocket

**This is a methods project, not a product.** The single question it tries to answer:

> Can a computationally designed receptor pocket **discriminate** a synthetic
> non-canonical amino acid (ncAA) from the canonical residues — i.e. bind the ncAA
> with a measurable **affinity gap** over its closest natural look-alike, rather than
> merely "binding" it?

The win condition is the *gap*, not raw affinity. A pocket that grips the ncAA tightly
but grips a canonical residue just as tightly has proven nothing. Everything below is
*in silico* (Phase 1 = ncAA physical chemistry, Phase 2 = pocket design + validation).

Terminology is the old "Lock and Key": the **Lock** is the designed receptor pocket;
the **Key** is the short peptide whose central recognition "tooth" is the ncAA. No
grander claims are made — the prior "orthogonal signalling axis / platform" framing has
been retired as marketing that outran the evidence.

---

## Target: 2-naphthylalanine (primary) · 4-fluorophenylalanine (control)

| Role | Residue | Why |
|---|---|---|
| **Primary target** | **L-2-naphthylalanine (2-Nal)** | A fused **bicyclic** aromatic side chain — substantially larger and differently shaped than any canonical residue, so a complementary pocket has real geometry to grip and there is no canonical look-alike to cross-react with. Commercially available (keeps eventual wet-lab work feasible). |
| **Negative control / baseline** | **L-4-fluorophenylalanine (L4F)** | Kept *precisely* as the documented example of an ncAA whose distinguishing feature (a single para-fluorine) is **too subtle** for the design tools to act on (see "What didn't work"). It is the project's calibration that "barely different from PHE" reads as noise. |
| **Canonical comparator** | **Phenylalanine (PHE)** | The F-absent / unsubstituted aromatic. Differs from `ligand_L4F` by exactly one atom and from 2-Nal by one fused ring — the natural canonical baseline to measure the gap against. |

**Why the pivot off L4F.** The earlier thesis was "discriminate L4F from PHE." A
LigandMPNN dry-run showed a single para-F moves sequence choice by ~**0.003 bits** —
essentially noise, ~50–100× smaller than whole-ligand context. One fluorine is too small
a feature to discriminate on. So the primary target moved to a side chain with genuine
shape/size (2-Nal), and L4F was demoted to the negative control. Full numbers in
[PROGRESS_LOG.md](PROGRESS_LOG.md) and in "What didn't work" below.

> **ACTION ITEM (on the books — also logged in [PROGRESS_LOG.md](PROGRESS_LOG.md)):** the
> most stringent *canonical* look-alike for a bicyclic ncAA is **tryptophan** — the only
> fused-bicyclic aromatic among the 20 — so **Trp, not PHE, is the real cross-reactivity
> threat**. Once a pocket exists, run the **2-Nal-vs-Trp** affinity-gap test explicitly: it
> is the single test most likely to falsify the 2-Nal story, so it is recorded now.

---

## The pipeline (current — ligand-aware)

The design is conditioned on the **actual ncAA ligand** at every learnable step, so the
models can "see" the feature we want discriminated.

| # | Stage | Tool | Script(s) | Status |
|---|---|---|---|---|
| 1 | **ncAA ligand anchor** | RDKit / CCD coords | [src/prep_ligand.py](src/prep_ligand.py) | L4F + PHE built; **2-Nal in progress** |
| 2 | **Pocket backbone** (diffuse a protein *around* the ligand) | RFdiffusionAA | [src/prep_rfdiffusion_aa_input.py](src/prep_rfdiffusion_aa_input.py), [src/run_rfdiffusion_aa.sh](src/run_rfdiffusion_aa.sh), [cloud/rfdiffusion_aa/](cloud/rfdiffusion_aa/) | **Staged, cloud-only — not yet run** |
| 3 | **Pocket sequence** (conditioned on ligand atoms) | LigandMPNN | [src/prep_ligandmpnn_input.py](src/prep_ligandmpnn_input.py), [src/analyze_ligandmpnn_dryrun.py](src/analyze_ligandmpnn_dryrun.py) | Env installed; **dry-run done** |
| 4 | **Interface validation** (the working instrument) | ColabFold AF2-multimer | [src/prep_multimer.py](src/prep_multimer.py), [src/run_multimer.sh](src/run_multimer.sh), [src/parse_interface.py](src/parse_interface.py) | **Working / validated** |
| 5 | **Affinity gap** (the discrimination test) | Psi4 / RESP / AmberTools → MD | [src/generate_conformer.py](src/generate_conformer.py), [src/run_psi4_dft.py](src/run_psi4_dft.py) | **Partial / deferred** |

### Stage 1 — the ncAA anchor (and why it's side-chain-only)
`src/prep_ligand.py` builds the ligand from authoritative **RCSB Chemical Component
Dictionary (CCD)** ideal coordinates and geometry-self-checks it. We use the
**side-chain-centric representation ("Option A")**: the side chain + the CA cap only,
with **zero free-terminus oxygens** and no hydrogens. Rationale: in the real Key the
ncAA is a *mid-peptide* residue (no free carboxylate/amino terminus), and a buried
charged carboxylate would dominate the energetics over the subtle shape/electronics we
actually want the pocket to read. Outputs (under `data/raw/`, git-ignored):

- `ligand_L4F.pdb` — 4-fluoro-Phe side chain, resname `PFF`, **9 heavy atoms incl. the para-F**.
- `ligand_PHE.pdb` — the same minus the F (**8 atoms**) — the canonical control; the two files differ by exactly one atom.
- `ligand_2NAL.pdb` — 2-naphthylalanine side chain *(to be built and verified in this pivot — Task B)*.

### Stage 2 — RFdiffusionAA pocket (cloud-only)
`prep_rfdiffusion_aa_input.py` rigid-places (Kabsch) the anchor into a **carrier protein**
and writes the cloud input. The carrier is required only so RFdiffusionAA's
`process_target()` can build a contig frame from protein CA atoms; with a de-novo contig
the carrier is **parsed then discarded** and only the ligand is kept as the motif the
fresh pocket is built around (the upstream 7v11/OQO pattern). RFdiffusionAA reads the
ligand via OpenBabel (atom identity from the element column, bonds from CONECT — **no CCD
lookup**), so the ncAA's distinguishing atoms are seen directly.

This step is **cloud-only by decision**: all-atom RFdiffusionAA wants ≥12–16 GB VRAM and
is Apptainer-only (~13 GB of container + weights). The local box is an 8 GB RTX 4060
Laptop GPU — fine for the LigandMPNN and AF2 stages, not for this one.
[cloud/rfdiffusion_aa/](cloud/rfdiffusion_aa/) is a self-contained copy-paste package
(setup + exact run command + verification stamp). **No cloud run has been launched.**

### Stage 4 — the validation instrument (this part works)
The interface is scored by **co-folding** the designed Lock with the frozen Key as a
two-chain complex and reading cross-chain metrics — not monomer fold quality.
`prep_multimer.py` pairs each designed Lock (chain A) sequence with the Key sequence read
straight from the scaffold (chain B) into a `LockSeq:KeySeq` FASTA; `run_multimer.sh`
folds it with `alphafold2_multimer_v3` (`--msa-mode single_sequence`, up to 20 recycles
with 0.5 Å early-stop); `parse_interface.py` ranks designs by **interface-PAE** (mean of
the two cross-chain PAE blocks) with **iPTM** as the tie-break. The 8 GB box only runs
this stably after switching JAX to the `platform` allocator — see the header of
`run_multimer.sh` for that hard-won detail.

### Stage 5 — the affinity gap (the actual discrimination test)
The headline result is a **thermodynamic gap**: score the designed pocket against the
ncAA vs. against the canonical comparators — **PHE and, as the decisive cross-reactivity
test, Trp** (the action item above) — and show the ncAA wins
by a margin. For 2-Nal this is the primary readout; the **L4F-vs-PHE** gap is the
documented baseline that the *sequence* lever could not produce. This uses the Phase 1
physics path (`generate_conformer.py`, `run_psi4_dft.py`, then RESP/AmberTools charges →
MD), which is **only partially built** (DFT SCF runs; ESP-grid export and the AmberTools
RESP/`frcmod` step are not implemented) and is **deferred** this pivot.

---

## What didn't work (and why)

This section is institutional memory — keep it. Newest first.

1. **A single fluorine is not a discriminable feature (for sequence design).**
   LigandMPNN, run on the docked pocket with the ligand atom-context ON, changed its
   sequence choice between L4F (para-F present) and PHE (absent) by a Jensen–Shannon
   divergence of **mean ~0.0034 bits** at F-proximal positions (max 0.0099) — versus up
   to **0.60 bits** for whole-ligand context on/off. At the ring-adjacent position the
   model picked Phe deterministically with *or* without the F. **Conclusion:** MPNN
   *sequence* design cannot discriminate one para-F; discrimination, if achievable, must
   come from pocket **shape** (RFdiffusionAA) and/or explicit **physics** scoring. This
   is why L4F is now the negative control and the primary target is the much larger,
   differently shaped 2-Nal.

2. **A helical Key docked into a barrel gave a frustrated interface ("era 2").**
   Docking the 15-mer α-helical Key into the rigid 4QAE Anticalin β-barrel and
   soft-melting the loops (RFdiffusion partial diffusion) never produced a clean grip.
   The input dock was shallow/one-sided (44/82 Key atoms buried; centre-of-mass ~13 Å
   off the pocket centre). AF2 co-folds buried the Key more but **deformed the helix**
   (helical windows dropped 12→2–8) and introduced hard clashes (down to 0.5 Å), with
   Key pLDDT stuck ~30 even while the Lock folded fine (complex pLDDT 68–76). The
   interface secondary structure (helix vs. a barrel-loop face) was fundamentally
   mismatched, and partial diffusion *cannot* change loop length or relieve a hard
   clash.

3. **The original anchor was chemically malformed.** The grafted anchor baked into
   `scaffold_clean.pdb` (from `graft_residue.py`) auto-named atoms by element count
   (`O1/O2/C3…C9`), kept a **spurious second carboxylate oxygen**, and kept explicit
   hydrogens — fine-ish for a free amino acid, wrong for a mid-peptide residue and
   unusable as a clean ligand. Fixed by `prep_ligand.py`, which rebuilds from verified
   CCD coordinates as a side-chain-only ligand with zero free-terminus oxygens.

4. **De-novo hallucination gave "pool-noodle" helices ("era 1").** Keeping the Key as a
   fixed motif and hallucinating 100–200 brand-new residues around it
   (`contigmap.contigs=[A1-15/0 150-200]`, with `ppi.hotspot_res`) collapsed to long
   α-helices rather than concave pockets; multi-point "claw" hotspots helped but never
   reliably produced a pocket, and it needed ~1000 designs to fish for one. This
   motivated the scaffold pivot (era 2), which in turn gave way to the ligand-aware
   approach (current).

**Smaller documented gotchas (kept so they aren't rediscovered):**
- RFdiffusion/AF2 crash on `L4F` (`KeyError: 'L4F'`) — deep-learning steps present a
  plain **PHE** "Trojan-Horse" stand-in; the fluorine returns only at the physics stage.
- **Tyr was rejected** as the stand-in: its hydroxyl makes MPNN build a polar H-bond
  network that later repels the hydrophobic fluorine on swap-back. PHE, not TYR.
- AF2 on the 8 GB box: the BFC allocator fragmented and segfaulted at recycle 2; the
  `platform` allocator fixed it. `TF_FORCE_UNIFIED_MEMORY` did nothing (ColabFold infers
  with JAX/XLA, not the TF GPU runtime).
- Synthetic proteins have no MSA → `--msa-mode single_sequence` is correct, but treat
  confidence accordingly.

---

## Repository map

Verified against the tree (2026-06-07). Vendored tool trees and `outputs/` are
git-ignored — install them per **Setup**.

```
src/
  # -- Phase 1: Key physical chemistry (feeds the Stage-5 affinity scoring) --
  generate_conformer.py          RDKit 3D conformer of the ncAA
  run_psi4_dft.py                Psi4 B3LYP/6-31G* electrostatics (SCF only; ESP grid NOT exported)
  build_scaffold.py              PeptideBuilder poly-Ala helix scaffold for the Key
  graft_residue.py               graft the L4F anchor at position 8 -> synthetic_key.pdb
  Dummy/graft_dummmy_residue.py  graft the PHE stand-in -> CANON.pdb   (note: 3 m's in the filename)
  Dummy/generate_dummy_conformer.py

  # -- Phase 2 (CURRENT): ligand-aware design --
  prep_ligand.py                 build the side-chain-centric ncAA ligand(s)  -> data/raw/ligand_*.pdb
  prep_rfdiffusion_aa_input.py   build the RFdiffusionAA cloud input (carrier protein + ligand)
  run_rfdiffusion_aa.sh          cloud-only RFdiffusionAA runner (Apptainer; preflight-gated)
  prep_ligandmpnn_input.py       Kabsch-place the ligand in the pocket + find pocket residues
  analyze_ligandmpnn_dryrun.py   4-condition (F-present/absent x context on/off) discrimination analysis

  # -- The validation instrument (AF2-multimer co-fold) --
  prep_multimer.py               pair designed Lock (chain A) with frozen Key (chain B) -> LockSeq:KeySeq FASTAs
  run_multimer.sh                AF2-multimer co-fold (alphafold2_multimer_v3, platform allocator)
  parse_interface.py             rank complexes by interface-PAE, then iPTM

  # -- Earlier strategy eras (superseded; kept for history, see "What didn't work") --
  prep_partial_diffusion.py      era 2: scaffold clean + length-checked partial-diffusion contig
  run_partial_diffusion.sh       era 2: RFdiffusion partial_T loop soft-melt
  run_rfdiffusion.sh             era 1: de-novo hallucination
  run_mpnn.sh                    eras 1-2: ProteinMPNN sequence design (current path uses LigandMPNN instead)
  run_alphafold.sh               eras 1-2: ColabFold *monomer* validation (interface test is run_multimer.sh)
  parse_scores.py                eras 1-2: rank monomer AF2 outputs by pLDDT
  pocket_geometry.py             pocket I/O helpers (calculate_volume is a NotImplementedError stub)
  setup_colabfold.sh             install localcolabfold/ (WSL/Linux)

cloud/rfdiffusion_aa/            self-contained cloud package: README.md + baked input PDB (TRACKED)
scripts/download_weights.py      fetch RFdiffusion checkpoints -> models/weights/
data/raw/                        inputs & keys (git-ignored except scaffold_complex.pdb)
data/processed/{rfdiffusion,mpnn,multimer,ligandmpnn,...}   stage outputs (git-ignored)
models/weights/                  RFdiffusion checkpoints (git-ignored)
tests/test_pocket.py             pytest
PROGRESS_LOG.md                  the authoritative decision log (read this for the "why")
```

**Key artifacts in `data/raw/` — be precise about which one you use:**
- `scaffold_complex.pdb` *(tracked)* — the assembled complex: **chain A** = 172-residue
  4QAE Anticalin β-barrel scaffold (residues 6–177), **chain B** = the 15-mer poly-Ala
  Key with the central anchor. The repo's one runnable-out-of-the-box input.
- `scaffold_clean.pdb` — `prep_partial_diffusion.py`'s cleaned version (waters/ANISOU/`CRYST1` stripped). The pocket/instrument scripts read the frozen Key from its chain B.
- `CANON.pdb` — poly-Ala 15-mer with **PHE8**: the authoritative canonical (Trojan-Horse) Key.
- `synthetic_key.pdb` — has **L4F8**; physics steps only — feeding it to a DL tool crashes (`KeyError: 'L4F'`).
- `synthetic_key_rfd.pdb` — has **TYR8**; **stale/wrong** (Tyr was rejected), candidate for deletion.
- `ligand_L4F.pdb` / `ligand_PHE.pdb` — the Stage-1 design ligands built by `prep_ligand.py`; **`ligand_2NAL.pdb` is pending (Task B)**.

## Setup

```bash
# 1. Python deps for the helper scripts
pip install -r requirements.txt          # numpy scipy biopython pytest (+ rdkit/psi4 via conda for Phase 1)

# 2. Vendored design tools (each clones separately; git-ignored)
git clone https://github.com/RosettaCommons/RFdiffusion.git          # eras 1-2 (deprecated path)
git clone https://github.com/dauparas/ProteinMPNN.git                # eras 1-2
git clone https://github.com/dauparas/LigandMPNN.git                 # current sequence stage (+ get_model_params.sh)
bash src/setup_colabfold.sh                                          # installs localcolabfold/ (WSL/Linux)
# RFdiffusionAA (Stage 2) is cloud-only — see cloud/rfdiffusion_aa/README.md, not installed locally.

# 3. RFdiffusion model weights -> models/weights/
python scripts/download_weights.py
```

**Environment (observed working):** local GPU is an **8 GB RTX 4060 Laptop**; conda + GPU
live in **WSL2** (the agent's shell is git-bash/MINGW). LigandMPNN runs CPU-only here
(trivial at these sizes); AF2-multimer runs under WSL with the `platform` allocator;
RFdiffusionAA is offloaded to a ≥16 GB cloud GPU.

## Scaffold provenance
The Lock scaffold is **PDB 4QAE** — an engineered lipocalin (Anticalin) in complex with
human hepcidin. The current ligand-aware path does **not** depend on its loops being
re-sculpted (it diffuses a fresh pocket around the ligand); 4QAE survives as the carrier
protein for the cloud input and as the source of the frozen Key chain for the instrument.

## Source of truth
**`PROGRESS_LOG.md`** is the authoritative decision log — what was tried, why, and what
didn't work. `CLAUDE.md` is the quick-orientation brief for future agents; both were brought
in line with this README during the 2-Nal pivot.
