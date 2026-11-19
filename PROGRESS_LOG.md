# OrthCom — Ligand-Aware Rebuild: Progress Log

**Purpose.** One place to review every decision in this effort: *what* was chosen,
*why*, and its *effect on the codebase*. Newest entries at the bottom.

- **Branch:** `ligand-aware-rebuild`
- **Started:** 2026-06-06
- **Goal (narrow, methodological):** show a designed receptor pocket can
  **discriminate L-4-fluorophenylalanine (L4F) from canonical PHE** — i.e. an
  *affinity gap*, not merely "binds L4F." Every result before this effort was
  against the PHE stand-in, so it only says the old scaffold doesn't grip the
  helix; it says **nothing yet** about L4F-vs-PHE discrimination.
  - **PIVOT (2026-06-07):** primary target is now **L-2-naphthylalanine (2-Nal)**; **L4F is
    demoted to the negative control**. The thesis is unchanged (discriminate an ncAA from
    canonical residues by an affinity gap), only the ncAA changes. See the 2026-06-07 pivot
    entry in the decision log.
- **Note on tracking:** `data/` is gitignored, so the ligand PDBs below are *not*
  tracked by git — only the generator scripts are. This file is tracked and can be
  committed when you want.

---

## Decision log

### 2026-06-06 — Diagnosis: why the partial-diffusion designs failed
- **Finding.** All prior loop-melt designs show uniformly low Key pLDDT (~26–33).
- **Evidence.** Input dock: Key is a clean α-helix (12/12 helical windows), only 1
  marginal clash (2.08 Å), but shallow/one-sided (44/82 Key atoms buried, COM 13.3 Å
  off pocket center, octant 6/8, the Phe "tooth" touches only 12 scaffold atoms
  <5 Å). AF2 co-folds bury the Key *more* (COM ~9 Å, octant 8/8, 225–550 contacts)
  but **deform the helix** (helical windows drop to 2–8/12; CA(i,i+3) compresses to
  4.3–4.9 Å) and introduce hard clashes (down to 0.5 Å), while Key pLDDT stays ~30.
  Complex pLDDT 68–76 → the Lock folds fine; the **interface** is frustrated.
- **Conclusion.** Combination: **(b)** wrong interface secondary structure (helix vs
  a barrel-loop face) is primary; **(c)** the dock is genuinely too shallow/one-sided;
  **(a)** bad loop sequences are secondary. Root cause traced to a malformed anchor.
- **Effect.** Diagnostic only. A throwaway `_diag_pose.py` was created and deleted
  (no committed change).

### 2026-06-06 — Strategy: full ligand-aware redesign
- **Decision.** Move off partial diffusion to a **ligand-aware** path: RFdiffusionAA
  (backbone around a ligand) + LigandMPNN (sequence conditioned on the ligand).
- **Why.** The PHE stand-in cannot test discrimination — the models never see the
  fluorine. Ligand-aware tools condition on the actual fluoroaromatic.
- **Effect.** Created branch `ligand-aware-rebuild`. No existing scripts/data deleted.

### 2026-06-06 — RFdiffusionAA → cloud offload
- **Decision.** Do **not** install RFdiffusionAA locally; prepare a cloud-ready
  package instead.
- **Why.** All-atom RFdiffusionAA needs ~12 GB+ VRAM and is Apptainer-only; this box
  has 8188 MiB total / ~5758 MiB free. "Moving a broken tooth to a bigger GPU helps
  nothing" — so the anchor is fixed first (see the next entry).
- **Effect.** No local install. Cloud package to be assembled later.

### 2026-06-06 — Decided to fix the malformed anchor before anything downstream
- **Decision.** Rebuild the Key's fluoroaromatic anchor as chemically correct, with
  standard PDB conventions, and confirm it before building downstream.
- **Why.** The frozen anchor in `data/raw/scaffold_clean.pdb` (from
  `src/graft_residue.py`) auto-names atoms by element count (`O1/O2/C3..C9`), keeps a
  spurious 2nd carboxylate O, and keeps explicit hydrogens — chemically wrong for a
  mid-peptide residue and unusable as a clean ligand.
- **Effect.** New generator `src/prep_ligand.py` (does not touch `graft_residue.py`,
  which is left in place).

### 2026-06-06 — Anchor source: PFF CCD ideal coordinates (verified)
- **Decision.** Build from RCSB CCD entry **PFF = 4-fluoro-L-phenylalanine**.
- **Why.** Authoritative atom names + geometry. A WebSearch guess ("F4P") was
  **rejected** after curl-verifying both `.cif` files: PFF is 4-fluoro-L-Phe; F4P is
  an unrelated drug.
- **Effect.** Verified heavy-atom coords embedded in `src/prep_ligand.py`; geometry
  self-check passes (C–F 1.339 Å; para CG···CZ 2.772 Å; F nearest ring C = CZ).

### 2026-06-06 — Anchor representation: Option A (side-chain-centric)
- **Decision.** Represent the ligand as the **side chain + CA cap only**:
  `CA CB CG CD1 CD2 CE1 CE2 CZ F` (9 heavy atoms, **0 oxygens**). Control = same minus
  F (8 atoms). Chosen over B (mid-peptide residue, 1 carbonyl O) and C (full free
  amino acid, 2 carboxylate O).
- **Why.** In the real Key, L4F is residue 8 of 15 — **mid-peptide**, so it has no
  free carboxylate/amino terminus. Including those charges would let RFdiffusionAA
  bury a carboxylate and let that dominate the energetics over the subtle C–F
  discrimination we're actually testing. Zero oxygens also satisfies "no extra O."
- **Effect.** `data/raw/ligand_L4F.pdb` (9 atoms, para-F) and
  `data/raw/ligand_PHE.pdb` (8 atoms, control) — the pair differs by exactly one atom.
  `src/prep_ligand.py --mode residue` retained for later all-atom complex assembly.

### 2026-06-06 — Environment facts (shape the plan)
- My shell is **git-bash (MINGW)**, not WSL; conda + GPU live in
  **WSL2**. GPU: 8188 MiB total / **5758 MiB free**. `conda` is not on the default
  WSL PATH (existing scripts call tools by absolute path).
- **Effect.** Confirms RFdiffusionAA stays cloud-only; LigandMPNN runs locally
  (CPU-capable, light). Conda base to be located at install time.

### 2026-06-06 — LigandMPNN local install approved — IN PROGRESS
- **Decision.** Proceed with clone + new conda env + the <1 GB weights download.
- **Plan.** (1) clone `dauparas/LigandMPNN` into the repo tree (vendored, to be
  gitignored); (2) new env from its `requirements.txt`; (3) `get_model_params.sh`
  weights; (4) dry-run on a dummy backbone + `ligand_L4F.pdb` HETATM; (5) restore
  diversity (`--temperature 0.2 --number_of_batches 8` vs the old
  `sampling_temp=0.0001`); (6) discrimination control (F-present vs F-absent on the
  same backbone) as a first-class output.
- **Effect (so far).**
  - Cloned `dauparas/LigandMPNN` → `LigandMPNN/` (MIT); added `LigandMPNN/` to
    `.gitignore` (vendored, like ProteinMPNN/RFdiffusion).
  - conda env **`ligandmpnn`** (python 3.10) created. First attempt at the pinned
    GPU stack (`torch 2.2.1` + `nvidia-*-cu12`, ~2.5 GB) **failed mid-download**
    (`BrokenPipeError` on the 410 MB cuBLAS wheel — flaky connection, not a real
    dependency error).
  - **Decision: install CPU-only torch instead** (`torch==2.2.1` from the cpu index,
    ~190 MB) + LigandMPNN deps (numpy 1.23.5, scipy, ProDy, biopython, ml-collections,
    dm-tree). **Why:** LigandMPNN inference on our sizes (~180-res Lock + 9-atom
    ligand) is trivial on CPU; this avoids a fragile multi-GB CUDA download for no
    practical benefit. GPU torch can be added later if throughput ever demands it.
    — done.
  - **Env ready (verified):** `torch 2.2.1+cpu`, ProDy 2.4.1, numpy 1.23.5 all import
    OK; `run.py --help` loads (RC=0). One fix needed: ProDy 2.4.1 imports
    `pkg_resources`, which setuptools 82 dropped → pinned `setuptools<81`.
    **Steps 1–3 (clone + env + weights) complete.**
  - Weights via `get_model_params.sh ./model_params` — **done: 118 MB**, 15
    checkpoints incl. `ligandmpnn_v_32_010_25.pt` (the ligand_mpnn default).
  - CLI facts for our use: `--model_type ligand_mpnn` reads ligand HETATM as context;
    **`--ligand_mpnn_use_atom_context 0` ablates it = a clean built-in control**;
    diversity via `--temperature 0.2 --number_of_batches 8`.
  - **Checkpoint.** I stopped to decide here, and went with **A+B**.

### 2026-06-06 — LigandMPNN dry-run (A + B) — DONE
- **Plan.** Part A: install sanity on `inputs/1BC8.pdb`. Part B: build a *throwaway*
  dummy (old Lock chain A + our `ligand_L4F`/`ligand_PHE` placed into the docked
  pocket as HETATM) and design the pocket residues four ways — **F-present vs
  F-absent × atom-context on vs off** — at `--temperature 0.2 --number_of_batches 8
  --save_stats 1 --seed 111`, then compare per-position AA distributions at pocket
  positions. This is the **discrimination-control mechanics** on a throwaway
  backbone, *not* a real design.
- **Effect (results).**
  - **Mechanics verified.** Part A: `run.py` on `inputs/1BC8.pdb` parsed ligand
    context (C/P/Zn) and wrote designs with `ligand_confidence` (RC=0). Part B:
    `src/prep_ligandmpnn_input.py` Kabsch-placed the anchor onto docked PHE8
    (8-atom fit RMSD **0.71 Å**; the two ligands differ by exactly the F atom —
    0.0000 Å over the 8 shared atoms) and found a **6-residue pocket** (A52 A68 A70
    A79 A80 A81; F-proximal: A52 A68 A70 A79). LigandMPNN read **9** ligand atoms
    for L4F (8 C + 1 **F**) and **8** for PHE — the fluorine is seen only in L4F.
  - **Control passes.** With `--ligand_mpnn_use_atom_context 0` the F-present and
    F-absent inputs gave **identical** designs (`L4F_off == PHE_off`, whole chain).
    Same seed → identical decoding order across all 4 conditions → clean A/B test.
  - **Ligand context is strong.** Context ON vs OFF substantially changes the pocket
    sequence (JS divergence up to **0.60 bits** at A79, 0.36 at A68): the
    ligand-aware machinery is functioning (e.g. A79 L/R/K-mix → pure L when the
    ligand is present).
  - **But the fluorine signal is ~negligible.** L4F-on vs PHE-on JS divergence:
    mean **0.0034** bits F-proximal, 0.0016 non-proximal (max **0.0099** at A68).
    At ring-adjacent A52 the model picks Phe deterministically (FFFFFFFF) with or
    without the F. So the single para-F barely moves LigandMPNN's sequence choice —
    **~50–100× smaller than the whole-ligand context effect.**
  - **Methodological takeaway.** LigandMPNN *sequence* design is unlikely to be the
    source of L4F-vs-PHE discrimination — one para-F is too small a perturbation for
    it. Discrimination, if achievable, must come from pocket **shape** (RFdiffusionAA
    building a backbone that complements the F) and/or explicit **physics** affinity
    scoring — not from MPNN AA preferences. *Caveat:* this used the **old shallow**
    scaffold pocket with 0.71 Å placement; a purpose-built F-complementary pocket
    could differ. The direction, not the magnitude, is the takeaway.
  - **Files added:** `src/prep_ligandmpnn_input.py` (place ligand + find pocket),
    `src/analyze_ligandmpnn_dryrun.py` (4-condition comparison). Throwaway I/O under
    `data/processed/ligandmpnn/` and `LigandMPNN/outputs/_dryB/` (both gitignored).
    No existing scripts/data/instrument/physics touched.

### 2026-06-07 — Checkpoint commit (housekeeping)
- **Commit `9d8ed8a`** on `ligand-aware-rebuild` ("Add corrected F-anchor + LigandMPNN
  F-vs-PHE discrimination dry-run"): `src/prep_ligand.py`, `src/prep_ligandmpnn_input.py`,
  `src/analyze_ligandmpnn_dryrun.py`, `PROGRESS_LOG.md`, `.gitignore` (+`LigandMPNN/`).
  651 insertions; diffs shown first. Data/ligand PDBs stay gitignored (untracked).

### 2026-06-07 — RFdiffusionAA: staged for a (small) local-GPU run, not launched
- **Decision.** Pivot RFdiffusionAA from "cloud-only" toward a **small local first run**,
  per your "let's try and get to that point" steer. Stage everything reversible; the heavy
  downloads + system install are held behind your explicit go (below).
- **Repo cloned (code only).** `git clone --recurse-submodules baker-laboratory/rf_diffusion_all_atom`
  → `rf_diffusion_all_atom/` (+ submodule `lib/rf2aa`); added to `.gitignore`. No `.sif`/weights yet.
- **Ligand input VERIFIED against the source — nothing to change.** Traced `make_indep()` →
  `filter_het()` (matches HETATM by `resname == inference.ligand`) → `parse_mol()` which uses
  **OpenBabel**: atom types from the **element column**, bonds from **CONECT/geometry**, and
  **no CCD/template lookup**. So `inference.ligand=PFF` is just a label, and our
  `data/raw/ligand_L4F.pdb` (8 C + para-F, CONECT) drops in unchanged — the para-F is read as
  fluorine; the side-chain-only form (no free termini) is exactly the pocket target we want.
- **Setup probe (sizes are the cost of "local").** Apptainer/Singularity **not installed**
  (and upstream is **Apptainer-only** — no conda/pip path). Container `.sif` = **~11.8 GB**,
  weights = **~1.27 GB** (~13 GB total). WSL disk free 891 GB (fine). VRAM: 8 GB box; NVIDIA's
  ~12 GB floor is for large designs — a **small** de-novo pocket (NRES≈80 + 9 ligand atoms) is
  tiny and expected to fit 8 GB. Hence start small.
- **Design mode (first run): de-novo small pocket**, not motif-scaffold. Why: minimal, directly
  tests "can RFdiffusionAA wrap a backbone around the para-F," and fits 8 GB. Motif-scaffold
  (reuse the beta-barrel framework, rebuild loops around the docked anchor) is the heavier
  follow-up once the toolchain is proven.
- **Effect.** New `src/run_rfdiffusion_aa.sh` (tracked): documents the one-time setup (Apptainer
  + `.sif`/weights URLs and sizes), stages `ligand_L4F.pdb` into the container bind mount, has a
  preflight that fails clearly if setup is missing, and runs the exact upstream command form
  (`contigmap.contigs=['80-80'] inference.ligand=PFF diffuser.T=100`, `--nv`, editable NRES/T/NUM).
  Same command serves a cloud ≥16 GB GPU (just a different host).
- **Before spending anything here I worked out the cost:** (1) `sudo` install Apptainer in WSL2, (2) download the
  ~11.8 GB `.sif`, (3) download the ~1.27 GB weights. These are heavy + reverse the earlier
  cloud-only call. No GPU job will be launched regardless until you say so.

### 2026-06-07 — Decided to stay cloud-only; cloud package finalized
- **Decision.** Do **not** touch the 8 GB box for RFdiffusionAA. Finalize a cloud-ready package;
  also explicitly **confirm the corrected F-ligand from the anchor fix is baked into it**.
  (Also informed by a code finding: `sudo` needs a password, so I couldn't have done the local
  Apptainer install non-interactively anyway.)
- **Carrier-protein finding (changes the input).** Read the source: a **ligand-only** input would
  break `inference/utils.py process_target()` (it builds the contig frame from protein **CA**
  atoms). The upstream examples always pass a real protein + ligand; with a **de-novo** contig
  the carrier protein is parsed then **discarded** by `insert_contig` (only the ligand is kept as
  motif). So the cloud input must be **carrier protein + ligand**, not the bare anchor — my first
  staging script had assumed ligand-only; corrected.
- **Re-checked before packaging.** Re-ran `src/prep_ligand.py`: anchor self-check
  passes (C–F **1.339 Å**, CG···CZ para **2.772 Å**, F nearest ring C = CZ **PARA**, **O = 0**,
  9 atoms). The cloud input's ligand block was machine-checked: resname **PFF**, **1 F**, **0 O**,
  9 heavy atoms → **matches the corrected ligand**. RFdiffusionAA reads it via OpenBabel (element +
  CONECT, no CCD lookup), so the para-F is seen as fluorine; nothing to edit.
- **Effect (files).**
  - `src/prep_rfdiffusion_aa_input.py` (new, tracked): builds the cloud input by reusing the
    verified Kabsch placement; bakes the corrected L4F anchor (resname PFF, CONECT) onto scaffold
    chain A; asserts 1 F / 0 O / 9 atoms.
  - `cloud/rfdiffusion_aa/input_L4F_pocket.pdb` (new, **tracked** — the baked-in package input).
  - `cloud/rfdiffusion_aa/README.md` (new, tracked): self-contained setup (clone + `.sif` +
    weights + Apptainer), the exact copy-paste run command, the verification stamp, tunables, and
    the downstream discrimination plan (LigandMPNN → AF2 co-fold → physics L4F-vs-PHE gap; control
    = `data/raw/ligand_PHE.pdb`).
  - `src/run_rfdiffusion_aa.sh` (updated): reframed cloud-only; now uses the carrier+ligand input
    and `inference.ligand=PFF`, default `contigmap.contigs=['120-120']` (no longer the 8 GB-driven
    `['80-80']`).
  - `rf_diffusion_all_atom/` stays cloned-but-light (no `.sif`/weights downloaded locally).

### 2026-06-07 — 2-Nal ligand anchor built and verified  *(completed)*
- **Decision.** Build `data/raw/ligand_2NAL.pdb` — the 2-Nal Option-A anchor — from verified
  RCSB CCD NAL ideal coordinates. Standalone `src/prep_ligand_2nal.py`; does **not** touch
  `src/prep_ligand.py` or the L4F/PHE controls (confirmed byte-identical by `git hash-object`
  before and after: `ligand_L4F.pdb` = `ee44ce78`, `ligand_PHE.pdb` = `fc4bdc78`).
- **Why.** 2-Nal is the primary ncAA target; the ligand anchor is the required input for
  RFdiffusionAA pocket design. CCD code NAL confirmed: "BETA-(2-NAPHTHYL)-ALANINE", L-peptide,
  (2S) stereocenter, naphthalen-2-yl attachment via CB methylene.
- **Asserts passed (machine-checked):** 12 atoms (all C) · O=0, N=0, F=0 · 11 ring C-C bonds
  (naphthalene) · degree-3 fusion carbons = {C4A, C8A}, mutually bonded · CB↔{CA, C2} only ·
  CA↔{CB} only · ring planarity max-dev 0.005 Å · sha1 = `ac9f19f7...` (raw) / git =
  `f1f6ede6...`.
- **Effect.** `src/prep_ligand_2nal.py` (tracked, reproduces `ligand_2NAL.pdb` byte-for-byte).
  Scratch `data/processed/_scratch_2nal.cif` deleted. `data/raw/ligand_2NAL.pdb` remains
  gitignored (generated).
- **Checked.** I confirmed the atom names, ring connectivity and O count of 0 before committing.

### 2026-06-07 — 2-Nal pivot: rewrote the README and my notes
- **Decision (yours).** Pivot the **primary target** from 4-fluoro-Phe (L4F) to
  **L-2-naphthylalanine (2-Nal)** — a fused bicyclic aromatic with genuine shape/size — and
  **demote L4F to the negative control / baseline** (the documented ncAA too subtle to
  discriminate). Thesis unchanged: discriminate an ncAA from canonical residues by an
  **affinity gap**. New work branch **`ncaa-naphthyl`**.
- **Why.** The LigandMPNN dry-run: a single para-F moved sequence choice ~0.003 bits
  (≈ noise; ~50–100× below whole-ligand context). One fluorine is too small a feature; a
  bicyclic side chain gives a pocket real geometry to grip and no canonical look-alike.
- **Effect — documentation overhaul (this commit).**
  - `README.md`: full rewrite. Dropped the "orthogonal communication axis / platform"
    framing; states the narrowed discrimination thesis, the **ligand-aware** pipeline
    (anchor → RFdiffusionAA pocket [cloud-only] → LigandMPNN → AF2-multimer co-fold
    instrument → physics affinity gap), the 2-Nal/L4F/PHE roles, a verified repo map, and a
    "What didn't work" section. Fixed a stale path ref (`graft_dummmy_residue.py`, 3 m's).
  - My notes file: full rewrite mirroring the README. Marked the **partial-diffusion** era and
    its contig mechanics **SUPERSEDED**; replaced "assume the partial-diffusion strategy"
    with "assume the ligand-aware strategy"; refreshed the artifact list and the still-true
    hardware/cloud-only facts; points here as source of truth.
  - This log: added the **"What didn't work (and why)"** section below and the
    **2-Nal-vs-Trp** action item under *Next*.
  - No scripts/data touched; instrument and physics path untouched; deprecated scripts kept.
- **Checked.** I read both files in full before committing.

---

## What didn't work (and why)

Institutional memory — keep it. Newest first. (Detailed evidence is in the dated decision-log
entries above; this is the durable index, mirrored in `README.md`.)

1. **A single fluorine is not a discriminable feature (for sequence design).** LigandMPNN with
   ligand atom-context ON moved its sequence choice between L4F and PHE by Jensen–Shannon
   divergence **mean ~0.0034 bits** F-proximal (max 0.0099) vs up to **~0.60 bits** for
   whole-ligand context — ~50–100× smaller. → MPNN *sequence* design cannot discriminate one
   para-F; discrimination must come from pocket **shape** (RFdiffusionAA) and/or **physics**.
   Hence **L4F → negative control** and the primary target → **2-Nal**.
2. **Helical Key in a barrel = frustrated interface (era 2, partial diffusion).** Shallow,
   one-sided dock (44/82 Key atoms buried, COM ~13 Å off centre); AF2 co-folds buried it more
   but **deformed the helix** (helical windows 12→2–8), introduced hard clashes (to 0.5 Å),
   Key pLDDT stuck ~30 while the Lock folded fine. Helix-vs-barrel-loop mismatch; `partial_T`
   cannot change loop length or relieve a hard clash. → abandoned partial diffusion.
3. **The original anchor was chemically malformed.** `graft_residue.py` auto-named atoms by
   element count (`O1/O2/C3…C9`), kept a spurious 2nd carboxylate O and explicit H — wrong for
   a mid-peptide residue, unusable as a clean ligand. → rebuilt by `prep_ligand.py` from
   verified CCD coords (side-chain-centric, 0 free O).
4. **De-novo hallucination = "pool-noodle" helices (era 1).** Hallucinating 100–200 residues
   around the fixed Key gave long α-helices, not pockets (~1000 designs to fish for a
   concavity). → motivated the scaffold pivot (era 2), then the ligand-aware approach (current).

*Smaller gotchas:* DL tools crash on `L4F` (`KeyError`) → PHE "Trojan-Horse" stand-in for DL
steps; Tyr rejected (hydroxyl polar penalty repels the swap-back F); AF2's BFC allocator
segfaults at recycle 2 → JAX `platform` allocator; synthetic proteins have no MSA →
`--msa-mode single_sequence`.

---

## Preserved / out of scope (guardrails for this effort)
- Physics path (Psi4 / RESP / AmberTools, `generate_conformer.py`, `graft_residue.py`)
  — **not touched**.
- Working validation instrument (`src/parse_interface.py`, `src/run_multimer.sh`,
  `src/prep_multimer.py`) — **not touched**.
- No existing scripts/data/outputs deleted; deprecate by renaming, not removing.
- Long GPU jobs require showing the command and getting a yes first.

## Next
- **Active effort — 2026-06-07 pivot (branch `ncaa-naphthyl`):** rewrote the README and my
  notes ✔; built the 2-Nal anchor ✔; re-pointing the scripts and cloud package at 2-Nal
  **paused, see the next bullet**.
- **Paused (2026-06-07).** I stopped to weigh up a reframe toward
  immunogenicity screening of ncAA-bearing peptides (structure-based ncAA-pMHC scoring
  vs. the canonical-surrogate shortcut). The verified 2-Nal anchor and the physics path
  carry over to that reframe; RFdiffusionAA pocket design is shelved, not deleted,
  pending that decision.
- **KNOWN PENDING (not done):** `cloud/rfdiffusion_aa/` still contains only
  `input_L4F_pocket.pdb` (the L4F build). The **2-Nal sibling input + the README re-point
  are still to do** — the cloud package is not re-pointed at 2-Nal yet.
- **ACTION ITEM — 2-Nal-vs-Trp cross-reactivity test (decisive).** Trp is the only
  fused-bicyclic canonical aromatic, so it — not PHE — is the real cross-reactivity threat to
  the 2-Nal story. Once a pocket exists, run the **2-Nal-vs-Trp** affinity gap explicitly; it
  is the single test most likely to falsify the thesis. On the books now.
- **Strategic note from the dry-run:** the discriminator is not LigandMPNN sequence
  design. Weight the plan toward (a) RFdiffusionAA pocket *shape* that complements
  the para-F, and (b) the physics affinity-gap scoring as the actual L4F-vs-PHE test.
- **Cloud package** (`cloud/rfdiffusion_aa/` — currently the **L4F** build; being re-pointed to **2-Nal**): take it to a ≥16 GB GPU host,
  do the one-time setup (clone + `.sif` + weights + Apptainer), run the command in the
  README. First design is `NRES=120, NUM=1`; raise `num_designs` to fan out.
- After a backbone exists: LigandMPNN sequence (env already installed) → AF2 co-fold
  instrument → physics affinity-gap (**L4F vs PHE**; control `data/raw/ligand_PHE.pdb`)
  = the real discrimination test.
- Possible follow-up: a **motif-scaffold** variant (keep beta-barrel strands, rebuild
  loops around the docked anchor) if de-novo pockets don't grip the aromatic.


---

## 2026-11-19 - Cloud package re-pointed at 2-Nal; picking the project back up

**Resuming.** The 2026-06-07 pause stopped the re-pointing work while a reframe toward
immunogenicity screening was evaluated. That reframe exists as its own live repo
(`orthcom-immuno`, which has measured Step 2 results). Decision: OrthCom continues on
2-Nal as a separate active project rather than shelving it. The "paused" and
"RFdiffusionAA pocket design is shelved" lines above are superseded. `orthcom-immuno`
said OrthCom was shelved; I corrected its README and notes to call it active and
to frame the no-DL-pipeline rule as a scope boundary between two live projects.

### What / why / effect

**1. Measured why L4F failed, in volume rather than in bits.**
New `src/analyze_shape_gap.py` (CPU, no weights) measures vdW volume on a grid after
superposing each anchor on the canonical it must be told apart from. PHE 97.6 A^3,
L4F 104.3, 2-Nal 136.6. The important column is the PHE-shaped *void*: 0.0 A^3 for L4F
and 1.0 A^3 for 2-Nal. PHE fits entirely inside both ncAA envelopes, so a pocket
complementary to either still admits PHE with no steric penalty. The L4F signal was
~0.2 to 0.3 kcal/mol, which is why 0.0034 bits was not a model failure. Effect: the
project can no longer claim discrimination against PHE on steric grounds.

**2. Built the Trp anchor and measured the comparison that matters.**
New `src/prep_ligand_trp.py` builds `data/raw/ligand_TRP.pdb` from RCSB CCD entry TRP
(TRYPTOPHAN, C11 H12 N2 O2; source cif sha256 4cea7e28...6e36d6), Option A, 11 atoms =
10 C + 1 N, 0 O, sha1-pinned 4b7ba79f447896636da9bd6ade07cd4b0cdd4db5, asserts fusion at
CD2/CE2 and NE1 present.

`analyze_shape_gap.py` now fixes CA/CB/ipso and searches chi2 for the Trp rotamer that
MAXIMISES overlap with 2-Nal (Trp's best case, the conservative direction) because indole
is 6-5 fused where naphthalene is 6-6 and there is no atom-to-atom ring correspondence to
Kabsch-fit. At chi2 180 deg: shared 103.9 A^3 = **76% of 2-Nal**, Trp surplus 19.7 A^3,
2-Nal void 32.5 A^3, and **NE1 lands 0.56 A from a naphthalene carbon position, inside
the 2-Nal envelope**. Across the scan overlap ranges 72.6 to 103.9 A^3.

Effect, and this is the main outcome: **the discrimination lever is the buried unsatisfied
N-H donor, not volume.** A pocket selected for being large admits Trp at 76% overlap and
discriminates nothing, which is the L4F mistake in a different form.

**3. Cloud package re-pointed at 2-Nal.**
`src/prep_rfdiffusion_aa_input.py` generalized with `--residue {2NAL,L4F}`, default 2NAL.
New `cloud/rfdiffusion_aa/input_2NAL_pocket.pdb` (carrier chain A 1391 atoms + NAL
ligand, 12 C, 0 O, 11 ring bonds, kabsch RMSD 0.720 A). `--residue L4F` reproduces the
tracked `input_L4F_pocket.pdb` byte-for-byte, verified by an empty `git diff`, so the
baseline stays reproducible. The 2-Nal correspondence composes NAL->PHE with the existing
PHE->docked map, since the docked anchor carries the malformed `C3..C9` names from
`graft_residue.py`; its RMSD of 0.720 A against L4F's 0.707 A on the same distorted
docked ring is the evidence the map is right.

The builder now reports carrier clashes. **L4F: 0 contacts under 3.0 A, min 3.14 A.
2-Nal: 8 contacts, min 1.76 A.** Harmless for the de-novo contig (the carrier is
discarded by `insert_contig`) but it settles an open question: the motif-scaffold variant
listed as a possible follow-up is **not available on 4QAE for 2-Nal**. De novo is the
only route. Contig raised from 120 to 130 because 2-Nal reaches 6.46 A from CB against
PHE's 4.28 A.

**4. The N1 criterion is a filter, not an objective.**
RFdiffusionAA conditions on ligand atoms and has no term for "apolar here", so the lever
cannot be passed on the command line. New `src/check_n1_pocket.py` applies it per design:
place the best-mimicking Trp, then require NE1 to be buried (at least 8 protein heavy
atoms within 5.0 A) AND no H-bond acceptor within 3.5 A. Exits 0 on PASS, 1 on FAIL, so
it can gate a loop. Run on the cloud input it FAILS with 3 atoms within 5 A of NE1, which
is correct: 4QAE was never designed around 2-Nal. Keep that as the filter's negative
control. Its shared-volume figure (104.0 A^3 at 0.2 A spacing) agrees with
`analyze_shape_gap.py` (103.9 at 0.15 A), which cross-checks both.

### Still not done
- No cloud run launched. No pocket exists, so no affinity gap for any target.
- The Trp result is shape, not affinity. The 1 to 3 kcal/mol donor penalty is a
  literature range for buried unsatisfied donors, not a calculation on this system.
- 2-Nal is held at its CCD ideal chi2 throughout; its own rotamer freedom is unexplored.
- Physics path still exports no ESP grid and runs no AmberTools RESP step.
