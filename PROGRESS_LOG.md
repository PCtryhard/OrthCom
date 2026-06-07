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

### 2026-06-06 — Strategy: full ligand-aware redesign  *(your decision)*
- **Decision.** Move off partial diffusion to a **ligand-aware** path: RFdiffusionAA
  (backbone around a ligand) + LigandMPNN (sequence conditioned on the ligand).
- **Why.** The PHE stand-in cannot test discrimination — the models never see the
  fluorine. Ligand-aware tools condition on the actual fluoroaromatic.
- **Effect.** Created branch `ligand-aware-rebuild`. No existing scripts/data deleted.

### 2026-06-06 — RFdiffusionAA → cloud offload  *(your decision)*
- **Decision.** Do **not** install RFdiffusionAA locally; prepare a cloud-ready
  package instead.
- **Why.** All-atom RFdiffusionAA needs ~12 GB+ VRAM and is Apptainer-only; this box
  has 8188 MiB total / ~5758 MiB free. "Moving a broken tooth to a bigger GPU helps
  nothing" — so the anchor is fixed first (GATE below).
- **Effect.** No local install. Cloud package to be assembled later.

### 2026-06-06 — GATE: fix the malformed anchor before anything downstream  *(your decision)*
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

### 2026-06-06 — Anchor representation: Option A (side-chain-centric)  *(your decision)*
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
- The agent's Bash tool is **git-bash (MINGW)**, not WSL; conda + GPU live in
  **WSL2**. GPU: 8188 MiB total / **5758 MiB free**. `conda` is not on the default
  WSL PATH (existing scripts call tools by absolute path).
- **Effect.** Confirms RFdiffusionAA stays cloud-only; LigandMPNN runs locally
  (CPU-capable, light). Conda base to be located at install time.

### 2026-06-06 — LigandMPNN local install approved  *(your decision)* — IN PROGRESS
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
  - **Gate:** stopped for OK; you approved **A+B**.

### 2026-06-06 — LigandMPNN dry-run (A + B)  *(your decision: A+B)* — DONE
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

### 2026-06-07 — Checkpoint commit + permission-audit (housekeeping)
- **Commit `9d8ed8a`** on `ligand-aware-rebuild` ("Add corrected F-anchor + LigandMPNN
  F-vs-PHE discrimination dry-run"): `src/prep_ligand.py`, `src/prep_ligandmpnn_input.py`,
  `src/analyze_ligandmpnn_dryrun.py`, `PROGRESS_LOG.md`, `.gitignore` (+`LigandMPNN/`).
  651 insertions; diffs shown first. Data/ligand PDBs stay gitignored (untracked).
- **`/fewer-permission-prompts`:** scanned all 6 session transcripts. **No changes made** —
  every ≥3× command is either auto-allowed (`cd`/`echo`/`ls`/`git` reads), arbitrary-exec
  (`wsl`/`python`/`bash`, unsafe to allowlist), or mutating (`rm`). The WSL-wrapper workflow
  can't be safely allowlisted without a blanket `Bash(wsl -e bash -lc *)` (= arbitrary exec),
  which I declined to add. `.claude/settings.json` untouched.

### 2026-06-07 — RFdiffusionAA: staged for a (small) local-GPU run — GATE PENDING  *(your steer: get to the run point, don't launch)*
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
- **GATE (awaiting your go):** before I (1) `sudo` install Apptainer in WSL2, (2) download the
  ~11.8 GB `.sif`, (3) download the ~1.27 GB weights. These are heavy + reverse the earlier
  cloud-only call. No GPU job will be launched regardless until you say so.

---

## Preserved / out of scope (guardrails for this effort)
- Physics path (Psi4 / RESP / AmberTools, `generate_conformer.py`, `graft_residue.py`)
  — **not touched**.
- Working validation instrument (`src/parse_interface.py`, `src/run_multimer.sh`,
  `src/prep_multimer.py`) — **not touched**.
- No existing scripts/data/outputs deleted; deprecate by renaming, not removing.
- Long GPU jobs require showing the command and getting a yes first.

## Next
- **Strategic note from the dry-run:** the discriminator is not LigandMPNN sequence
  design. Weight the plan toward (a) RFdiffusionAA pocket *shape* that complements
  the para-F, and (b) the physics affinity-gap scoring as the actual L4F-vs-PHE test.
- **RFdiffusionAA is staged** (`src/run_rfdiffusion_aa.sh`; ligand verified drop-in).
  **GATE — awaiting your go** on the ~13 GB + Apptainer install (see entry above):
  1) `sudo` install Apptainer in WSL2, 2) download `.sif` (~11.8 GB), 3) download
  weights (~1.27 GB). Then a small first run (`NRES=80, NUM=1`) — shown before launch.
- After a backbone exists: LigandMPNN sequence (env already installed) → AF2 co-fold
  instrument → physics affinity-gap (L4F vs PHE) = the real discrimination test.
