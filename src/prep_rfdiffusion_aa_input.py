#!/usr/bin/env python3
"""Build the RFdiffusionAA cloud input: carrier protein + a verified ncAA anchor.

WHY a carrier protein (not the bare ligand)
-------------------------------------------
RFdiffusionAA's inference/utils.py `process_target()` builds the contig frame
from protein ATOM/CA records (and centers on the CA centroid); a ligand-only PDB
has no CA and breaks it. For a DE-NOVO contig (e.g. ['120-120'], no chain motif
refs) the carrier protein is parsed and then DISCARDED by `insert_contig` --
only the ligand (is_sm) is kept as the motif the fresh pocket is built around.
This mirrors the upstream README example, which passes the full 7v11 protein with
a ['150-150'] contig to design a de-novo binder to the OQO ligand.

We reuse scaffold_clean.pdb chain A (a real 172-res protein) as the carrier and
place the anchor in its pocket. For a de-novo contig the placement is cosmetic,
since the carrier is discarded; it matters only if someone later attempts a
motif-scaffold run, which is why this script reports carrier clashes instead of
staying silent about them.

TARGETS (--residue)
-------------------
  2NAL  L-2-naphthylalanine, the PRIMARY target. resname NAL, 12 carbons.
  L4F   4-fluoro-phenylalanine, now the NEGATIVE CONTROL / baseline. resname PFF.

Both anchors come from src/prep_ligand_2nal.py and src/prep_ligand.py, which build
from verified RCSB CCD ideal coordinates and self-check their geometry. OpenBabel
reads element + CONECT with no CCD lookup (verified in
rf_diffusion_all_atom/aa_model.py parse_mol), so the resname is only a label that
has to match `inference.ligand` on the run command.

OUTPUTS (tracked, so the verified ligand is baked into the cloud package)
  cloud/rfdiffusion_aa/input_2NAL_pocket.pdb
  cloud/rfdiffusion_aa/input_L4F_pocket.pdb
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from prep_ligandmpnn_input import (parse_pdb_atoms, kabsch, het_line,
                                   REPO, SCAFFOLD, LIG_L4F, CORRESP)

OUTDIR = os.path.join(REPO, "cloud", "rfdiffusion_aa")
LIGCHAIN = "B"
LIGRESSEQ = 1
CLASH_CUTOFF = 3.0      # A, heavy-atom contact closer than this is a clash

# The docked anchor in scaffold_clean.pdb chain B resseq 8 carries the MALFORMED
# auto-generated atom names from graft_residue.py (CA, C3..C9), not standard PHE
# names. CORRESP (imported) maps L4F atom names onto those. The 2-Nal map is the
# composition of NAL->PHE with that same PHE->docked map:
#   NAL   PHE    docked
#   CA    CA     CA
#   CB    CB     C3
#   C2    CG     C4     (ipso, bonded to CB)
#   C3    CD1    C5
#   C1    CD2    C9
#   C4    CE1    C6
#   C8A   CE2    C8
#   C4A   CZ     C7     (para; the fused ring grows off C4A/C8A)
CORRESP_2NAL = [
    ("CA", "CA"),
    ("CB", "C3"),
    ("C2", "C4"),
    ("C3", "C5"),
    ("C1", "C9"),
    ("C4", "C6"),
    ("C8A", "C8"),
    ("C4A", "C7"),
]

# Side-chain connectivity, subsets of each CCD graph.
BONDS_L4F = [
    ("CA", "CB"), ("CB", "CG"),
    ("CG", "CD1"), ("CG", "CD2"),
    ("CD1", "CE1"), ("CD2", "CE2"),
    ("CE1", "CZ"), ("CE2", "CZ"),
    ("CZ", "F"),
]
BONDS_2NAL = [
    ("CA", "CB"), ("CB", "C2"),
    ("C2", "C1"), ("C2", "C3"),
    ("C3", "C4"), ("C4", "C4A"),
    ("C1", "C8A"), ("C4A", "C8A"),
    ("C4A", "C5"), ("C5", "C6"),
    ("C6", "C7"), ("C7", "C8"),
    ("C8", "C8A"),
]

SPECS = {
    "2NAL": {
        "label": "2-Nal",
        "ligfile": os.path.join(REPO, "data", "raw", "ligand_2NAL.pdb"),
        "ligcode": "NAL",
        "out": "input_2NAL_pocket.pdb",
        "corresp": CORRESP_2NAL,
        "bonds": BONDS_2NAL,
        "n_atoms": 12,
        "elements": {"C": 12},
        "ring_bonds": 11,
        "nonring": ("CA", "CB"),
    },
    "L4F": {
        "label": "L4F",
        "ligfile": LIG_L4F,
        "ligcode": "PFF",
        "out": "input_L4F_pocket.pdb",
        "corresp": CORRESP,
        "bonds": BONDS_L4F,
        "n_atoms": 9,
        "elements": {"C": 8, "F": 1},
        "ring_bonds": 6,
        # F hangs off the ring rather than closing it, so it is not a ring bond.
        "nonring": ("CA", "CB", "F"),
    },
}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--residue", choices=sorted(SPECS), default="2NAL",
                    help="which anchor to bake in (default 2NAL, the primary target)")
    args = ap.parse_args()
    spec = SPECS[args.residue]
    out = os.path.join(OUTDIR, spec["out"])

    os.makedirs(OUTDIR, exist_ok=True)
    if not os.path.exists(spec["ligfile"]):
        sys.exit(f"missing {spec['ligfile']}\n  build it first: "
                 f"python src/prep_ligand{'_2nal' if args.residue == '2NAL' else ''}.py")

    scaf = parse_pdb_atoms(SCAFFOLD)
    chainA = [a for a in scaf if a["chain"] == "A" and a["record"] == "ATOM"]
    docked = {a["name"]: a["xyz"] for a in scaf
              if a["chain"] == "B" and a["resseq"] == 8}

    lig = parse_pdb_atoms(spec["ligfile"])
    lig_by_name = {a["name"]: a for a in lig}
    missing = [ln for ln, _ in spec["corresp"] if ln not in lig_by_name]
    missing += [pn for _, pn in spec["corresp"] if pn not in docked]
    if missing:
        sys.exit(f"correspondence refers to absent atoms: {missing}")

    P = np.array([lig_by_name[ln]["xyz"] for ln, _ in spec["corresp"]])
    Q = np.array([docked[pn] for _, pn in spec["corresp"]])
    R, t, rmsd = kabsch(P, Q)
    print(f"[kabsch] {spec['label']} anchor -> pocket on {len(spec['corresp'])} "
          f"atoms, RMSD = {rmsd:.3f} A")
    if rmsd > 0.5:
        print("[kabsch] WARNING: RMSD > 0.5 A, check the ring correspondence",
              file=sys.stderr)
    placed = [{**a, "xyz": (R @ a["xyz"]) + t} for a in lig]

    # serials: ligand HETATM continue after the protein to avoid id collisions
    base = len(chainA)
    serial = {a["name"]: base + i + 1 for i, a in enumerate(placed)}
    names = [a["name"] for a in placed]

    with open(out, "w", newline="\n") as fh:
        fh.write("REMARK  RFdiffusionAA cloud input: carrier protein (scaffold chain A,\n")
        fh.write(f"REMARK  discarded by a de-novo contig) + corrected {spec['label']} "
                 f"anchor (resname\n")
        fh.write(f"REMARK  {spec['ligcode']}). Generated by "
                 f"src/prep_rfdiffusion_aa_input.py.\n")
        for a in chainA:
            fh.write(a["raw"] + "\n")
        fh.write("TER\n")
        for a in placed:
            fh.write(het_line(serial[a["name"]], a["name"], spec["ligcode"], LIGCHAIN,
                              LIGRESSEQ, a["xyz"], a["element"]))
        present = set(names)
        for n in names:
            neigh = [b for (x, y) in spec["bonds"]
                     for (a_, b) in [(x, y), (y, x)]
                     if a_ == n and b in present]
            if neigh:
                refs = "".join(f"{serial[m]:>5}" for m in neigh)
                fh.write(f"CONECT{serial[n]:>5}{refs}\n")
        fh.write("END\n")

    counts = {}
    for a in placed:
        counts[a["element"]] = counts.get(a["element"], 0) + 1
    print(f"[write] {out}")
    print(f"        carrier: chain A, {len(chainA)} protein atoms")
    print(f"        ligand : resname {spec['ligcode']}, {len(placed)} atoms, "
          f"elements {counts}")
    print(f"        atoms  : {names}")

    assert len(placed) == spec["n_atoms"], \
        f"expected {spec['n_atoms']} ligand atoms, got {len(placed)}"
    assert counts == spec["elements"], \
        f"expected elements {spec['elements']}, got {counts}"
    assert counts.get("O", 0) == 0, "side-chain-centric anchor must have zero oxygens"
    nonring = spec["nonring"]
    ring = {(min(a, b), max(a, b)) for a, b in spec["bonds"]
            if a not in nonring and b not in nonring}
    assert len(ring) == spec["ring_bonds"], \
        f"expected {spec['ring_bonds']} ring bonds, got {len(ring)}"
    print(f"        CHECK: {spec['n_atoms']} atoms, elements {spec['elements']}, "
          f"0 O, {spec['ring_bonds']} ring bonds  PASS")

    # Carrier clashes. Harmless for a de-novo contig (the carrier is discarded),
    # but a motif-scaffold run would have to resolve them, so report rather than hide.
    lig_xyz = np.array([a["xyz"] for a in placed])
    car_xyz = np.array([a["xyz"] for a in chainA])
    d = np.linalg.norm(car_xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    n_clash = int((d < CLASH_CUTOFF).sum())
    print(f"        carrier contacts closer than {CLASH_CUTOFF} A: {n_clash} "
          f"(min {d.min():.2f} A)")
    if n_clash:
        print(f"        NOTE: {spec['label']} overlaps the 4QAE carrier. This is "
              f"EXPECTED and harmless")
        print("              for the de-novo contig below, where insert_contig keeps only")
        print("              the ligand. It does block a motif-scaffold run on this "
              "scaffold.")


if __name__ == "__main__":
    main()
