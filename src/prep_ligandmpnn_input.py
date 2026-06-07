#!/usr/bin/env python3
"""
prep_ligandmpnn_input.py  —  build a *throwaway* LigandMPNN dry-run input.

Purpose (discrimination-control mechanics, NOT a real design):
  Take the Lock (chain A of data/raw/scaffold_clean.pdb) and drop our corrected
  fluoroaromatic anchor into the docked pocket as a HETATM ligand, so LigandMPNN
  reads it as atom context. Emit two inputs that differ by exactly one atom:
    - dummy_L4F.pdb : ligand WITH the para-F   (the discriminating ligand)
    - dummy_PHE.pdb : ligand WITHOUT the F     (the F-absent control)

  The anchor is placed by rigid (Kabsch) superposition of the ligand's shared
  heavy atoms onto the side chain of the docked PHE8 in chain B. The two ligands
  share identical coords for the 8 shared atoms, so after the SAME transform the
  only difference between the two dummies is the single F row.

  Ligand resname is set to LIG (non-standard) so the PDB parser treats every
  ligand atom as ligand context, never as a protein residue.

Outputs (gitignored throwaway): data/processed/ligandmpnn/
  dummy_L4F.pdb, dummy_PHE.pdb, pocket_residues.txt
"""
import os
import sys
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCAFFOLD = os.path.join(REPO, "data", "raw", "scaffold_clean.pdb")
LIG_L4F = os.path.join(REPO, "data", "raw", "ligand_L4F.pdb")
LIG_PHE = os.path.join(REPO, "data", "raw", "ligand_PHE.pdb")
OUTDIR = os.path.join(REPO, "data", "processed", "ligandmpnn")

POCKET_CUTOFF = 6.0   # A, chain-A residue within this of any ligand atom
FPROX_CUTOFF = 5.0    # A, chain-A residue within this of the para-F

# ligand atom name  ->  docked PHE8 atom name  (CA, CB, ring; F excluded from fit)
CORRESP = [
    ("CA", "CA"),
    ("CB", "C3"),
    ("CG", "C4"),
    ("CD1", "C5"),
    ("CE1", "C6"),
    ("CZ", "C7"),
    ("CE2", "C8"),
    ("CD2", "C9"),
]


def parse_pdb_atoms(path):
    """Return list of dicts for ATOM/HETATM lines, preserving the raw line."""
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append({
                    "record": line[0:6].strip(),
                    "name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "resseq": int(line[22:26]),
                    "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                    "element": line[76:78].strip(),
                    "raw": line.rstrip("\n"),
                })
    return atoms


def kabsch(P, Q):
    """Rotation+translation mapping P onto Q (both Nx3). Returns R, t, rmsd."""
    Pc, Qc = P.mean(0), Q.mean(0)
    A = (P - Pc).T @ (Q - Qc)
    V, _, Wt = np.linalg.svd(A)
    d = np.sign(np.linalg.det(Wt.T @ V.T))
    D = np.diag([1.0, 1.0, d])
    R = Wt.T @ D @ V.T
    t = Qc - R @ Pc
    rmsd = np.sqrt(np.mean(np.sum(((R @ P.T).T + t - Q) ** 2, axis=1)))
    return R, t, rmsd


def het_line(serial, name, resname, chain, resseq, xyz, element):
    name4 = name.ljust(4) if len(name) == 4 else (" " + name).ljust(4)
    x, y, z = xyz
    return (f"HETATM{serial:5d} {name4} {resname:>3s} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  1.00          {element:>2s}\n")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    scaf = parse_pdb_atoms(SCAFFOLD)
    chainA = [a for a in scaf if a["chain"] == "A" and a["record"] == "ATOM"]
    phe8 = {a["name"]: a["xyz"] for a in scaf
            if a["chain"] == "B" and a["resseq"] == 8}
    print(f"[parse] chain A atoms: {len(chainA)}  | PHE8 atoms: {len(phe8)}")

    lig_l4f = parse_pdb_atoms(LIG_L4F)
    lig_phe = parse_pdb_atoms(LIG_PHE)
    lig_by_name = {a["name"]: a for a in lig_l4f}

    # build Kabsch correspondence
    P = np.array([lig_by_name[ln]["xyz"] for ln, _ in CORRESP])   # ligand frame
    Q = np.array([phe8[pn] for _, pn in CORRESP])                  # docked frame
    R, t, rmsd = kabsch(P, Q)
    print(f"[kabsch] superpose ligand->PHE8 on {len(CORRESP)} atoms, RMSD = {rmsd:.3f} A")
    if rmsd > 0.5:
        print("[kabsch] WARNING: RMSD > 0.5 A — check ring correspondence", file=sys.stderr)

    def place(atoms):
        out = []
        for a in atoms:
            out.append({**a, "xyz": (R @ a["xyz"]) + t})
        return out

    placed_l4f = place(lig_l4f)
    placed_phe = place(lig_phe)

    # sanity: shared atoms identical between the two placed ligands
    shared = ["CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    pl = {a["name"]: a["xyz"] for a in placed_l4f}
    pp = {a["name"]: a["xyz"] for a in placed_phe}
    maxdev = max(np.linalg.norm(pl[n] - pp[n]) for n in shared)
    print(f"[check] max coord deviation on 8 shared atoms (L4F vs PHE): {maxdev:.4f} A")

    def write_dummy(path, placed):
        with open(path, "w") as fh:
            fh.write("REMARK  throwaway LigandMPNN dry-run input (src/prep_ligandmpnn_input.py)\n")
            for a in chainA:
                fh.write(a["raw"] + "\n")
            fh.write("TER\n")
            for i, a in enumerate(placed, start=1):
                fh.write(het_line(i, a["name"], "LIG", "L", 1, a["xyz"], a["element"]))
            fh.write("END\n")

    out_l4f = os.path.join(OUTDIR, "dummy_L4F.pdb")
    out_phe = os.path.join(OUTDIR, "dummy_PHE.pdb")
    write_dummy(out_l4f, placed_l4f)
    write_dummy(out_phe, placed_phe)
    print(f"[write] {out_l4f}  (chain A {len(chainA)} atoms + {len(placed_l4f)} HETATM)")
    print(f"[write] {out_phe}  (chain A {len(chainA)} atoms + {len(placed_phe)} HETATM)")

    # pocket residues: chain-A residues near any placed ligand atom
    lig_xyz = np.array([a["xyz"] for a in placed_l4f])
    f_xyz = pl["F"]
    res_atoms = {}
    for a in chainA:
        res_atoms.setdefault(a["resseq"], []).append(a["xyz"])
    pocket, fprox = [], []
    for rs in sorted(res_atoms):
        coords = np.array(res_atoms[rs])
        dmin = np.min(np.linalg.norm(coords[:, None, :] - lig_xyz[None, :, :], axis=2))
        if dmin <= POCKET_CUTOFF:
            pocket.append(rs)
            dF = np.min(np.linalg.norm(coords - f_xyz[None, :], axis=1))
            if dF <= FPROX_CUTOFF:
                fprox.append(rs)

    pocket_str = " ".join(f"A{r}" for r in pocket)
    fprox_str = " ".join(f"A{r}" for r in fprox)
    print(f"[pocket] {len(pocket)} residues within {POCKET_CUTOFF} A of ligand:")
    print(f"         {pocket_str}")
    print(f"[F-prox] {len(fprox)} residues within {FPROX_CUTOFF} A of the para-F:")
    print(f"         {fprox_str}")

    with open(os.path.join(OUTDIR, "pocket_residues.txt"), "w") as fh:
        fh.write(pocket_str + "\n")
        fh.write("# F-proximal (<= %.1f A of para-F):\n" % FPROX_CUTOFF)
        fh.write(fprox_str + "\n")
    print(f"[write] {os.path.join(OUTDIR, 'pocket_residues.txt')}")


if __name__ == "__main__":
    main()
