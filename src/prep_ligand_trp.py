#!/usr/bin/env python3
"""src/prep_ligand_trp.py — build the L-tryptophan (Trp) Option-A anchor.

WHY this exists
---------------
Trp is the decisive canonical comparator for 2-Nal, not PHE. It is the only
fused-bicyclic aromatic among the 20, so it is the one canonical residue that can
actually fill a pocket shaped for naphthalene. The shape measurement in
src/analyze_shape_gap.py shows the PHE-shaped void inside the 2-Nal envelope is about
1 A^3, meaning a 2-Nal pocket does not sterically exclude PHE either. Against Trp the
discrimination lever is not volume at all: indole carries an N1-H donor (atom NE1
here) where naphthalene has a bare carbon, and burying an unsatisfied donor against
apolar surface costs 1-3 kcal/mol.

Writes data/raw/ligand_TRP.pdb: 11 HETATM (10 C + 1 N), resname TRP, chain B, resseq 1.
Atoms and ideal coordinates from RCSB CCD entry TRP (TRYPTOPHAN, C11 H12 N2 O2),
pdbx_model_Cartn_{x,y,z}_ideal columns, fetched and parsed 2026-11-19.
Source file TRP.cif sha256 4cea7e281ee404f3e72fc20d3ba98064373cdd9309a7f5c75479a04c4b6e36d6

Option A = side chain + CA cap only; 0 free-terminus oxygens; no H. Same convention as
src/prep_ligand.py (L4F, PHE) and src/prep_ligand_2nal.py (2-Nal), so all four anchors
are directly comparable.

Run:
    python src/prep_ligand_trp.py
"""
import hashlib
import os
import sys
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "data", "raw", "ligand_TRP.pdb")

# CCD TRP ideal coordinates, heavy atoms only. N / C / O / OXT are dropped by Option A:
# in the real peptide Trp is a mid-chain residue with no free terminus, and a buried
# charged carboxylate would dominate the energetics over the shape and polarity
# differences being tested.
ATOMS = [
    ("CA",  "C",  -0.008,  0.417,  1.970),
    ("CB",  "C",   0.168, -0.868,  1.161),
    ("CG",  "C",   0.650, -0.526, -0.225),
    ("CD1", "C",   1.928, -0.418, -0.622),
    ("CD2", "C",  -0.186, -0.256, -1.396),
    ("NE1", "N",   1.978, -0.095, -1.951),
    ("CE2", "C",   0.701,  0.014, -2.454),
    ("CE3", "C",  -1.564, -0.210, -1.615),
    ("CZ2", "C",   0.190,  0.314, -3.712),
    ("CZ3", "C",  -2.044,  0.086, -2.859),
    ("CH2", "C",  -1.173,  0.348, -3.907),
]

# Heavy-atom connectivity by atom name, from _chem_comp_bond in the same CCD entry.
# Indole: 5-ring CG-CD1-NE1-CE2-CD2, 6-ring CD2-CE2-CZ2-CH2-CZ3-CE3, fused at CD2/CE2.
BONDS = [
    ("CA", "CB"),
    ("CB", "CG"),
    ("CG", "CD1"), ("CG", "CD2"),
    ("CD1", "NE1"),
    ("CD2", "CE2"), ("CD2", "CE3"),
    ("NE1", "CE2"),
    ("CE2", "CZ2"),
    ("CE3", "CZ3"),
    ("CZ2", "CH2"),
    ("CZ3", "CH2"),
]

CHAIN = "B"
RESSEQ = 1
RESNAME = "TRP"

# Raw SHA-1 (hashlib) of the reference file, pinned after first verified build.
REF_SHA1 = "4b7ba79f447896636da9bd6ade07cd4b0cdd4db5"


def _hetatm(serial, name, element, x, y, z):
    # PDB fixed columns. A 1-3 char atom name gets one leading pad space.
    name4 = name.ljust(4) if len(name) == 4 else (" " + name).ljust(4)
    return (f"HETATM{serial:>5} {name4}{'':1}{RESNAME:>3} {CHAIN:1}{RESSEQ:>4}{'':1}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{1.00:6.2f}          {element:>2}")


def build():
    serial = {name: i + 1 for i, (name, *_rest) in enumerate(ATOMS)}
    lines = [
        "REMARK   L-tryptophan (Trp) side-chain anchor -- Option A",
        "REMARK   atoms+ideal coords from RCSB CCD entry TRP; heavy atoms only",
        "REMARK   side chain + CA cap; 0 free-terminus O; no H.",
        "REMARK   decisive canonical comparator for 2-Nal (only fused-bicyclic canonical).",
    ]
    for name, element, x, y, z in ATOMS:
        lines.append(_hetatm(serial[name], name, element, x, y, z))

    neigh = defaultdict(list)
    for a, b in BONDS:
        neigh[a].append(b)
        neigh[b].append(a)
    for name, *_rest in ATOMS:
        if neigh[name]:
            refs = "".join(f"{serial[m]:>5}" for m in neigh[name])
            lines.append(f"CONECT{serial[name]:>5}{refs}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def verify(content):
    het = [l for l in content.splitlines() if l.startswith("HETATM")]
    assert len(het) == 11, f"expected 11 atoms, got {len(het)}"

    elements = {l[12:16].strip(): l[76:78].strip() for l in het}
    n_c = sum(1 for e in elements.values() if e == "C")
    n_n = sum(1 for e in elements.values() if e == "N")
    n_o = sum(1 for e in elements.values() if e == "O")
    assert (n_c, n_n, n_o) == (10, 1, 0), f"expected 10 C / 1 N / 0 O, got {n_c}/{n_n}/{n_o}"
    assert elements.get("NE1") == "N", "NE1 must be the indole nitrogen"
    assert "CA" in elements and "CB" in elements, "missing CA or CB"

    names = [n for n, *_r in ATOMS]
    bonds = {(min(a, b), max(a, b)) for a, b in BONDS}
    assert all(a in names and b in names for a, b in bonds), "bond to a dropped atom"

    # Ring system = everything except the CA cap and the CB methylene.
    ring_bonds = {(a, b) for a, b in bonds if a not in ("CA", "CB") and b not in ("CA", "CB")}
    assert len(ring_bonds) == 10, f"expected 10 ring bonds, got {len(ring_bonds)}"

    deg = defaultdict(int)
    for a, b in ring_bonds:
        deg[a] += 1
        deg[b] += 1
    fusions = sorted(k for k, v in deg.items() if v == 3)
    assert fusions == ["CD2", "CE2"], f"expected fusion at CD2/CE2, got {fusions}"

    # The N1-H donor position is what a 2-Nal pocket has to penalise.
    ne1_partners = sorted({b for a, b in bonds if a == "NE1"} |
                          {a for a, b in bonds if b == "NE1"})
    assert ne1_partners == ["CD1", "CE2"], f"NE1 should bond CD1 and CE2, got {ne1_partners}"
    assert deg["NE1"] == 2, f"NE1 should have 2 ring bonds, got {deg['NE1']}"

    cb_partners = sorted({b for a, b in bonds if a == "CB"} |
                         {a for a, b in bonds if b == "CB"})
    assert cb_partners == ["CA", "CG"], f"CB should bond CA and CG, got {cb_partners}"
    ca_partners = sorted({b for a, b in bonds if a == "CA"} |
                         {a for a, b in bonds if b == "CA"})
    assert ca_partners == ["CB"], f"CA should bond only CB, got {ca_partners}"

    print("  geometry: 10 C + 1 N, 0 O, 10 ring bonds, fusion={CD2,CE2}, "
          "NE1 donor present, CB->{CA,CG}, CA->{CB}  PASS")


def main():
    content = build()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="\n") as fh:
        fh.write(content)
    sha1 = hashlib.sha1(open(OUT, "rb").read()).hexdigest()
    print(f"wrote  {OUT}")
    print(f"sha1   {sha1}")
    verify(content)
    if REF_SHA1 == "REPLACE_AFTER_FIRST_BUILD":
        print(f"NOTE: REF_SHA1 is unpinned. Set it to {sha1} in this file.")
    elif sha1 != REF_SHA1:
        print(f"ERROR: sha1 mismatch, expected {REF_SHA1}", file=sys.stderr)
        sys.exit(1)
    else:
        print("sha1 matches reference  PASS")


if __name__ == "__main__":
    main()
