#!/usr/bin/env python3
"""
src/prep_ligand_2nal.py — build the L-2-naphthylalanine (2-Nal) Option-A anchor.

Writes data/raw/ligand_2NAL.pdb: 12 HETATM (all C), resname NAL, chain B, resseq 1.
Atoms and ideal coordinates from RCSB CCD entry NAL (verified 2026-06-07).
CCD atom C9 is renamed CB (side-chain beta carbon); CA is the alpha-carbon cap.

Option A = side chain + CA cap only; 0 free-terminus oxygens; no H.
2-Nal is the PRIMARY ncAA target (pivot from L4F, 2026-06-07).
L4F and PHE controls remain in src/prep_ligand.py — that file is NOT modified here.
This script is 2-Nal specific; the cloud-input builder is the one that takes --residue.
"""
import hashlib
import os
import sys
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(REPO, "data", "raw", "ligand_2NAL.pdb")

# CCD NAL ideal coordinates, heavy atoms only.
# Source: pdbx_model_Cartn_{x,y,z}_ideal columns, NAL.cif (RCSB, fetched 2026-06-07).
# CCD atom C9 (the beta methylene) is renamed CB per standard PDB protein convention.
ATOMS = [
    ("CA",   2.446,  0.058,  0.429),
    ("CB",   1.670, -0.304, -0.839),
    ("C1",  -0.616,  0.617, -0.580),
    ("C2",   0.209, -0.463, -0.507),
    ("C3",  -0.277, -1.713, -0.129),
    ("C4",  -1.592, -1.890,  0.178),
    ("C4A", -2.474, -0.799,  0.113),
    ("C5",  -3.837, -0.943,  0.422),
    ("C6",  -4.662,  0.137,  0.348),
    ("C7",  -4.176,  1.387, -0.030),
    ("C8",  -2.862,  1.564, -0.337),
    ("C8A", -1.978,  0.474, -0.266),
]

# CONECT records (1-based serial indices matching ATOMS order above).
CONECTS = [
    (1,  [2]),
    (2,  [4, 1]),
    (3,  [12, 4]),
    (4,  [5, 3, 2]),
    (5,  [6, 4]),
    (6,  [7, 5]),
    (7,  [8, 6, 12]),
    (8,  [9, 7]),
    (9,  [8, 10]),
    (10, [9, 11]),
    (11, [10, 12]),
    (12, [7, 11, 3]),
]

# Raw SHA-1 (hashlib) of the reference file verified and approved 2026-06-07.
# (git hash-object of the same file is f1f6ede6c33d3edb0b2aa1116d6989aa8ef93ee9)
REF_SHA1 = "ac9f19f7b913fe5130be11ecd0f3a6d2bafa9bc3"


def _hetatm(serial, name, x, y, z):
    return (f"HETATM{serial:5d} " + f" {name:<3}" +
            f" NAL B   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  1.00           C")


def _conect(serial, partners):
    return "CONECT" + "".join(f"{s:5d}" for s in [serial] + list(partners))


def build():
    lines = [
        "REMARK   L-2-naphthylalanine (2-Nal) side-chain anchor -- Option A",
        "REMARK   atoms+ideal coords from RCSB CCD entry NAL; heavy atoms only",
        "REMARK   side chain + CA cap; 0 free-terminus O; no H. PRIMARY ncAA target.",
    ]
    for i, (name, x, y, z) in enumerate(ATOMS, 1):
        lines.append(_hetatm(i, name, x, y, z))
    for serial, partners in CONECTS:
        lines.append(_conect(serial, partners))
    lines.append("END")
    return "\n".join(lines) + "\n"


def verify(content):
    hetlines = [l for l in content.splitlines() if l.startswith("HETATM")]
    n = len(hetlines)
    assert n == 12, f"expected 12 atoms, got {n}"

    elements = {l[12:16].strip(): l[76:78].strip() for l in hetlines}
    assert all(e == "C" for e in elements.values()), f"non-C elements: {elements}"
    assert "CA" in elements and "CB" in elements, "missing CA or CB"

    bond_set = set()
    for l in content.splitlines():
        if not l.startswith("CONECT"):
            continue
        nfields = (len(l) - 6) // 5
        cols = [int(l[6 + 5*i: 6 + 5*(i+1)])
                for i in range(nfields)
                if l[6 + 5*i: 6 + 5*(i+1)].strip()]
        for p in cols[1:]:
            bond_set.add((min(cols[0], p), max(cols[0], p)))

    ring_cc = {(a, b) for a, b in bond_set if a >= 3 and b >= 3}
    assert len(ring_cc) == 11, f"expected 11 ring C-C bonds, got {len(ring_cc)}"

    deg = defaultdict(int)
    for a, b in ring_cc:
        deg[a] += 1
        deg[b] += 1
    fusions = sorted(k for k, v in deg.items() if v == 3)
    assert fusions == [7, 12], (
        f"expected fusion carbons at serial indices 7 (C4A) and 12 (C8A), got {fusions}")
    assert (7, 12) in ring_cc, "C4A(7)–C8A(12) bond missing"

    cb_bonds = ({b for a, b in bond_set if a == 2} |
                {a for a, b in bond_set if b == 2})
    assert cb_bonds == {1, 4}, (
        f"CB should bond only CA(1) and C2(4), got {cb_bonds}")

    ca_bonds = ({b for a, b in bond_set if a == 1} |
                {a for a, b in bond_set if b == 1})
    assert ca_bonds == {2}, (
        f"CA should bond only CB(2), got {ca_bonds}")

    o_count = sum(1 for l in hetlines if l[76:78].strip() == "O")
    assert o_count == 0, f"O count must be 0, got {o_count}"

    print("  geometry: 12 C, 0 O, 11 ring bonds, "
          "fusion={C4A,C8A}, CB→{CA,C2}, CA→{CB}  PASS")


def main():
    content = build()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="\n") as fh:
        fh.write(content)
    sha1 = hashlib.sha1(open(OUT, "rb").read()).hexdigest()
    print(f"wrote  {OUT}")
    print(f"sha1   {sha1}")
    if REF_SHA1 and sha1 != REF_SHA1:
        print(f"ERROR: sha1 mismatch — expected {REF_SHA1}", file=sys.stderr)
        sys.exit(1)
    elif REF_SHA1:
        print(f"sha1 matches reference  PASS")
    verify(content)


if __name__ == "__main__":
    main()
