#!/usr/bin/env python3
"""src/check_n1_pocket.py — does this pocket actually discriminate 2-Nal from Trp?

WHY this exists
---------------
src/analyze_shape_gap.py established that volume is not the lever. Against the
best-mimicking Trp rotamer, Trp fills 76% of the 2-Nal envelope, so a pocket selected
for being large admits Trp and discriminates nothing. What Trp cannot change is that
its indole N-H (atom NE1) lands 0.56 A from a naphthalene carbon position, inside the
2-Nal envelope. A pocket complementary to 2-Nal therefore has apolar surface there, and
burying that donor unsatisfied costs 1 to 3 kcal/mol.

RFdiffusionAA has no objective term for this. It conditions on ligand atoms and nothing
else, so the criterion has to be applied as a filter after designs exist. That is what
this script is: a per-design PASS/FAIL on the only lever that survived measurement.

Two conditions must both hold for the lever to be real:
  1. NE1 is BURIED. If the donor points at solvent there is no desolvation penalty to
     collect, however apolar the surroundings are.
  2. NO HYDROGEN-BOND ACCEPTOR sits within reach of NE1. A backbone carbonyl, or a Ser,
     Thr, Tyr, Asn, Gln, Asp or Glu oxygen, satisfies the N-H and throws the gap away.
     This is the failure mode worth filtering for, because a sequence-design step will
     happily place a polar residue there.

The Trp rotamer is chosen to MAXIMISE overlap with the 2-Nal ligand, the same definition
analyze_shape_gap.py uses, so the numbers are comparable and the choice is conservative.

Run:
    python src/check_n1_pocket.py --pocket cloud/rfdiffusion_aa/input_2NAL_pocket.pdb
    python src/check_n1_pocket.py --pocket output/2nal_pocket/des_0.pdb

Exit status is 0 on PASS and 1 on FAIL, so it can gate a loop over many designs.

Run it on the cloud input and it FAILS, reporting only 3 protein atoms within 5 A of
NE1. That is correct and worth keeping as the filter's negative control: the 4QAE
carrier was never designed around 2-Nal, the anchor was rigid-placed into it with
clashes, and the carrier is discarded by the de-novo contig anyway. A real test needs a
design that came out of RFdiffusionAA.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from analyze_shape_gap import (read_anchor, kabsch, apply_fit, set_chi2, occupancy,
                               TRP_TO_NAL, TRP_DISTAL)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRP_ANCHOR = os.path.join(REPO, "data", "raw", "ligand_TRP.pdb")

HBOND_CUTOFF = 3.5      # A, N-H...acceptor heavy-atom distance
BURIAL_CUTOFF = 5.0     # A, protein heavy atoms counted around NE1
BURIAL_MIN = 8          # fewer neighbours than this and NE1 is not buried

# Side-chain oxygens and His ring nitrogens can accept from an N-H. Backbone N is a
# donor, not an acceptor, so a bare "N" is not counted.
ACCEPTOR_NAMES = {
    "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OXT",
    "ND1", "NE2",
}


def parse_pocket(path, ligand_resname):
    """Split a pocket PDB into protein atoms and the ncAA ligand."""
    protein, ligand = [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element == "H":
                continue
            rec = {
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resseq": int(line[22:26]),
                "element": element,
                "xyz": np.array([float(line[30:38]), float(line[38:46]),
                                 float(line[46:54])]),
            }
            if rec["resname"] == ligand_resname:
                ligand.append((rec["name"], rec["element"], rec["xyz"]))
            elif line.startswith("ATOM"):
                protein.append(rec)
    return protein, ligand


def best_trp_rotamer(ligand, trp, spacing=0.2):
    """Superpose Trp on the ligand attachment, then turn chi2 for maximum overlap."""
    lig_by = {n: x for n, _e, x in ligand}
    trp_by = {n: x for n, _e, x in trp}
    missing = [a for a in TRP_TO_NAL if a not in trp_by]
    missing += [b for b in TRP_TO_NAL.values() if b not in lig_by]
    if missing:
        sys.exit(f"cannot place Trp, missing mapped atoms: {missing}\n"
                 f"  the ligand must be the 2-Nal anchor (CA, CB, C2 present)")

    pairs = list(TRP_TO_NAL.items())
    mob = np.array([trp_by[a] for a, _ in pairs])
    tgt = np.array([lig_by[b] for _, b in pairs])
    R, mc, tc = kabsch(mob, tgt)
    rmsd = float(np.sqrt((((mob - mc) @ R + tc - tgt) ** 2).sum(1).mean()))
    trp_fit = apply_fit(trp, R, mc, tc)

    cb = {n: x for n, _e, x in trp_fit}["CB"]
    rmax = max(float(np.linalg.norm(x - cb)) for _n, _e, x in trp_fit)
    ligxyz = np.array([x for _n, _e, x in ligand])
    lo = np.minimum(ligxyz.min(0), cb - rmax) - 2.0
    hi = np.maximum(ligxyz.max(0), cb + rmax) + 2.0
    axes = [np.arange(lo[i], hi[i] + spacing, spacing) for i in range(3)]
    mesh = np.meshgrid(*axes, indexing="ij")
    grid = np.stack([m.ravel() for m in mesh], axis=1)
    occ_lig = occupancy(ligand, grid)

    def overlap(theta):
        rot = set_chi2(trp_fit, np.radians(theta), "CG", "CB", TRP_DISTAL)
        return int((occupancy(rot, grid) & occ_lig).sum()), rot

    coarse = max((overlap(t)[0], t) for t in range(0, 360, 10))
    best = max((overlap(t)[0], t) for t in range(coarse[1] - 10, coarse[1] + 11))
    vox = spacing ** 3
    return best[1] % 360, overlap(best[1])[1], best[0] * vox, rmsd


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pocket", required=True,
                    help="PDB with protein ATOM records and the 2-Nal HETATM ligand")
    ap.add_argument("--ligand-resname", default="NAL")
    ap.add_argument("--trp-anchor", default=TRP_ANCHOR)
    args = ap.parse_args()

    protein, ligand = parse_pocket(args.pocket, args.ligand_resname)
    if not ligand:
        sys.exit(f"no HETATM with resname {args.ligand_resname} in {args.pocket}")
    if not protein:
        sys.exit(f"no protein ATOM records in {args.pocket}")
    print(f"pocket : {args.pocket}")
    print(f"         {len(protein)} protein heavy atoms, "
          f"{len(ligand)} ligand atoms (resname {args.ligand_resname})")

    trp = read_anchor(args.trp_anchor)
    chi2, trp_placed, shared, rmsd = best_trp_rotamer(ligand, trp)
    print(f"Trp    : attachment RMSD {rmsd:.3f} A, best-mimicking chi2 {chi2} deg, "
          f"shared volume {shared:.1f} A^3")

    ne1 = {n: x for n, _e, x in trp_placed}["NE1"]
    pxyz = np.array([a["xyz"] for a in protein])
    d = np.linalg.norm(pxyz - ne1, axis=1)

    n_near = int((d <= BURIAL_CUTOFF).sum())
    buried = n_near >= BURIAL_MIN
    print(f"\nindole NE1 (the donor the pocket must leave unsatisfied):")
    print(f"  protein heavy atoms within {BURIAL_CUTOFF} A : {n_near} "
          f"({'buried' if buried else 'NOT buried'}, threshold {BURIAL_MIN})")

    acceptors = [(protein[i], float(d[i])) for i in np.argsort(d)
                 if d[i] <= HBOND_CUTOFF
                 and protein[i]["name"] in ACCEPTOR_NAMES]
    if acceptors:
        print(f"  hydrogen-bond acceptors within {HBOND_CUTOFF} A : {len(acceptors)}")
        for a, dist in acceptors:
            print(f"    {a['resname']}{a['resseq']}{a['chain']} {a['name']:<4} "
                  f"{dist:.2f} A")
    else:
        print(f"  hydrogen-bond acceptors within {HBOND_CUTOFF} A : none")

    lining = {}
    for i in np.where(d <= BURIAL_CUTOFF)[0]:
        a = protein[i]
        key = (a["chain"], a["resseq"], a["resname"])
        lining[key] = min(lining.get(key, 99.0), float(d[i]))
    if lining:
        print(f"  residues lining the N1 position (within {BURIAL_CUTOFF} A):")
        for (ch, rs, rn), dist in sorted(lining.items(), key=lambda kv: kv[1]):
            print(f"    {rn}{rs}{ch} at {dist:.2f} A")

    ok = buried and not acceptors
    print()
    if ok:
        print("VERDICT: PASS. NE1 is buried with no acceptor in reach, so this pocket")
        print("         would bury Trp's N-H unsatisfied. The discrimination lever is")
        print("         intact. Shape alone does not make this pocket good; verify the")
        print("         affinity gap against data/raw/ligand_TRP.pdb before believing it.")
    else:
        print("VERDICT: FAIL.")
        if not buried:
            print("         NE1 is not buried, so there is no desolvation penalty to")
            print("         collect and Trp pays nothing for binding here.")
        if acceptors:
            print("         An acceptor is within reach of NE1, which satisfies Trp's")
            print("         N-H and discards the only lever that survived measurement.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
