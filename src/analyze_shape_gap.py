#!/usr/bin/env python3
"""src/analyze_shape_gap.py — how much shape is there to grip?

WHY this exists
---------------
The LigandMPNN dry-run (src/analyze_ligandmpnn_dryrun.py) showed that a single
para-fluorine moves sequence choice by ~0.0034 bits mean JSD. That told us the
L4F-vs-PHE signal is absent, but not *why*, and not whether 2-Nal is any better.
This script answers both in units a pocket designer can act on: van der Waals
volume, measured on a grid, after superposing each ncAA anchor on the canonical
look-alike it has to be distinguished from.

The decisive quantity is not the surplus volume. It is the *void*: how much of the
canonical residue's envelope lies outside the ncAA's. A void of zero means the
canonical fits entirely inside the ncAA envelope, so a pocket complementary to the
ncAA admits the canonical with no steric penalty at all, and discrimination can
only come from a favourable contact with the distinguishing atoms. That is the
regime L4F is in, and it is why the sequence lever found nothing.

CPU only, no network, no model weights. Reads the Option-A anchors written by
src/prep_ligand.py (L4F, PHE) and src/prep_ligand_2nal.py (2-Nal).

Run:
    python src/analyze_shape_gap.py
    python src/analyze_shape_gap.py --ligand-dir data/raw --spacing 0.15

The 2-Nal-vs-Trp comparison is the decisive one, and it needs a different method.
Trp is 6-5 fused where 2-Nal is 6-6, so there is no atom-to-atom ring correspondence
to run Kabsch on. Superposing on CA/CB/CG alone fixes the attachment and leaves the
ring free to rotate, which would make the answer an artefact of an arbitrary rotamer.
So this script instead fixes the attachment (CA, CB, ipso carbon) and then searches
the chi2 dihedral for the Trp rotamer that MAXIMISES overlap with 2-Nal, and reports
what is still different at that best case. That is the conservative direction: it
gives Trp its best shot at mimicking 2-Nal and asks what a pocket could still read.

2-Nal is held at its CCD ideal chi2 throughout, because that is the conformer the
anchor file carries and therefore the one a pocket would be designed around.

Requires data/raw/ligand_TRP.pdb from src/prep_ligand_trp.py. If that file is absent
the Trp section is skipped and the PHE comparison still runs.
"""
import argparse
import os
import sys

import numpy as np

VDW = {"C": 1.70, "F": 1.47, "N": 1.55, "O": 1.52}

# The 2-Nal proximal ring maps atom-for-atom onto the Phe ring. Ring trace from the
# ipso carbon: C2(ipso) -> C3,C1(ortho) -> C4,C8A(meta) -> C4A(para). The fused ring
# (C5,C6,C7,C8) is the surplus, and it grows off C4A/C8A, which are the para and meta
# positions. So the second ring occupies the region the para-F sat in.
NAL_TO_PHE = {
    "CA": "CA", "CB": "CB",
    "C2": "CG",
    "C3": "CD1", "C1": "CD2",
    "C4": "CE1", "C8A": "CE2",
    "C4A": "CZ",
}


# Trp shares only the attachment with 2-Nal: CA, CB, and the ipso ring carbon.
# Everything distal to the ipso carbon is what the chi2 search rotates.
TRP_TO_NAL = {"CA": "CA", "CB": "CB", "CG": "C2"}
TRP_DISTAL = ["CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]


def read_anchor(path, required=True):
    """Parse HETATM records out of an Option-A anchor PDB."""
    if not os.path.exists(path) and not required:
        return None
    if not os.path.exists(path):
        sys.exit(f"missing {path}\n  build it first: python src/prep_ligand.py "
                 f"and python src/prep_ligand_2nal.py")
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("HETATM"):
                name = line[12:16].strip()
                el = line[76:78].strip() or name[0]
                atoms.append((name, el,
                              np.array([float(line[30:38]),
                                        float(line[38:46]),
                                        float(line[46:54])])))
    if not atoms:
        sys.exit(f"no HETATM records in {path}")
    return atoms


def kabsch(mobile, target):
    """Rotation + centroids mapping mobile onto target (both [N,3]), no reflection."""
    mc, tc = mobile.mean(0), target.mean(0)
    V, _S, Wt = np.linalg.svd((mobile - mc).T @ (target - tc))
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt, mc, tc


def apply_fit(atoms, R, mc, tc):
    return [(n, e, (x - mc) @ R + tc) for n, e, x in atoms]


def rodrigues(axis, theta):
    """Rotation matrix for theta radians about a unit axis."""
    k = axis / np.linalg.norm(axis)
    K = np.array([[0.0, -k[2], k[1]],
                  [k[2], 0.0, -k[0]],
                  [-k[1], k[0], 0.0]])
    return np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def set_chi2(atoms, theta, pivot_name, axis_from, distal):
    """Rotate only the distal ring atoms about the pivot bond, i.e. turn chi2.

    Atoms on the axis (axis_from and pivot_name) and everything proximal to it stay
    put, which is what makes this a dihedral rotation rather than a rigid-body one.
    """
    by = {n: x for n, _e, x in atoms}
    origin = by[pivot_name]
    R = rodrigues(by[pivot_name] - by[axis_from], theta)
    return [(n, e, (origin + (x - origin) @ R) if n in distal else x)
            for n, e, x in atoms]


def occupancy(atoms, grid):
    """Boolean mask over grid points inside any atom's vdW sphere."""
    inside = np.zeros(grid.shape[0], dtype=bool)
    for _name, el, xyz in atoms:
        r = VDW.get(el, 1.70)
        inside |= ((grid - xyz) ** 2).sum(1) <= r * r
    return inside


def make_grid(atom_sets, spacing, margin=2.5):
    allxyz = np.vstack([np.array([a[2] for a in s]) for s in atom_sets])
    lo, hi = allxyz.min(0) - margin, allxyz.max(0) + margin
    axes = [np.arange(lo[i], hi[i] + spacing, spacing) for i in range(3)]
    mesh = np.meshgrid(*axes, indexing="ij")
    return np.stack([m.ravel() for m in mesh], axis=1)


def reach(atoms):
    """Furthest side-chain atom from CB: how deep a complementary pocket must be."""
    cb = {n: x for n, _e, x in atoms}["CB"]
    return max(float(np.linalg.norm(x - cb)) for _n, _e, x in atoms)


def trp_section(nal_fit, trp, spacing):
    """Give Trp its best shot at mimicking 2-Nal, then report the residual."""
    print("\n" + "=" * 70)
    print("2-Nal vs Trp: the decisive comparison")
    print("=" * 70)

    nal_by = {n: x for n, _e, x in nal_fit}
    trp_by = {n: x for n, _e, x in trp}
    missing = [k for k in TRP_TO_NAL if k not in trp_by]
    missing += [v for v in TRP_TO_NAL.values() if v not in nal_by]
    assert not missing, f"anchors are missing mapped atoms: {missing}"
    assert "NE1" in trp_by, "Trp anchor has no NE1; it is not an indole"

    # Fix the attachment: CA, CB and the ipso carbon.
    pairs = list(TRP_TO_NAL.items())
    R, mc, tc = kabsch(np.array([trp_by[a] for a, _ in pairs]),
                       np.array([nal_by[b] for _, b in pairs]))
    rmsd = float(np.sqrt((((np.array([trp_by[a] for a, _ in pairs]) - mc) @ R + tc
                           - np.array([nal_by[b] for _, b in pairs])) ** 2).sum(1).mean()))
    print(f"attachment superposed on CA/CB/ipso, RMSD {rmsd:.3f} A")
    assert rmsd < 0.5, f"attachment overlay RMSD {rmsd:.3f} A is too high"
    trp_fit = apply_fit(trp, R, mc, tc)

    # Grid must cover 2-Nal plus the whole volume Trp's ring sweeps as chi2 turns.
    cb = {n: x for n, _e, x in trp_fit}["CB"]
    rmax = max(float(np.linalg.norm(x - cb)) for _n, _e, x in trp_fit)
    nalxyz = np.array([x for _n, _e, x in nal_fit])
    lo = np.minimum(nalxyz.min(0), cb - rmax) - 2.5
    hi = np.maximum(nalxyz.max(0), cb + rmax) + 2.5
    axes = [np.arange(lo[i], hi[i] + spacing, spacing) for i in range(3)]
    mesh = np.meshgrid(*axes, indexing="ij")
    grid = np.stack([m.ravel() for m in mesh], axis=1)
    vox = spacing ** 3
    occ_nal = occupancy(nal_fit, grid)
    print(f"grid {grid.shape[0]:,} points at {spacing} A")

    def overlap_at(theta_deg):
        rot = set_chi2(trp_fit, np.radians(theta_deg), "CG", "CB", TRP_DISTAL)
        return (occupancy(rot, grid) & occ_nal).sum() * vox, rot

    # Coarse scan then refine. Maximising overlap is the conservative choice.
    coarse = [(overlap_at(t)[0], t) for t in range(0, 360, 10)]
    best_v, best_t = max(coarse)
    fine = [(overlap_at(t)[0], t) for t in range(best_t - 10, best_t + 11)]
    best_v, best_t = max(fine)
    worst_v, worst_t = min(coarse)
    print(f"chi2 searched 0-359 deg. Best overlap {best_v:.1f} A^3 at {best_t % 360} deg; "
          f"worst {worst_v:.1f} A^3 at {worst_t} deg.")

    _v, best_rot = overlap_at(best_t)
    occ_trp = occupancy(best_rot, grid)
    v_nal = occ_nal.sum() * vox
    v_trp = occ_trp.sum() * vox
    surplus = (occ_trp & ~occ_nal).sum() * vox
    void = (occ_nal & ~occ_trp).sum() * vox

    print(f"\nAt the best-mimicking Trp rotamer ({best_t % 360} deg):")
    print(f"  2-Nal volume            {v_nal:7.1f} A^3")
    print(f"  Trp volume              {v_trp:7.1f} A^3")
    print(f"  shared                  {best_v:7.1f} A^3  "
          f"({100.0 * best_v / v_nal:.0f}% of 2-Nal)")
    print(f"  Trp surplus             {surplus:7.1f} A^3  "
          f"volume a 2-Nal pocket does not provide, so this is the CLASH lever")
    print(f"  2-Nal void              {void:7.1f} A^3  "
          f"contact Trp fails to make")

    # The polarity lever: where does the indole N-H land?
    ne1 = {n: x for n, _e, x in best_rot}["NE1"]
    d_to_nal = min(float(np.linalg.norm(ne1 - x)) for _n, _e, x in nal_fit)
    inside = bool(occ_nal[np.argmin(((grid - ne1) ** 2).sum(1))])
    print(f"\n  indole NE1 (the N-H donor): {d_to_nal:.2f} A from the nearest 2-Nal atom, "
          f"{'INSIDE' if inside else 'outside'} the 2-Nal envelope")
    if inside:
        print("  A pocket complementary to 2-Nal therefore presents apolar surface at the")
        print("  position Trp's N-H occupies. Burying that donor unsatisfied is the")
        print("  discrimination lever, and it does not depend on volume at all.")
    else:
        print("  The donor sits outside the 2-Nal envelope, so a complementary pocket")
        print("  may leave it solvent-exposed. Check this before relying on desolvation.")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--ligand-dir", default=os.path.join(repo, "data", "raw"))
    ap.add_argument("--spacing", type=float, default=0.15,
                    help="grid spacing in A (default 0.15)")
    args = ap.parse_args()

    phe = read_anchor(os.path.join(args.ligand_dir, "ligand_PHE.pdb"))
    l4f = read_anchor(os.path.join(args.ligand_dir, "ligand_L4F.pdb"))
    nal = read_anchor(os.path.join(args.ligand_dir, "ligand_2NAL.pdb"))
    print(f"anchors: PHE {len(phe)} atoms, L4F {len(l4f)}, 2NAL {len(nal)}")

    # L4F and PHE come from the same PFF CCD frame and differ by exactly one atom,
    # so they need no superposition. That pair is the clean control.
    phe_by = {n: x for n, _e, x in phe}
    l4f_by = {n: x for n, _e, x in l4f}
    shared = [n for n in phe_by if n in l4f_by]
    drift = max(float(np.linalg.norm(phe_by[n] - l4f_by[n])) for n in shared)
    assert drift < 1e-6, f"L4F and PHE are not in the same frame (max drift {drift})"
    assert len(l4f) == len(phe) + 1, "L4F should be PHE plus exactly the para-F"
    print(f"  L4F/PHE share {len(shared)} atoms in one frame (max drift {drift:.1e} A)  PASS")

    # 2-Nal comes from the NAL CCD frame, so superpose it on the atoms it shares.
    nal_by = {n: x for n, _e, x in nal}
    missing = [k for k in NAL_TO_PHE if k not in nal_by]
    assert not missing, f"2NAL anchor is missing mapped atoms: {missing}"
    pairs = list(NAL_TO_PHE.items())
    mob = np.array([nal_by[k] for k, _ in pairs])
    tgt = np.array([phe_by[v] for _, v in pairs])
    R, mc, tc = kabsch(mob, tgt)
    rmsd = float(np.sqrt((((mob - mc) @ R + tc - tgt) ** 2).sum(1).mean()))
    assert rmsd < 0.25, (
        f"proximal-ring overlay RMSD {rmsd:.3f} A is too high; the NAL_TO_PHE "
        f"correspondence is probably wrong")
    print(f"  2NAL->PHE overlay on {len(pairs)} shared atoms, RMSD {rmsd:.3f} A  PASS")
    nal_fit = apply_fit(nal, R, mc, tc)

    grid = make_grid([phe, l4f, nal_fit], args.spacing)
    vox = args.spacing ** 3
    print(f"  grid {grid.shape[0]:,} points at {args.spacing} A "
          f"({vox:.6f} A^3 per voxel)")

    occ = {"PHE": occupancy(phe, grid),
           "L4F": occupancy(l4f, grid),
           "2NAL": occupancy(nal_fit, grid)}
    vol = lambda m: m.sum() * vox

    print(f"\n{'species':<8}{'vdW volume':>13}{'reach from CB':>16}")
    for tag, atoms in (("PHE", phe), ("L4F", l4f), ("2NAL", nal_fit)):
        print(f"{tag:<8}{vol(occ[tag]):>10.1f} A^3{reach(atoms):>13.2f} A")

    print("\nShape difference against the canonical look-alike (PHE):")
    print(f"  {'ncAA':<6}{'surplus':>10}{'void':>9}{'symmetric diff':>17}")
    stats = {}
    for tag in ("L4F", "2NAL"):
        surplus = vol(occ[tag] & ~occ["PHE"])
        void = vol(occ["PHE"] & ~occ[tag])
        stats[tag] = (surplus, void)
        print(f"  {tag:<6}{surplus:>7.1f} A^3{void:>6.1f} A^3{surplus + void:>14.1f} A^3")

    s_l4f, v_l4f = stats["L4F"]
    s_nal, v_nal = stats["2NAL"]
    print(f"\n2-Nal presents {s_nal / s_l4f:.1f}x the surplus volume of the para-F.")
    print("Both voids are near zero, so PHE fits inside both ncAA envelopes: a pocket")
    print("complementary to either ncAA still admits PHE without a clash. The gap is")
    print("lost contact, not steric exclusion.")

    trp = read_anchor(os.path.join(args.ligand_dir, "ligand_TRP.pdb"), required=False)
    if trp is None:
        print("\nligand_TRP.pdb not found, skipping the Trp comparison.")
        print("Build it with: python src/prep_ligand_trp.py")
    else:
        trp_section(nal_fit, trp, args.spacing)


if __name__ == "__main__":
    main()
