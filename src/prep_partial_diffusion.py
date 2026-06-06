"""Prepare the Anticalin+Key complex for an RFdiffusion partial-diffusion run.

Does two things, deterministically and without touching the network:

  1. Cleans ``data/raw/scaffold_complex.pdb`` -> ``data/raw/scaffold_clean.pdb``:
     keeps only ATOM records for the scaffold and Key chains, dropping
     crystallographic waters, ANISOU records, hetero-atoms and alt-loc copies.
     Original residue numbering is preserved (RFdiffusion contigs reference it).

  2. Emits a length-checked RFdiffusion contig to ``data/raw/contig.txt``.
     Default ("loop melt") freezes the beta-strand framework and the Key as
     fixed motifs and diffuses ONLY the scaffold residues that contact the Key
     (the binding loops). ``--whole`` instead melts the whole scaffold gently.

Partial diffusion requires the contig to be exactly the same length as the
input, so this script asserts that the emitted contig tiles every residue.
No RFdiffusion is run here -- see src/run_partial_diffusion.sh for that.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


class CleanSelect(Select):
    """Keep standard-residue ATOM records for the two chains, first alt-loc only."""

    def __init__(self, keep_chains: set[str]):
        self.keep_chains = keep_chains

    def accept_chain(self, chain):
        return chain.id in self.keep_chains

    def accept_residue(self, residue):
        # hetfield is residue.id[0]; ' ' means a standard polymer residue
        return residue.id[0] == " " and residue.resname.strip() in STD_AA

    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")


def ordered_residues(chain):
    """Polymer residues of a chain, sorted by sequence number."""
    res = [r for r in chain if r.id[0] == " " and r.resname.strip() in STD_AA]
    return sorted(res, key=lambda r: r.id[1])


def min_distance_to_key(residue, key_xyz: np.ndarray) -> float:
    coords = np.array([a.coord for a in residue.get_atoms()], dtype=float)
    # pairwise min euclidean distance, numpy-only
    d = np.sqrt(((coords[:, None, :] - key_xyz[None, :, :]) ** 2).sum(-1))
    return float(d.min())


def _runs(flags):
    """Yield (value, start_idx, end_idx_inclusive) for each maximal constant run."""
    i, n = 0, len(flags)
    while i < n:
        j = i
        while j < n and flags[j] == flags[i]:
            j += 1
        yield flags[i], i, j - 1
        i = j


def clean_melt_mask(flags, close_gap: int, min_run: int, pad: int):
    """Tidy a per-residue melt mask into coherent loop runs.

    Raw distance thresholding fragments a single loop: antiparallel strands put
    spatially-adjacent contacts at distant sequence positions, and sub-cutoff
    residues inside a loop punch frozen holes through it. Three morphological
    passes fix that:

      close_gap  fill interior frozen gaps this short or shorter (merge a loop
                 split by a couple of sub-cutoff residues),
      min_run    drop melt runs shorter than this (lone strand side-chain
                 contacts that are not loops),
      pad        grow each surviving run by this many residues per side (slack
                 for RFdiffusion to re-stitch the loop<->strand junctions).

    Residue count never changes -- only frozen<->melt labels move -- so the
    exact-length contig invariant is preserved.
    """
    flags, n = list(flags), len(flags)
    for val, a, b in list(_runs(flags)):  # closing
        if not val and 0 < a and b < n - 1 and (b - a + 1) <= close_gap:
            flags[a:b + 1] = [True] * (b - a + 1)
    for val, a, b in list(_runs(flags)):  # opening
        if val and (b - a + 1) < min_run:
            flags[a:b + 1] = [False] * (b - a + 1)
    if pad:  # dilation
        grown = list(flags)
        for val, a, b in list(_runs(flags)):
            if val:
                grown[max(0, a - pad):b + pad + 1] = [True] * (
                    min(n, b + pad + 1) - max(0, a - pad))
        flags = grown
    return flags


def build_loop_contig(scaffold_res, key_res, threshold, close_gap, min_run, pad):
    key_xyz = np.array(
        [a.coord for r in key_res for a in r.get_atoms()], dtype=float
    )
    melt_flags = [min_distance_to_key(r, key_xyz) < threshold for r in scaffold_res]
    melt_flags = clean_melt_mask(melt_flags, close_gap, min_run, pad)

    segments, melt_total, loops = [], 0, []
    chain_id = scaffold_res[0].get_parent().id
    for val, a, b in _runs(melt_flags):
        run = scaffold_res[a:b + 1]
        if val:
            segments.append(f"{len(run)}-{len(run)}")  # diffused
            melt_total += len(run)
            loops.append(f"{run[0].id[1]}-{run[-1].id[1]}")
        else:
            segments.append(f"{chain_id}{run[0].id[1]}-{run[-1].id[1]}")  # frozen
    return segments, melt_total, loops


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="data/raw/scaffold_complex.pdb")
    ap.add_argument("--clean-out", default="data/raw/scaffold_clean.pdb")
    ap.add_argument("--contig-out", default="data/raw/contig.txt")
    ap.add_argument("--scaffold-chain", default="A")
    ap.add_argument("--key-chain", default="B")
    ap.add_argument("--contact-threshold", type=float, default=6.0,
                    help="Angstrom cutoff: scaffold residues within this of the "
                         "Key are melted (loop-melt mode).")
    ap.add_argument("--close-gap", type=int, default=2,
                    help="Loop-melt: merge melt fragments separated by a frozen "
                         "gap this short or shorter (default 2).")
    ap.add_argument("--min-run", type=int, default=3,
                    help="Loop-melt: drop melt runs shorter than this -- lone "
                         "strand side-chain contacts, not loops (default 3).")
    ap.add_argument("--pad", type=int, default=1,
                    help="Loop-melt: grow each loop by this many residues per side "
                         "for strand-junction freedom (default 1).")
    ap.add_argument("--whole", action="store_true",
                    help="Melt the whole scaffold gently instead of just loops.")
    args = ap.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"Input not found: {in_path}")

    structure = PDBParser(QUIET=True).get_structure("complex", str(in_path))
    model = structure[0]
    for cid in (args.scaffold_chain, args.key_chain):
        if cid not in model:
            raise SystemExit(f"Chain '{cid}' not in {in_path}")

    scaffold_res = ordered_residues(model[args.scaffold_chain])
    key_res = ordered_residues(model[args.key_chain])
    if any(r.id[2] != " " for r in scaffold_res + key_res):
        print("WARNING: insertion codes present; verify contig ranges by hand.")

    n_scaf, n_key = len(scaffold_res), len(key_res)
    key_lo, key_hi = key_res[0].id[1], key_res[-1].id[1]
    kc = args.key_chain

    # 1. Write the cleaned PDB.
    Path(args.clean_out).parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(args.clean_out,
            CleanSelect({args.scaffold_chain, args.key_chain}))

    # 2. Build the contig.
    if args.whole:
        segments = [f"{n_scaf}-{n_scaf}"]
        melt_total = n_scaf
        loops = None
    else:
        segments, melt_total, loops = build_loop_contig(
            scaffold_res, key_res, args.contact_threshold,
            args.close_gap, args.min_run, args.pad)

    contig = f"[{'/'.join(segments)}/0 {kc}{key_lo}-{key_hi}]"

    # Length invariant: partial diffusion demands an exact-length contig.
    def seg_len(s: str) -> int:
        lo, hi = s.split("-")
        if lo[:1].isdigit():  # diffused "L-L": both numbers are the run length
            return int(hi)
        # frozen motif "A{start}-{end}": numbers are residue positions
        return int(hi) - int(lo.lstrip(args.scaffold_chain)) + 1

    tiled = sum(seg_len(s) for s in segments)
    assert tiled == n_scaf, f"contig tiles {tiled} scaffold residues, expected {n_scaf}"

    Path(args.contig_out).write_text(contig + "\n")

    mode = "whole-scaffold melt" if args.whole else "loop melt"
    print(f"Cleaned complex -> {args.clean_out}")
    print(f"  scaffold chain {args.scaffold_chain}: {n_scaf} residues "
          f"({scaffold_res[0].id[1]}-{scaffold_res[-1].id[1]})")
    print(f"  key chain {kc}: {n_key} residues ({key_lo}-{key_hi})")
    print(f"Mode: {mode}  |  melted {melt_total}/{n_scaf} scaffold residues")
    if loops:
        print(f"  loops melted ({len(loops)}): {', '.join(loops)}")
    print(f"Contig ({n_scaf + n_key} residues total) -> {args.contig_out}")
    print(f"  {contig}")


if __name__ == "__main__":
    main()
