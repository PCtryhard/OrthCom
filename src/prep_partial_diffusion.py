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


def build_loop_contig(scaffold_res, key_res, threshold: float):
    key_xyz = np.array(
        [a.coord for r in key_res for a in r.get_atoms()], dtype=float
    )
    melt_flags = [min_distance_to_key(r, key_xyz) < threshold for r in scaffold_res]

    segments, i, n = [], 0, len(scaffold_res)
    chain_id = scaffold_res[0].get_parent().id
    melt_total = 0
    while i < n:
        j = i
        while j < n and melt_flags[j] == melt_flags[i]:
            j += 1
        run = scaffold_res[i:j]
        length = len(run)
        if melt_flags[i]:
            segments.append(f"{length}-{length}")  # diffused
            melt_total += length
        else:
            segments.append(f"{chain_id}{run[0].id[1]}-{run[-1].id[1]}")  # frozen motif
        i = j
    return segments, melt_total


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
    else:
        segments, melt_total = build_loop_contig(
            scaffold_res, key_res, args.contact_threshold)

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
    print(f"Contig ({n_scaf + n_key} residues total) -> {args.contig_out}")
    print(f"  {contig}")


if __name__ == "__main__":
    main()
