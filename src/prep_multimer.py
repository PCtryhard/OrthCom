#!/usr/bin/env python3
"""Prepare ColabFold AF2-multimer FASTAs: pair each designed Lock (chain A)
sequence with the frozen Key (chain B) so the *interface* (iPTM / interface-PAE)
can be scored, not the monomer fold.

The Key sequence is read programmatically from the canonical complex
(scaffold_clean.pdb chain B) so it can never drift from what was actually
diffused. ColabFold separates chains inside one FASTA record with a ':'
delimiter -- verified against the installed colabfold/input.py
(query_sequence.upper().split(":")).

Usage:
    python src/prep_multimer.py --lock-seq-dir data/processed/top10
    python src/prep_multimer.py --lock-seqs a.fa b.fa --out-dir data/processed/multimer_input
"""

import argparse
import glob
import os

AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
          "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
          "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
          "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def read_chain_sequence(pdb_path, chain_id):
    """1-letter sequence of one chain, ordered by CA-atom appearance."""
    seq = []
    with open(pdb_path) as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA" \
                    and ln[21] == chain_id:
                seq.append(AA3TO1.get(ln[17:20].strip(), "X"))
    if not seq:
        raise ValueError(f"No CA atoms for chain {chain_id} in {pdb_path}")
    return "".join(seq)


def read_lock_sequence(fa_path):
    """Lock (chain A) sequence from an MPNN .fa file (header skipped)."""
    with open(fa_path) as fh:
        body = "".join(line.strip() for line in fh if not line.startswith(">"))
    body = body.replace(" ", "")
    # Defensive: these files are already chain-A-only, but if a chain break ever
    # leaks in, keep only the first (Lock) chain.
    for delim in (":", "/"):
        if delim in body:
            body = body.split(delim)[0]
    return body


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--lock-seq-dir",
                     help="Directory of Lock .fa files (e.g. your isolated top-N).")
    src.add_argument("--lock-seqs", nargs="+",
                     help="Explicit list of Lock .fa files.")
    ap.add_argument("--out-dir", default="data/processed/multimer_input")
    ap.add_argument("--key-pdb", default="data/raw/scaffold_clean.pdb",
                    help="Structure to read the frozen Key sequence from.")
    ap.add_argument("--key-chain", default="B")
    args = ap.parse_args()

    key_seq = read_chain_sequence(args.key_pdb, args.key_chain)
    print(f"Key  (chain {args.key_chain} of {args.key_pdb}): "
          f"{key_seq}  [{len(key_seq)} res]")

    if args.lock_seq_dir:
        lock_files = sorted(glob.glob(os.path.join(args.lock_seq_dir, "*.fa")))
    else:
        lock_files = args.lock_seqs
    if not lock_files:
        raise SystemExit("No Lock .fa files found.")

    os.makedirs(args.out_dir, exist_ok=True)
    written = 0
    for fa in lock_files:
        name = os.path.splitext(os.path.basename(fa))[0]
        try:
            lock_seq = read_lock_sequence(fa)
        except OSError as exc:
            print(f"  [warn] skip {name}: {exc}")
            continue
        if not lock_seq:
            print(f"  [warn] skip {name}: empty Lock sequence")
            continue
        out = os.path.join(args.out_dir, f"{name}.fa")
        with open(out, "w") as fh:
            # One record, two chains: Lock ':' Key  (filename stem == jobname).
            fh.write(f">{name}\n{lock_seq}:{key_seq}\n")
        written += 1
        print(f"  wrote {name}.fa   Lock={len(lock_seq)} : Key={len(key_seq)}")

    print(f"\n{written} multimer FASTA(s) -> {args.out_dir}")
    print("Next: bash src/run_multimer.sh")


if __name__ == "__main__":
    main()
