#!/usr/bin/env python3
"""Rank co-folded Lock+Key complexes by INTERFACE quality -- not monomer fold.

Sort key: interface PAE (ascending; field success threshold < 10), then iPTM
(descending) as the tie-break.

interface PAE = mean of the PAE matrix restricted to the two CROSS-CHAIN blocks
    rows[Lock] x cols[Key]   AND   rows[Key] x cols[Lock]
i.e. both off-diagonal blocks only -- NOT the full matrix, and explicitly NOT
monomer pLDDT. (PAE is asymmetric, so both directions are averaged.)

Operates on ColabFold AF2-multimer outputs from src/run_multimer.sh
(data/processed/multimer/*_scores_*.json), which carry 'pae', 'iptm', 'ptm',
'plddt'. The Lock/Key boundary is recovered from the LockSeq:KeySeq FASTA that
was folded, and printed so it can be sanity-checked.
"""

import argparse
import glob
import json
import os


def lock_key_lengths(name, multimer_dir):
    """(L_lock, K_key) from the folded FASTA (LockSeq:KeySeq), or None."""
    fa = os.path.join(multimer_dir, f"{name}.fa")
    if not os.path.exists(fa):
        return None
    with open(fa) as fh:
        body = "".join(line.strip() for line in fh if not line.startswith(">"))
    if ":" not in body:
        return None
    parts = body.split(":")
    return len(parts[0]), len(parts[1])


def interface_pae(pae, lock_len, n):
    """Mean PAE over both cross-chain blocks given Lock=res[0:lock_len]."""
    total, count = 0.0, 0
    for i in range(n):
        row = pae[i]
        if i < lock_len:                 # Lock row -> Key columns
            for j in range(lock_len, n):
                total += row[j]
                count += 1
        else:                            # Key row -> Lock columns
            for j in range(lock_len):
                total += row[j]
                count += 1
    return (total / count) if count else float("nan")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="data/processed/multimer",
                    help="Dir of AF2-multimer *_scores_*.json complex outputs.")
    ap.add_argument("--multimer-dir", default="data/processed/multimer_input",
                    help="Dir of the LockSeq:KeySeq FASTAs (for the chain boundary).")
    ap.add_argument("--key-len", type=int, default=15,
                    help="Fallback Key length if the FASTA is missing (Lock = N - key_len).")
    ap.add_argument("--threshold", type=float, default=10.0,
                    help="interface-PAE success threshold (Angstrom).")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.results_dir, "*_scores_*.json")))
    if not files:
        print(f"No *_scores_*.json in {args.results_dir}. Run src/run_multimer.sh first.")
        return

    rows = []
    boundaries_seen = set()
    for fp in files:
        name = os.path.basename(fp).split("_scores_")[0]
        try:
            with open(fp) as fh:
                data = json.load(fh)
        except (OSError, ValueError) as exc:
            print(f"  [warn] skip {os.path.basename(fp)}: {exc}")
            continue

        pae = data.get("pae")
        iptm = data.get("iptm")
        plddt = data.get("plddt")
        if pae is None or iptm is None:
            missing = "pae" if pae is None else "iptm"
            print(f"  [warn] skip {name}: no '{missing}' "
                  f"(folded with a multimer model?)")
            continue
        n = len(pae)

        lk = lock_key_lengths(name, args.multimer_dir)
        if lk is not None and (lk[0] + lk[1] == n):
            lock_len, key_len, via = lk[0], lk[1], "fasta"
        else:
            if lk is not None:
                print(f"  [warn] {name}: FASTA Lock+Key={lk[0]}+{lk[1]} != PAE N={n}; "
                      f"using N-key_len fallback")
            lock_len, key_len, via = n - args.key_len, args.key_len, "fallback"

        if (lock_len, n) not in boundaries_seen:
            boundaries_seen.add((lock_len, n))
            print(f"[boundary] {name}: Lock=res[0:{lock_len}]  Key=res[{lock_len}:{n}]  "
                  f"(L={lock_len}, K={key_len}, N={n}, via {via})")

        ipae = interface_pae(pae, lock_len, n)
        complex_plddt = sum(plddt) / len(plddt) if plddt else float("nan")
        key_plddt = (sum(plddt[lock_len:]) / key_len
                     if plddt and key_len else float("nan"))
        rows.append({"name": name, "ipae": ipae, "iptm": iptm,
                     "cplddt": complex_plddt, "kplddt": key_plddt})

    if not rows:
        print("No usable complex scores parsed.")
        return

    # Primary: interface PAE ascending. Tie-break: iPTM descending.
    rows.sort(key=lambda r: (r["ipae"], -r["iptm"]))

    print("\n" + "=" * 72)
    print(" INTERFACE RANKING  (sort: interface-PAE asc, then iPTM desc)")
    print("=" * 72)
    print(f"{'#':>2}  {'design':<22}{'iPAE':>7}{'iPTM':>7}{'cplDDT':>8}{'keyplDDT':>10}   ")
    print("-" * 72)
    for i, r in enumerate(rows, 1):
        flag = "PASS" if r["ipae"] < args.threshold else ""
        print(f"{i:>2}  {r['name']:<22}{r['ipae']:>7.2f}{r['iptm']:>7.3f}"
              f"{r['cplddt']:>8.1f}{r['kplddt']:>10.1f}   {flag}")
    n_pass = sum(1 for r in rows if r["ipae"] < args.threshold)
    print("-" * 72)
    print(f"{n_pass}/{len(rows)} designs below interface-PAE {args.threshold} "
          f"(lower iPAE = better; higher iPTM = better)")


if __name__ == "__main__":
    main()
