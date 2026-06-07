#!/usr/bin/env python3
"""
analyze_ligandmpnn_dryrun.py  —  read the 4-condition LigandMPNN dry-run stats
and report the discrimination-control mechanics.

Conditions (out_folder subdirs under outputs/_dryB):
  L4F_on, L4F_off, PHE_on, PHE_off   =  {F-present, F-absent} x {atom-context on/off}

Reports, at the redesigned pocket positions:
  - which positions were designed (chain_mask), mapped to chain-A resseq
  - per-condition: AA sampled across the 8 batches + mean P(AA)
  - CONTROL  : L4F_off vs PHE_off should be identical (ligand ablated)
  - SIGNAL   : L4F_on  vs PHE_on  divergence  (does the model sense the F?)
  - CONTEXT  : L4F_on  vs L4F_off divergence  (does the ligand matter at all?)
F-proximal positions (<=5 A of the para-F) are flagged with *.
"""
import os
import sys
import glob
import numpy as np
import torch

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "LigandMPNN", "outputs", "_dryB")
RESSEQ0 = 6                     # chain A starts at resseq 6 -> index 0
FPROX = {52, 68, 70, 79}       # from prep_ligandmpnn_input.py


def load(tag):
    f = glob.glob(os.path.join(BASE, tag, "stats", "*.pt"))[0]
    return torch.load(f, map_location="cpu")


def js_div(p, q, eps=1e-9):
    p = p + eps; q = q + eps
    p = p / p.sum(); q = q / q.sum()
    m = 0.5 * (p + q)
    kl = lambda a, b: np.sum(a * np.log2(a / b))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def main():
    cond = {t: load(t) for t in ("L4F_on", "L4F_off", "PHE_on", "PHE_off")}
    ref = cond["L4F_on"]
    designed = [i for i in range(ref["chain_mask"].shape[0]) if int(ref["chain_mask"][i]) == 1]
    resseq = [RESSEQ0 + i for i in designed]
    print(f"designed positions (chain_mask==1): {len(designed)}")
    print(f"  index : {designed}")
    print(f"  resseq: A{', A'.join(map(str, resseq))}")

    # decoding-order identical across conditions? (same seed)
    same_order = all(torch.equal(cond[t]["decoding_order"], ref["decoding_order"])
                     for t in cond)
    print(f"decoding_order identical across all 4 conditions: {same_order}")

    def counts(d, idx):
        seqs = d["generated_sequences"][:, idx].numpy()        # (8,)
        aas = [ALPHABET[s] for s in seqs]
        uniq = sorted(set(aas), key=lambda a: -aas.count(a))
        return "".join(aas), " ".join(f"{a}:{aas.count(a)}" for a in uniq)

    def meanprob(d, idx):
        return d["sampling_probs"][:, idx, :].mean(0).numpy()   # (20,)

    # ---- CONTROL: off conditions must be identical ----
    off_ident = torch.equal(cond["L4F_off"]["generated_sequences"],
                            cond["PHE_off"]["generated_sequences"])
    print(f"\n[CONTROL]  L4F_off == PHE_off generated_sequences (whole chain): {off_ident}")
    print("           (ablating atom context should erase any F-vs-noF difference)")

    # ---- per-position table ----
    print("\n[PER-POSITION]  sampled AA across 8 batches  (counts);  * = F-proximal")
    hdr = f"{'pos':>5} {'F?':>2}  {'L4F_on':>14} {'PHE_on':>14} {'L4F_off':>14} {'PHE_off':>14}"
    print(hdr); print("-" * len(hdr))
    for idx, rs in zip(designed, resseq):
        flag = "*" if rs in FPROX else " "
        row = f"A{rs:<4}{flag:>1}  "
        for t in ("L4F_on", "PHE_on", "L4F_off", "PHE_off"):
            s, _ = counts(cond[t], idx)
            row += f"{s:>14} "
        print(row)

    # ---- SIGNAL & CONTEXT divergences ----
    print("\n[DIVERGENCE]  Jensen-Shannon (bits, 0..1) of mean P(AA) over 8 batches")
    hdr2 = f"{'pos':>5} {'F?':>3} {'SIGNAL L4F_on|PHE_on':>22} {'CONTEXT L4F_on|L4F_off':>24}"
    print(hdr2); print("-" * len(hdr2))
    sig_fprox, sig_other = [], []
    for idx, rs in zip(designed, resseq):
        flag = "*" if rs in FPROX else " "
        sig = js_div(meanprob(cond["L4F_on"], idx), meanprob(cond["PHE_on"], idx))
        ctx = js_div(meanprob(cond["L4F_on"], idx), meanprob(cond["L4F_off"], idx))
        print(f"A{rs:<4}{flag:>3} {sig:>22.4f} {ctx:>24.4f}")
        (sig_fprox if rs in FPROX else sig_other).append(sig)
    print("-" * len(hdr2))
    if sig_fprox:
        print(f"mean SIGNAL  F-proximal (*): {np.mean(sig_fprox):.4f}")
    if sig_other:
        print(f"mean SIGNAL  non-proximal  : {np.mean(sig_other):.4f}")

    # ---- top-AA mean-prob detail at F-proximal positions ----
    print("\n[DETAIL]  mean P(AA) top-3 at F-proximal positions (L4F_on vs PHE_on)")
    for idx, rs in zip(designed, resseq):
        if rs not in FPROX:
            continue
        for t in ("L4F_on", "PHE_on"):
            p = meanprob(cond[t], idx)
            top = np.argsort(-p)[:3]
            s = ", ".join(f"{ALPHABET[j]}={p[j]:.2f}" for j in top)
            print(f"  A{rs:<4} {t:>8}: {s}")


if __name__ == "__main__":
    main()
