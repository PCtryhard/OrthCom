import json
import glob
import os


def prep_multimer_fastas(results_dir="data/processed/alphafold", seq_dir="data/processed/mpnn/seqs",
                         out_dir="data/processed/multimer_input"):
    # ==========================================
    # PASTE YOUR ACTUAL KEY SEQUENCE HERE
    # ==========================================
    KEY_SEQUENCE = "REPLACE_ME_WITH_YOUR_PEPTIDE_SEQUENCE"

    os.makedirs(out_dir, exist_ok=True)
    results = []

    # 1. Parse the scores (same logic as before)
    for filepath in glob.glob(f"{results_dir}/*_scores_*.json"):
        with open(filepath, 'r') as f:
            try:
                data = json.load(f)
                if 'plddt' in data:
                    avg_plddt = sum(data['plddt']) / len(data['plddt'])
                    name = os.path.basename(filepath).split('_scores_')[0]
                    results.append({"name": name, "plddt": avg_plddt})
            except Exception as e:
                pass

    # Sort and grab Top 10
    results.sort(key=lambda x: x["plddt"], reverse=True)
    top_10 = results[:10]

    print("🧬 Prepping Top 10 Locks for Multimer...")

    # 2. Extract original sequences and write Multimer FASTAs
    for rank, res in enumerate(top_10):
        lock_name = res["name"]
        seq_file = f"{seq_dir}/{lock_name}.fa"

        try:
            with open(seq_file, 'r') as f:
                # Read the sequence, skipping the >header line
                lock_seq = "".join([line.strip() for line in f.readlines() if not line.startswith(">")])

            # Multimer format: LockSequence:KeySequence
            multimer_fasta = f">multimer_rank{rank + 1}_{lock_name}\n{lock_seq}:{KEY_SEQUENCE}\n"

            out_file = f"{out_dir}/rank{rank + 1}_{lock_name}.fa"
            with open(out_file, 'w') as f:
                f.write(multimer_fasta)

            print(f"✅ Created: rank{rank + 1}_{lock_name}.fa (pLDDT: {res['plddt']:.1f})")

        except FileNotFoundError:
            print(f"⚠️ Could not find sequence file for {lock_name}")

    print(f"\nDone! Multimer files are waiting in {out_dir}/")


if __name__ == "__main__":
    prep_multimer_fastas()