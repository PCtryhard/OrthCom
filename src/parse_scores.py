import json
import glob
import os


def parse_alphafold_scores(results_dir="data/processed/alphafold"):
    results = []

    # Strictly search for the score files, ignoring configs and PAE files
    for filepath in glob.glob(f"{results_dir}/*_scores_*.json"):
        with open(filepath, 'r') as f:
            try:
                data = json.load(f)

                if 'plddt' in data:
                    avg_plddt = sum(data['plddt']) / len(data['plddt'])
                    ptm = data.get('ptm', 0)

                    # Clean up the massive filename to just get "lock_X_seqY"
                    name = os.path.basename(filepath).split('_scores_')[0]
                    results.append({"name": name, "plddt": avg_plddt, "ptm": ptm})
            except Exception as e:
                print(f"Skipping {filepath} due to error: {e}")

    # Sort the list from highest pLDDT to lowest
    results.sort(key=lambda x: x["plddt"], reverse=True)

    print("\n" + "=" * 40)
    print(" 🏆 TOP 10 SYNTHETIC LOCKS 🏆")
    print("=" * 40)

    if not results:
        print("No scores found yet! Let AlphaFold cook.")
        return

    # Print the top 10
    for i, res in enumerate(results[:10]):
        print(f"{i + 1}. {res['name']}")
        print(f"   pLDDT: {res['plddt']:.1f} | pTM: {res['ptm']:.3f}\n")


if __name__ == "__main__":
    parse_alphafold_scores()