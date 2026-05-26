#!/bin/bash

mkdir -p data/processed/mpnn

# 1. Run ProteinMPNN (Your command is already perfect for Chain Freezing)
for file in data/processed/rfdiffusion/*.pdb; do
    echo "Running ProteinMPNN on $file..."
    python ProteinMPNN/protein_mpnn_run.py --pdb_path "$file" --pdb_path_chains "B" --out_folder data/processed/mpnn --num_seq_per_target 10 --sampling_temp 0.0001 --batch_size 1
done

# 2. Automatically split FASTA files and isolate the Lock
echo "Splitting bundled FASTA files and extracting Chain B for AlphaFold..."
python -c '
import os, glob

seq_dir = "data/processed/mpnn/seqs"

for f in glob.glob(f"{seq_dir}/*.fa"):
    with open(f, "r") as file:
        entries = file.read().strip().split(">")[1:]

    base_name = f.replace(".fa", "")

    # Skip entry 0 (the native RFdiffusion poly-glycine sequence)
    for i, entry in enumerate(entries[1:]):
        lines = entry.strip().split("\n")
        header = lines[0]

        # ProteinMPNN separates chains with a slash. Grab the last part (Chain B)
        full_seq = "".join(lines[1:]).replace(" ", "")
        lock_seq = full_seq.split("/")[-1]

        # Save as seq1 through seq10
        with open(f"{base_name}_seq{i+1}.fa", "w") as out:
            out.write(">" + header + "\n" + lock_seq + "\n")

    # Delete the original bundled file
    os.remove(f)
'
echo "Done! 1,000 individual Lock sequences are isolated and ready."