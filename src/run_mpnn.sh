#!/bin/bash

mkdir -p data/processed/mpnn

# 1. Run ProteinMPNN: DESIGN the Lock (chain A), FREEZE the Key (chain B).
#    --pdb_path_chains lists the chains to DESIGN; every other chain in the PDB
#    is held fixed as context. NEW topology is scaffold/Lock=A, Key=B, so we
#    design A and freeze B. (Under the OLD topology Lock was B; hence the flip.)
for file in data/processed/rfdiffusion/*.pdb; do
    echo "Running ProteinMPNN on $file..."
    python ProteinMPNN/protein_mpnn_run.py --pdb_path "$file" --pdb_path_chains "A" --out_folder data/processed/mpnn --num_seq_per_target 3 --sampling_temp 0.0001 --batch_size 1
done

# 2. Automatically split FASTA files and isolate the Lock
echo "Splitting bundled FASTA files and extracting Chain A (the Lock) for AlphaFold..."
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

        # ProteinMPNN separates chains with a slash, in PDB chain order (A/B).
        # The Lock is chain A -> take the FIRST part. (Old topology had Lock=B,
        # which is why this used to grab the last part.)
        full_seq = "".join(lines[1:]).replace(" ", "")
        lock_seq = full_seq.split("/")[0]

        # Save as seq1 .. seqN (N = num_seq_per_target)
        with open(f"{base_name}_seq{i+1}.fa", "w") as out:
            out.write(">" + header + "\n" + lock_seq + "\n")

    # Delete the original bundled file
    os.remove(f)
'
echo "Done! Individual Lock (chain A) sequences are isolated and ready."