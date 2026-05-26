import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Atom, Residue
from Bio.PDB.Superimposer import Superimposer
from rdkit import Chem

# 1. Load the alpha helix via BioPython
parser = PDBParser(QUIET=True)
structure = parser.get_structure("helix", "data/raw/helix_scaffold.pdb")
model = structure[0]
chain = model["A"]
target_residue = chain[8]  # BioPython uses 1-indexing; 8 is the exact center

# 2. Load the canonical PHE via RDKit (keep Hydrogens)
mol = Chem.MolFromMolFile("data/raw/phenylalanine.mol", removeHs=False)
conf = mol.GetConformer()

# 3. Extract matching backbone atoms for alignment
# In BioPython
bp_atoms = [target_residue["N"], target_residue["CA"], target_residue["C"]]

# In RDKit: Find the backbone atoms
rd_n_idx, rd_ca_idx, rd_c_idx = None, None, None
for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    if atom.GetSymbol() == "N" and len([n for n in atom.GetNeighbors() if n.GetSymbol() == "C"]) == 1:
        rd_n_idx = idx
    elif atom.GetSymbol() == "C":
        nbors = [n.GetSymbol() for n in atom.GetNeighbors()]
        if "N" in nbors and nbors.count("C") == 2:
            rd_ca_idx = idx
        elif "O" in nbors and "C" in nbors:
            rd_c_idx = idx

rd_coords = conf.GetPositions()

# 4. Create dummy BioPython atoms for the Superimposer
moving_atoms = [
    Atom.Atom("N", rd_coords[rd_n_idx], 60.0, 1.0, " ", "N", 1, "N"),
    Atom.Atom("CA", rd_coords[rd_ca_idx], 60.0, 1.0, " ", "CA", 2, "C"),
    Atom.Atom("C", rd_coords[rd_c_idx], 60.0, 1.0, " ", "C", 3, "C")
]

# 5. Compute and apply alignment matrix
sup = Superimposer()
sup.set_atoms(bp_atoms, moving_atoms)
transformed_coords = np.dot(rd_coords, sup.rotran[0]) + sup.rotran[1]

# 6. Rebuild residue 8 into canonical PHE
new_residue = Residue.Residue((" ", 8, " "), "PHE", " ")

# Add transformed atoms
atom_counts = {}
for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()
    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    atom_name = f"{symbol}{atom_counts[symbol]}"

    # Standardize backbone names
    if atom.GetIdx() == rd_n_idx:
        atom_name = "N"
    elif atom.GetIdx() == rd_ca_idx:
        atom_name = "CA"
    elif atom.GetIdx() == rd_c_idx:
        atom_name = "C"

    coord = transformed_coords[atom.GetIdx()]
    bp_atom = Atom.Atom(atom_name, coord, 60.0, 1.0, " ", atom_name, atom.GetIdx(), symbol)
    new_residue.add(bp_atom)

# 7. Swap old residue 8 with our new mutated residue
chain.detach_child((" ", 8, " "))
chain.add(new_residue)

# 8. Write out the final physical Key PDB
io = PDBIO()
io.set_structure(structure)
output_path = Path("data/raw/canonical_key.pdb")
io.save(str(output_path))
print(f"Successfully generated {output_path} with a true physical Phenylalanine.")