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

# 2. Load the ncAA via RDKit (keep Hydrogens)
mol = Chem.MolFromMolFile("data/raw/L_4_fluorophenylalanine.mol", removeHs=False)
conf = mol.GetConformer()

# 3. Extract matching backbone atoms (N, CA, C) for alignment
# In BioPython
bp_atoms = [target_residue["N"], target_residue["CA"], target_residue["C"]]

# In RDKit: Find the backbone atoms by matching atomic numbers and environments
# L-4-fluorophenylalanine backbone map
rd_n_idx, rd_ca_idx, rd_c_idx = None, None, None
for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()
    # Find the amine Nitrogen
    if symbol == "N" and atom.GetTotalNumHs() >= 2:
        rd_n_idx = atom.GetIdx()
    elif symbol == "C":
        # Alpha Carbon is bonded to N and a carboxyl C
        neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
        if "N" in neighbors and "O" in [n.GetSymbol() for n in atom.GetNeighbors()]:
            # This is the carbonyl carbon or alpha carbon. Let's distinguish:
            pass

# To make this foolproof without fragile graph matching, we map by known SMILES order
# The SMILES used was: C(O)(=O)[C@H](CC1=CC=C(F)C=C1)N
# Let's map coordinates directly based on standard atom indexing
rd_coords = conf.GetPositions()

# For our specific generated molecule, let's find backbone by neighbor counts:
for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    if atom.GetSymbol() == "N" and len([n for n in atom.GetNeighbors() if n.GetSymbol() == "C"]) == 1:
        rd_n_idx = idx
    elif atom.GetSymbol() == "C":
        nbors = [n.GetSymbol() for n in atom.GetNeighbors()]
        # CA is bonded to N, C (carbonyl), and C (beta)
        if "N" in nbors and nbors.count("C") == 2:
            rd_ca_idx = idx
        # Carbonyl C is bonded to CA, O, O
        elif "O" in nbors and "C" in nbors:
            rd_c_idx = idx

# 4. Compute and apply alignment matrix
sup = Superimposer()
# BioPython atoms as fixed targets
fixed_coords = np.array([atom.get_coord() for atom in bp_atoms])
# RDKit moving atoms
moving_coords = np.array([rd_coords[rd_n_idx], rd_coords[rd_ca_idx], rd_coords[rd_c_idx]])

sup.set(fixed_coords, moving_coords)

# Rotate and translate all coordinates of the RDKit molecule
transformed_coords = np.dot(rd_coords, sup.rotran[0]) + sup.rotran[1]

# 5. Rebuild residue 8 into our custom L4F residue
new_residue = Residue.Residue((" ", 8, " "), "L4F", " ")

# Add transformed atoms into the BioPython residue object
atom_counts = {}
for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()
    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    atom_name = f"{symbol}{atom_counts[symbol]}"
    
    # Standardize backbone names so downstream tools don't break
    if atom.GetIdx() == rd_n_idx: atom_name = "N"
    elif atom.GetIdx() == rd_ca_idx: atom_name = "CA"
    elif atom.GetIdx() == rd_c_idx: atom_name = "C"
    
    coord = transformed_coords[atom.GetIdx()]
    bp_atom = Atom.Atom(atom_name, coord, 60.0, 1.0, " ", atom_name, atom.GetIdx(), symbol)
    new_residue.add(bp_atom)

# 6. Swap old residue 8 with our new mutated residue
chain.detach_child((" ", 8, " "))
chain.add(new_residue)

# 7. Write out the final physical Key PDB
io = PDBIO()
io.set_structure(structure)
output_path = Path("data/raw/synthetic_key.pdb")
io.save(str(output_path))
print("Successfully generated data/raw/synthetic_key.pdb")
