from Bio.PDB import PDBParser, PDBIO, Superimposer
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Load the PDB file
parser = PDBParser()
structure = parser.get_structure("helix_scaffold", "data/raw/helix_scaffold.pdb")

# Load the RDKit molecule
mol = Chem.MolFromMolFile("data/raw/L_4_fluorophenylalanine.mol", removeHs=False)

# Identify the backbone atoms (N, CA, C) in the 8th residue of the BioPython helix
helix_backbone_atoms = []
for atom in structure[0][7].get_atoms():
    if atom.get_name() in ["N", "CA", "C"]:
        helix_backbone_atoms.append(atom.get_coord())

# Identify the corresponding backbone atoms in the RDKit molecule
rdkit_backbone_atoms = []
for atom in mol.GetAtoms():
    if atom.GetSymbol() in ["N", "C", "CA"]:
        rdkit_backbone_atoms.append(atom.GetPos())

# Extract the 3D coordinates of both sets of backbone atoms
helix_coords = np.array(helix_backbone_atoms)
rdkit_coords = np.array(rdkit_backbone_atoms)

# Use Bio.PDB.Superimposer to calculate the rotation and translation matrix
super = Superimposer()
super.set(helix_coords, rdkit_coords)
rot, trans = super.get_rotran()

# Apply this transformation to all atom coordinates in the RDKit molecule
transformed_coords = np.dot(rdkit_coords, rot) + trans

# Rebuild the BioPython structure
new_structure = structure[0][:7] + structure[0][8:]  # Keep residues 1-7 and 9-15
new_structure[0][7] = structure[0][7].detached_copy()  # Detach residue 8

# Add the transformed RDKit atoms to residue 8
for i, atom in enumerate(mol.GetAtoms()):
    new_atom = new_structure[0][7].add_atom(
        " ", atom.GetSymbol(), "", 1, transformed_coords[i]
    )
    new_atom.set_name(atom.GetSymbol())

# Rename the residue to L4F
new_structure[0][7].id = " L4F"

# Save the final hybrid structure as a PDB file
io = PDBIO()
io.set_structure(new_structure)
io.save("data/raw/synthetic_key.pdb")
