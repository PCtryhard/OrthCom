from Bio.PDB import PDBIO, PPBuilder, Structure, Model, Chain, Residue
import numpy as np

# Define the phi and psi angles for an alpha-helix
phi_angles = [-57] * 15
psi_angles = [-47] * 15

# Create a new structure
structure = Structure.Structure("Ideal_15-mer_Poly-Alanine_Alpha-Helix")

# Create a new model
model = Model.Model(0)
structure.add(model)

# Create a new chain
chain = Chain.Chain("A")
model.add(chain)

# Create the poly-alanine residues
for i in range(15):
    # Create a new residue
    residue = Residue.Residue("ALA", i + 1, " ")
    chain.add(residue)

# Use the PPBuilder to assign phi and psi angles
ppb = PPBuilder.PPBuilder()
ppb.set_phi(phi_angles)
ppb.set_psi(psi_angles)
ppb.build_peptides(chain)

# Calculate the coordinates
ppb.ppbuild(chain)

# Save the structure to a PDB file
io = PDBIO()
io.set_structure(structure)
with open("data/raw/helix_scaffold.pdb", "w") as pdb_file:
    io.save(pdb_file)
