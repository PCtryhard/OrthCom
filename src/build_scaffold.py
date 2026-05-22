from PeptideBuilder import Geometry
import PeptideBuilder
from Bio.PDB import PDBIO
from pathlib import Path

# Set ideal alpha-helix geometry for Alanine
geo = Geometry.geometry("A")
geo.phi = -57
geo.psi_im1 = -47

# Initialize first residue, then loop to add 14 more (15 total)
structure = PeptideBuilder.initialize_res(geo)
for i in range(14):
    PeptideBuilder.add_residue(structure, geo)

# Save the PDB file safely to the data/raw directory
io = PDBIO()
io.set_structure(structure)
pdb_path = Path("data/raw/helix_scaffold.pdb")
io.save(str(pdb_path))
