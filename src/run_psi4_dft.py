import psi4
import rdkit
from rdkit import Chem
from pathlib import Path

def main() -> None:
    # Load the 3D molecule from the mol file
    mol_path = Path('data/raw/L_4_fluorophenylalanine.mol')
    mol = Chem.MolFromMolFile(mol_path, removeHs=False)

    # Extract the element symbols and their 3D coordinates
    geom_str = '0 1\n'
    for conf in mol.GetConformers():
        for i, atom in enumerate(mol.GetAtoms()):
            x, y, z = conf.GetAtomPosition(atom.GetIdx())
            geom_str += f'{atom.GetSymbol()} {x:.6f} {y:.6f} {z:.6f}\n'

    # Set the Psi4 output file
    psi4.core.set_output_file('data/raw/psi4_output.log')

    # Configure Psi4 to use the B3LYP functional and the 6-31G* basis set
    psi4.set_options({'basis': '6-31G*', 'reference': 'rhf'})

    # Create the Psi4 molecule
    mol = psi4.geometry(geom_str)

    # Run a single-point energy calculation
    energy, wfn = psi4.energy('b3lyp', return_wfn=True)

    # Print the final DFT energy to the console
    print(f'The final DFT energy is: {energy}')

if __name__ == '__main__':
    main()
