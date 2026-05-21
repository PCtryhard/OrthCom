from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def main() -> None:
    smiles = "C(O)(=O)[C@H](CC1=CC=C(F)C=C1)N"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    embed_result = AllChem.EmbedMolecule(mol, params)
    if embed_result == -1:
        raise RuntimeError("Failed to generate 3D coordinates")

    AllChem.MMFFOptimizeMolecule(mol)

    output_path = Path("data/raw/L_4_fluorophenylalanine.mol")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    Chem.MolToMolFile(mol, str(output_path))


if __name__ == "__main__":
    main()
