#!bin/usr/env python3

"""Tools for converting input data into standard data format for NP-MRD data exchange"""

from rdkit import Chem
from pathlib import Path


def generate_mol_index(smiles):
    """"Create list of mol line index numbers for molecule"""
    mol = Chem.MolFromSmiles(smiles)

    output_path = Path("/Users/roger/Git/rgl_utilities/data/np_mrd/test_mol.mol")
    mol_block = Chem.MolToMolBlock(mol)
    mol2 = Chem.MolFromMolFile("/Users/roger/Git/rgl_utilities/data/np_mrd/test_mol.mol")
    print(Chem.MolToSmiles(mol2))

    for atom in mol2.GetAtoms():
        try:
            print(f"Atomic number: {atom.GetAtomicNum()} "
                  f"rdkit index: {atom.GetIdx()} "
                  f"mol block index: {atom.GetProp('molFileAlias')}")
        except KeyError:
            print(f"MISSING MOL INDEX "
                  f"Atomic number: {atom.GetAtomicNum()} "
                  f"rdkit index: {atom.GetIdx()} ")


if __name__ == "__main__":
    smiles = "CCCCO"
    generate_mol_index(smiles)