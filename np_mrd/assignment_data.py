#!bin/usr/env python3

"""Tools for converting input data into standard data format for NP-MRD data exchange"""

from rdkit import Chem
from pathlib import Path
import json
import csv


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


def count_double_bonds(json_path: Path):
    """Scan all NP-MRD exchange jsons for NP-MRD Curator app data and count number of molecules containing double bonds.
    Optionally, count how many of these do not have their configurations defined.
    Note that by default double bonds in aromatic systems, and those in rings with eight atoms or fewer do not need
    stereochemical definition because there is only one possible configuration allowed based on ring strain"""

    compound_count = 0
    double_bond_count = 0
    undefined_double_bond_count = 0

    undefined_compounds = []

    for json_file in json_path.glob("*.json"):
        with open(json_file, "rb") as f:
            input_data = json.load(f)

        for compound in input_data:
            compound_count += 1
            double_bonds_undefined = False
            compound_smiles = compound["smiles"]
            if "=" in compound_smiles:
                double_bond_count += 1
                if "\\" not in compound_smiles and "/" not in compound_smiles:
                    mol = Chem.MolFromSmiles(compound_smiles)
                    si = Chem.FindPotentialStereo(mol)
                    for element in si:
                        if str(element.type) == "Bond_Double" and str(element.specified) == "Unspecified":
                            print(f"Type: {element.type}, "
                                  f"Which: {element.centeredOn}, "
                                  f"Specified: {element.specified}, "
                                  f"Descriptor: {element.descriptor}")
                            print(Chem.MolToSmiles(mol))
                            double_bonds_undefined = True
            if double_bonds_undefined:
                undefined_double_bond_count += 1
                undefined_compounds.append([compound["citation"]["doi"],
                                            compound["compound_name"],
                                            compound["smiles"]])

    print(f"Compound count: {compound_count}\n"
          f"Compounds with double bonds: {double_bond_count}\n"
          f"Compounds with undefined configurations in double bonds: {undefined_double_bond_count}")

    headers = [["doi", "compound_name", "undefined_smiles"]]
    output_path = Path("/Users/roger/Documents/NP_MRD/Backfilling/npmrd_curator/exchange_files",
                       "exchange_jsons_20240926_missing_stereochem.csv")
    with open(output_path, "w") as g:
        csv_g = csv.writer(g)
        csv_g.writerows(headers + undefined_compounds)


if __name__ == "__main__":

    json_path = Path("/Users/roger/Documents/NP_MRD/Backfilling/npmrd_curator/exchange_files/exchange_jsons_20240926")

    count_double_bonds(json_path)