#!/usr/bin/env python3

"""Class definition for compounds for NMR project"""
from pathlib import Path
import csv
import os

from rdkit import Chem
from rdkit.Chem import Draw


ATOMIC_NUMBERS = {0: "Unk",
                  1: "H",
                  5: "B",
                  6: "C",
                  7: "N",
                  8: "O",
                  9: "F",
                  11: "Na",
                  12: "Mg",
                  13: "Al",
                  14: "Si",
                  15: "P",
                  16: "S",
                  17: "Cl",
                  19: "K",
                  26: "Fe",
                  29: "Cu",
                  30: "Zn",
                  31: "Ga",
                  33: "As",
                  34: "Se",
                  35: "Br",
                  53: "I"}

BOND_TYPES = {0.0: "unknown",
              1.0: "single",
              1.5: "aromatic",
              2.0: "double",
              3.0: "triple"}

BOND_TYPE_ABBREVIATIONS = {0.0: "U",
                           1.0: "S",
                           1.5: "A",
                           2.0: "D",
                           3.0: "T"}


class Compound:
    def __init__(self, source, source_id, source_type, source_data, smiles, h_shifts=None, c_shifts=None, solvent=None):
        """Constructs a new molecule object

        Args:
            source (pathlib.Path): a Path object that defines the origin of the data contained in the compound object
            source_id (str): the original id code for the molecule from the source file
            source_type (str): description of source of data. Examples include 'jeol' and 'curator'
            source_data (dict): the original data from the source file
            smiles (str): a SMILES string that provides a machine-readable representation of the chemical structure
            mol (rdkit.Chem.Mol): rdkit mol object depicting the chemical structure
            h_shifts (dict): a dictionary containing the standardized np-mrd atom indices and h chemical shifts as
            key, value. Obtain np-mrd atom index using atom.GetAtomMapNum()
            c_shifts (dict): a dictionary containing the standardized np-mrd atom indices and c chemical shifts as
            key, value
            solvent (str): the NMR solvent in which the data were acquired
            c_shift_error (dict): a dictionary containing the standardized np-mrd atom indices and number of standard
            deviations the 13C shift is away from the mean value for carbons in the same direct environment
            (same adjacent atoms and bond types)
            direct_neighbor_labels (dict): a dictionary containing the standardized np-mrd atom indices and the hash
            describing the adjacent atoms and bond types (r=1)
            extended_neighbor_labels (dict): a dictionary containing the standardized np-mrd atom indices and the hash
            describing the adjacent atoms and bond types (r=2)
        """
        self.source = source
        self.source_id = source_id
        self.source_type = source_type
        self.source_data = source_data
        self.smiles = smiles
        self.mol = self._rdkit_atom_order(Chem.MolFromSmiles(self.smiles))
        self.h_shifts = h_shifts
        self.c_shifts = c_shifts
        self.solvent = solvent
        self.c_shift_error = None
        self.direct_neighbor_labels = self._add_direct_neighbor_labels()
        self.extended_neighbor_labels = self._add_extended_neighbor_labels()

        if self.c_shifts:
            self.add_13c_atom_labels()

    def _rdkit_atom_order(self, m, add_hs=False):
        """Canonicalize using RDKit SMILES export

        Args:
            m (rdkit.Chem.Mol): Mol object for RDKit

        Returns:
            rdkit.Chem.Mol: New canonicalized RDKit mol
        """
        m_renum = Chem.MolFromSmiles(Chem.MolToSmiles(m))

        # if self.source_type == "curator":
        #     add_hs = True
        if add_hs:
            m_canon = Chem.AddHs(m_renum)
        else:
            m_canon = m_renum
        for i, a in enumerate(m_canon.GetAtoms()):
            a.SetAtomMapNum(i + 1)

        return m_canon

    def add_13c_atom_labels(self):
        """Add NP-MRD standardized atom numbers and 13C chemical shifts to all carbons in molecule

        Args:
            mol (rdkit.Chem.Mol): Mol object for RDKit

        Returns:
            rdkit.Chem.Mol: New labeled RDKit mol
        """

        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                # NOTE: it is important to use 'GetAtomMapNum' rather than 'GetIdx' to retrieve index values that map to
                # the curator output because np-mrd indexing starts at 1, not 0 (as is standard for rdkit)
                atom_index = atom.GetAtomMapNum()
                try:
                    atom.SetProp("atomLabel", f"{atom_index}, {round(self.c_shifts[atom_index], 1)}")
                except KeyError:
                    atom.SetProp("atomLabel", f"{atom_index}, None")

    def _add_direct_neighbor_labels(self):
        """Create dictionary of np-mrd indices and direct neighbor labels"""

        direct_neighbor_labels = {}
        for atom in self.mol.GetAtoms():
            direct_neighbor_labels[atom.GetAtomMapNum()] = atom_fixed_hash(self.mol, atom.GetIdx())

        return direct_neighbor_labels

    def _add_extended_neighbor_labels(self):
        """Create dictionary of np-mrd indices and direct neighbor labels"""

        direct_neighbor_labels = {}
        for atom in self.mol.GetAtoms():
            direct_neighbor_labels[atom.GetAtomMapNum()] = atom_extended_full_hash(self.mol, atom.GetIdx())

        return direct_neighbor_labels

    def create_c_shift_sd_error(self, shift_error_dict: dict) -> dict:
        """Create dictionary of np-mrd index numbers and number of standard deviations away from bulk mean for each
        13C atom environment

        Args:
            shift_error_dict (dict): Dictionary of direct neigbor environments and mean and sd values for all carbons
            within that environment. Derives from analyze.atom_types.environment_shift_distributions

        Returns:
            dictionary containing np-mrd index numbers and number of standard deviations away from mean each 13C value
            is"""

        c_shift_error_dict = {}
        for atom in self.mol.GetAtoms():
            atom_label = self.direct_neighbor_labels[atom.GetAtomMapNum()]
            try:
                c_shift_value = self.c_shifts[atom.GetAtomMapNum()]
                if shift_error_dict[atom_label]["sd"] != 0:
                    c_shift_sd_error = abs(c_shift_value - shift_error_dict[atom_label]["mean"]) / \
                                       shift_error_dict[atom_label]["sd"]
                else:
                    c_shift_sd_error = 0
            except KeyError:
                c_shift_sd_error = None

            c_shift_error_dict[atom.GetAtomMapNum()] = c_shift_sd_error

        return c_shift_error_dict

    def create_c_shift_table(self, export_dir: Path):
        """Export 13C chemical shifts and standardized atom indices to file

        Args:
            export_dir (Path): Path to directory where shift tables are to be saved

        Returns:
            csv file containing np_mrd indices and 13C chemical shifts
        """

        headers = [["npmrd_index", "13C_shift"]]
        data = list(self.c_shifts.items())

        output_path = Path(export_dir, f"{self.source_id}_13c_shifts.csv")
        with open(output_path, "w", newline='', encoding='utf-8') as f:
            csv_f = csv.writer(f)
            csv_f.writerows(headers + data)

    def create_image(self, export_dir: Path):
        """Create an image file of a given molecule, including RDKit atom indices

        Args:
            export_dir (Path): Path to directory where image is to be stored

        Returns:
            Image file of structure with 13C labels added
            """

        output_path = Path(export_dir,
                           f"{self.source_id}.png")

        for atom in self.mol.GetAtoms():
            atom.SetProp("molAtomMapNumber", str(atom.GetIdx() + 1))

        Draw.MolToFile(self.mol, output_path, size=(2500, 2500))

    def create_shift_error_image(self, output_path: Path):
        """Create an image for each compound containing 13C shifts, annotating 13C shift with deviation from mean/sd

        Args:
            output_path (Path): Path to top level output directory, in which data for each dataset is stored

        Returns:
            Image file displaying structure with index numbers and 13C chemical shifts, with each carbon atom labeled
            for sd of value to bulk mean for that environment. Green = 0 sd, red = 3 sd (continuous scale), blue = no
            value available
        """

        if self.c_shifts and self.c_shift_error:
            atom_highlights = {}
            atom_radii = {}
            for atom in self.mol.GetAtoms():
                if atom.GetAtomicNum() == 6:
                    c_shift_sd_error = self.c_shift_error[atom.GetAtomMapNum()]
                    if c_shift_sd_error:
                        # Limit max. sd error value to set color gradient to fixed scale
                        if c_shift_sd_error > 3:
                            c_shift_sd_error = 3
                        atom_highlights[atom.GetIdx()] = [(round(c_shift_sd_error / 3, 2),
                                                           round(1 - (c_shift_sd_error / 3), 2),
                                                           0,
                                                           0.2)]
                        atom_radii[atom.GetIdx()] = 0.3
                        atom.SetProp("atomLabel",
                                     f"{atom.GetAtomMapNum()}, {round(self.c_shifts[atom.GetAtomMapNum()], 1)}")
                    else:
                        atom_highlights[atom.GetIdx()] = [(0, 0, 1, 0.2)]
                        atom_radii[atom.GetIdx()] = 0.3
                        atom.SetProp("atomLabel",
                                     f"{atom.GetAtomMapNum()}, {None}")
                else:
                    atom_highlights[atom.GetIdx()] = [(0, 0, 0, 0)]
                    atom_radii[atom.GetIdx()] = 0.3

            output_path = Path(output_path,
                               f"{self.source_type}_compounds",
                               "shift_error_images",
                               f"Compound_{self.source_id}.png")

            if not Path.exists(output_path.parent):
                os.mkdir(output_path.parent)

            d2d = Draw.rdMolDraw2D.MolDraw2DCairo(2500, 2500)
            d2d.DrawMoleculeWithHighlights(self.mol, "", atom_highlights, {}, atom_radii, {})
            d2d.FinishDrawing()
            d2d.GetDrawingText()
            d2d.WriteDrawingText(str(output_path))

        else:
            print(f"ERROR: No 13C shift data and/or sd data available for compound {self.source_id}")

def atom_fixed_hash(mol: Chem.Mol, atom_index: int) -> str:
    """Generate a fixed extended hash for a given atom that lists both the neighboring atoms and bond types, and the
    extended neighbors (i.e. radius = 2)

    The basic unit of the hash is a string containing information about the atom, and the atoms it is connected to.
    The first 3 characters are for the atom. These are: NNA where NN is the two digit atomic number, and A is a 0/1 bool
     denoting whether the atom is aromatic or not. The following characters are four sets of fields
     (one for each possible bond) separated by a dash, in the format NNAB where NN is the atomic number, A denotes
     aromatic atom, and B denotes bond type (S = single, D = double, T = triple, A = aromatic)"""

    target_atom = mol.GetAtomWithIdx(atom_index)
    atom_atomic_number = f"{target_atom.GetAtomicNum():02d}"
    atom_aromatic = 0
    neighbors_data = []
    for neighbor in target_atom.GetNeighbors():
        neighbor_aromatic = 0
        neighbor_atomic_number = f"{neighbor.GetAtomicNum():02d}"
        bond_type = BOND_TYPE_ABBREVIATIONS[mol.GetBondBetweenAtoms(atom_index,
                                                                    neighbor.GetIdx()).GetBondTypeAsDouble()]
        if bond_type == "A":
            atom_aromatic = 1
            neighbor_aromatic = 1
        neighbors_data.append(f"{neighbor_atomic_number}{neighbor_aromatic}{bond_type}")

    if len(neighbors_data) < 4:
        for i in range(len(neighbors_data), 4):
            neighbors_data.append("XXXX")
    atom_descriptor = [f"{atom_atomic_number}{atom_aromatic}"]
    full_hash = "-".join(atom_descriptor + sorted(neighbors_data))

    return full_hash


def atom_extended_full_hash(mol: Chem.Mol, atom_index: int) -> str:
    """Create extended has that includes both direct (r=1) atoms and the extended atoms (r=2) for each of these
    positions arranged in the same index order"""

    target_atom = mol.GetAtomWithIdx(atom_index)
    target_atom_hash = atom_fixed_hash(mol, atom_index)
    neighbors_data = []
    for neighbor in target_atom.GetNeighbors():
        neighbor_aromatic = 0
        neighbor_atomic_number = f"{neighbor.GetAtomicNum():02d}"
        bond_type = BOND_TYPE_ABBREVIATIONS[mol.GetBondBetweenAtoms(atom_index,
                                                                    neighbor.GetIdx()).GetBondTypeAsDouble()]
        if bond_type == "A":
            neighbor_aromatic = 1
        neighbors_data.append([f"{neighbor_atomic_number}{neighbor_aromatic}{bond_type}",
                               atom_fixed_hash(mol, neighbor.GetIdx())])
    if len(neighbors_data) < 4:
        for i in range(len(neighbors_data), 4):
            neighbors_data.append(["XXX", "XXX-XXXX-XXXX-XXXX-XXXX"])
    sorted_neighbors = sorted(neighbors_data, key=itemgetter(0, 1))
    output_data = [target_atom_hash]
    for neighbor_data in sorted_neighbors:
        output_data.append(neighbor_data[1])
    full_hash = "_".join(output_data)

    return full_hash
