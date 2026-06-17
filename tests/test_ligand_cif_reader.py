"""Regression tests for the ligand monomer-dictionary reader.

Guards the bond-order-column bug fixed by reading dicts via gemmi's ChemComp
parser: modern acedrg/PDBx dictionaries name the bond-order column
``_chem_comp_bond.value_order`` (with ``pdbx_aromatic_flag``), while older
refmac/grade dictionaries use ``_chem_comp_bond.type`` (with ``aromatic``). The
previous hand-rolled reader only looked for ``.type``, so on an acedrg dict the
bond loop came back empty and the molecule was built with NO bonds -- a bag of
disconnected atoms that RDKit then "embedded" all at the origin, collapsing the
modelled ligand to a single point. These tests assert both dialects yield a
fully-bonded, connected molecule.
"""
import os
import tempfile

import pytest
from rdkit import Chem

from pandda_gemmi.dataset.small import get_fragment_mol_from_dataset_cif_path

# A minimal acetate-like ligand: C1-C2, C2=O1 (double), C2-O2 (single).
_HEADER = """data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
LIG LIG
data_comp_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
LIG C1 C 0
LIG C2 C 0
LIG O1 O 0
LIG O2 O -1
"""

# Modern acedrg / PDBx dialect: value_order + pdbx_aromatic_flag.
_ACEDRG_DICT = _HEADER + """loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
LIG C1 C2 SINGLE N
LIG C2 O1 DOUBLE N
LIG C2 O2 SINGLE N
"""

# Older refmac / grade dialect: type + aromatic.
_REFMAC_DICT = _HEADER + """loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.aromatic
LIG C1 C2 single n
LIG C2 O1 double n
LIG C2 O2 single n
"""


def _read_dict(text):
    with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as f:
        f.write(text)
        path = f.name
    try:
        return get_fragment_mol_from_dataset_cif_path(path)
    finally:
        os.unlink(path)


@pytest.mark.parametrize(
    "dialect,text",
    [("acedrg_value_order", _ACEDRG_DICT), ("refmac_type", _REFMAC_DICT)],
)
def test_reader_builds_bonded_connected_mol(dialect, text):
    mol = _read_dict(text)
    assert mol.GetNumAtoms() == 4, dialect
    # The bug: acedrg dict -> 0 bonds. Must be fully bonded...
    assert mol.GetNumBonds() == 3, f"{dialect}: got {mol.GetNumBonds()} bonds, expected 3"
    # ...and a single connected fragment (not a bag of disconnected atoms).
    assert len(Chem.GetMolFrags(mol)) == 1, f"{dialect}: expected one connected fragment"
    # The declared C=O double bond order is preserved.
    n_double = sum(
        1 for b in mol.GetBonds()
        if b.GetBondType() == Chem.rdchem.BondType.DOUBLE
    )
    assert n_double == 1, f"{dialect}: expected the declared C=O double bond"


def test_inbuilt_reader_is_the_single_dataset_small_implementation():
    """The duplicate reader in autobuild.inbuilt is gone; it re-exports the one
    in dataset.small (so fixing the reader fixes autobuild too)."""
    from pandda_gemmi.autobuild import inbuilt
    from pandda_gemmi.dataset import small
    assert (inbuilt.get_fragment_mol_from_dataset_cif_path
            is small.get_fragment_mol_from_dataset_cif_path)
