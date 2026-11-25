"""
Constants used throughout KIMMDY
"""

import sys

nN_per_kJ_per_mol_nm = 0.001661
R = 8.31446261815324e-3  # [kJ K-1 mol-1]
KB = 1.380649e-23  # [J K-1]
H_PLANCK = 6.62607015e-34  # [J s]

OPTIONAL_CONFIG_PATHS = []
"""Paths that may or may not be defined in the input
config and are thus not meant to be checked for existence
or only checked if they are using in conjunction with a certain plugin/option.
"""

CONFIG_LOGS = {
    "infos": [],
    "warnings": [],
    "errors": [],
    "debugs": [],
}
"""The logger is set up with information from the config
the the config can't use the logger.
Instead it collects the logmessages and displays them at the end.
"""

FIELD_SIZE_LIMIT = sys.maxsize
"""Maximum size of a field when reading a RecipeCollection from csv file.
"""

MARK_STARTED = ".kimmdy_started"
MARK_DONE = ".kimmdy_done"
MARK_FINISHED = ".kimmdy_finished"
MARK_FAILED = ".kimmdy_failed"
MARK_REACION_TIME = ".kimmdy_reaction_time"
MARKERS = [MARK_STARTED, MARK_DONE, MARK_FAILED, MARK_FINISHED, MARK_REACION_TIME]

REACTION_GRO = ".kimmdy_reaction.gro"
REACTION_TRR = ".kimmdy_reaction.trr"
REACTION_EDR = ".kimmdy_reaction.edr"

ATOM_ID_FIELDS = {
    "atoms": [0, 5],  # atomnr, chargegroup
    "bonds": [0, 1],
    "angles": [0, 1, 2],
    "dihedrals": [0, 1, 2, 3],
    "pairs": [0, 1],
    "settles": [0],
    "exclusions": True,  # all fields at atom ids
    "position_restraints": [0],
    "dihedral_restraints": [0, 1, 2, 3],
}

FALLBACK_EDISSOC = 400.0  # [kJ mol-1]
"""Dissociation energy used by reaction plugins when the edissoc_kimmdy file
is missing or has no entry for a bond. Close to typical single bonds between
heavy atoms.
"""

# see https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html
FFFUNC = {
    "mult_proper_dihedral": "9",
    "mult_improper_dihedral": "4",
    "harmonic_bond": "1",
    "morse_bond": "3",
    "harmonic_angle": "1",
    "urey_bradley_angle": "5",
    "pair": "1",
    "coulomb_pair": "2",
}

RESNR_ID_FIELDS = {
    "atoms": [2],
}

REACTIVE_MOLECULEYPE = "Reactive"


SOLVENT_NAMES: list[str] = [
    "SOL",
    "WATER",
    "TIP3P",
    "TIP4P",
    "TIP4P-Ew",
    "SPC",
    "SPC/E",
]

ION_NAMES: list[str] = [
    "I",
    "F",
    "CA",
    "CL",
    "NA",
    "MG",
    "K",
    "RB",
    "CS",
    "LI",
    "ZN",
]

# compare to atom type perception paper (2006) doi:10.1016/j.jmgm.2005.12.005
ATOMTYPE_BONDORDER: dict[tuple, int]
"""
To determin if an atom is a radical.
Compare to atom type perception paper (2006) doi:10.1016/j.jmgm.2005.12.005
"""
ATOMTYPE_BONDORDER = {
    ("MG", "NA", "CO"): 0,
    (
        "H",
        "HW",
        "HO",
        "HS",
        "HA",
        "HC",
        "H1",
        "H2",
        "H3",
        "HP",
        "H4",
        "H5",
        "HO",
        "H0",
        "HP",
        "O",
        "O2",
        "Cl",
        "Na",
        "I",
        "F",
        "Br",
    ): 1,
    ("NB", "NC", "OW", "OH", "OS", "SH", "S"): 2,
    (
        "C",
        "CN",
        "CB",
        "CR",
        "CK",
        "CC",
        "CW",
        "CV",
        "C*",
        "CQ",
        "CM",
        "CA",
        "CD",
        "CZ",
        "N",
        "NA",
        "N*",
        "N2",
    ): 3,
    ("CT", "N3", "P", "SO"): 4,
}


ATOMTYPE_BONDORDER_FLAT: dict[str, frozenset[int]]
"""
To determin if an atom is a radical.
Maps an AMBER atom type to the numbers of bonds it can have when it is not a
radical. An atom with fewer bonds than the smallest of these is a radical.
Some types allow several numbers: N3 is used for neutral (3 bonds) and
protonated (4 bonds) amines, SH for deprotonated (CYM, 1 bond) and
protonated (CYS, 2 bonds) cysteine and N2 for imine crosslinks (LY5, LX5,
2 bonds) and arginine (3 bonds).
Types with 0 in their set are never detected as radicals, e.g. hydrogens and
ions, because water (settles instead of bonds) and ions have no bonds.
Compare to atom type perception paper (2006) doi:10.1016/j.jmgm.2005.12.005
"""
ATOMTYPE_BONDORDER_FLAT = {
    "MG": frozenset({0}),
    "CO": frozenset({0}),
    "H": frozenset({0, 1}),
    "HW": frozenset({0, 1}),
    "HO": frozenset({0, 1}),
    "HS": frozenset({0, 1}),
    "HA": frozenset({0, 1}),
    "HC": frozenset({0, 1}),
    "H1": frozenset({0, 1}),
    "H2": frozenset({0, 1}),
    "H3": frozenset({0, 1}),
    "HP": frozenset({0, 1}),
    "H4": frozenset({0, 1}),
    "H5": frozenset({0, 1}),
    "H0": frozenset({0, 1}),
    "O": frozenset({0, 1}),
    "O2": frozenset({0, 1}),
    "Cl": frozenset({0, 1}),
    "Na": frozenset({0, 1}),
    "I": frozenset({0, 1}),
    "F": frozenset({0, 1}),
    "Br": frozenset({0, 1}),
    "NB": frozenset({2}),
    "NC": frozenset({2}),
    "OW": frozenset({2}),
    "OH": frozenset({2}),
    "OS": frozenset({2}),
    "SH": frozenset({1, 2}),
    "S": frozenset({2}),
    "C": frozenset({3}),
    "CN": frozenset({3}),
    "CB": frozenset({3}),
    "CR": frozenset({3}),
    "CK": frozenset({3}),
    "CC": frozenset({3}),
    "CW": frozenset({3}),
    "CV": frozenset({3}),
    "C*": frozenset({3}),
    "CQ": frozenset({3}),
    "CM": frozenset({3}),
    "CA": frozenset({3}),
    "CD": frozenset({3}),
    "CZ": frozenset({3}),
    "N": frozenset({3}),
    "NA": frozenset({3}),
    "N*": frozenset({3}),
    "N2": frozenset({2, 3}),
    "N3": frozenset({3, 4}),
    "CT": frozenset({4}),
    "P": frozenset({4}),
    "SO": frozenset({4}),
}

AA3 = [
    "ALA",
    "CYS",
    "ASP",
    "GLU",
    "PHE",
    "GLY",
    "HIE",
    "ILE",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "HYP",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "VAL",
    "TRP",
    "TYR",
    "DOP",
    "ACE",
    "NME",
]
