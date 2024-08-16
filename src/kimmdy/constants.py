"""
Constants used throughout KIMMDY
"""

MARK_STARTED = ".kimmdy_started"
MARK_DONE = ".kimmdy_done"
MARK_FAILED = ".kimmdy_failed"
MARKERS = [MARK_STARTED, MARK_DONE, MARK_FAILED]

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

# see https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html
FFFUNC = {
    "mult_proper_dihedral": "9",
    "mult_improper_dihedral": "4",
    "harmonic_bond": "1",
    "morse_bond": "3",
    "harmonic_angle": "1",
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


ATOMTYPE_BONDORDER_FLAT: dict[str, int]
"""
To determin if an atom is a radical.
Compare to atom type perception paper (2006) doi:10.1016/j.jmgm.2005.12.005
"""
ATOMTYPE_BONDORDER_FLAT = {
    "MG": 0,
    "NA": 0,
    "CO": 0,
    "H": 0,
    "HW": 0,
    "HO": 0,
    "HS": 0,
    "HA": 0,
    "HC": 0,
    "H1": 0,
    "H2": 0,
    "H3": 0,
    "HP": 0,
    "H4": 0,
    "H5": 0,
    "HO": 0,
    "H0": 0,
    "HP": 0,
    "O": 0,
    "O2": 0,
    "Cl": 0,
    "Na": 0,
    "I": 0,
    "F": 0,
    "Br": 0,
    "NB": 2,
    "NC": 2,
    "OW": 2,
    "OH": 2,
    "OS": 2,
    "SH": 2,
    "S": 2,
    "C": 3,
    "CN": 3,
    "CB": 3,
    "CR": 3,
    "CK": 3,
    "CC": 3,
    "CW": 3,
    "CV": 3,
    "C*": 3,
    "CQ": 3,
    "CM": 3,
    "CA": 3,
    "CD": 3,
    "CZ": 3,
    "N": 3,
    "NA": 3,
    "N*": 3,
    "N2": 3,
    "CT": 4,
    "N3": 4,
    "P": 4,
    "SO": 4,
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
