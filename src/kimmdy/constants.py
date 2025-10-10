"""
Constants used throughout KIMMDY
"""

import sys

from gmx_top4py.constants import (
    ATOM_ID_FIELDS,
    FFFUNC,
    RESNR_ID_FIELDS,
    SOLVENT_NAMES,
    ION_NAMES,
    ATOMTYPE_BONDORDER,
    ATOMTYPE_BONDORDER_FLAT,
    AA3
)

nN_per_kJ_per_mol_nm = 0.001661
R = 8.31446261815324e-3  # [kJ K-1 mol-1]

OPTIONAL_CONFIG_PATHS = ["config.edissoc"]
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

# from gmx data dir
DEFAULT_EDISSOC: dict[tuple[str, str], float] = {  # type: ignore
    tuple(sorted(list(k))): v
    for k, v in {
        ("C", "N"): 500.0,
        ("CA", "C"): 341.0,
        ("CA", "N"): 379.0,
        ("CA", "CB"): 400.0,
        ("CB", "CG"): 400.0,
        ("CG", "CD"): 400.0,
        ("CD", "CE"): 400.0,
        ("CE", "NZ"): 400.0,
    }.items()
}

REACTIVE_MOLECULEYPE = "Reactive"
