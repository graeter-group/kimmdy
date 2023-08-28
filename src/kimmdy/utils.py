"""
Utilities for building plugins, shell convenience functions and GROMACS related functions
"""

import subprocess as sp
import numpy as np
import logging
from typing import Optional
from pathlib import Path

from kimmdy.topology.utils import get_protein_section

logger = logging.getLogger(__name__)

TopologyAtomAddress = str | tuple[str, str] | tuple[int, str]
"""Address to an atom in the topology.

One of (id: str), (moleculetype: str, id: str) or (moleculetype_ix: int, id).
"""


## input/output utility functions


def run_shell_cmd(s, cwd=None) -> sp.CompletedProcess:
    """Run command in shell."""
    return sp.run(s, shell=True, cwd=cwd)


def run_gmx(s: str, cwd=None) -> Optional[sp.CalledProcessError]:
    """Run GROMACS command in shell.

    Adds a '-quiet' flag to the command and checks the return code.
    """
    result = run_shell_cmd(f"{s} -quiet", cwd)
    if result.returncode != 0:
        logger.error(f"Gromacs process failed with exit code {result.returncode}.")
        result.check_returncode()


def get_shell_stdout(s):
    """Run command in shell and capture stdout."""
    process = sp.run(s, shell=True, capture_output=True, encoding="utf-8")
    return process.stdout


def backup_if_existing(f: Path) -> None:
    """Checks whether a file exists and if so, backs it up.

    Mainly used for files from previous KIMMDY runs to
    prevent overwriting these files.

    Parameters
    ----------
    f:
        Path to a file that will be backed up, if existing.
    """
    backup_file_prefix = "#"
    backup_file_suffix = "#"
    if f.exists():
        backup_count = 1
        backup_file = f"{backup_file_prefix}{f}_{backup_count}{backup_file_suffix}"
        while Path(backup_file).exists():
            backup_count += 1
            backup_file = f"{backup_file_prefix}{f}_{backup_count}{backup_file_suffix}"
        f.rename(backup_file)


## reaction plugin building blocks


def get_atominfo_from_plumedid(
    plumedid: str, plumed: dict, top: dict
) -> tuple[frozenset[str], list[str]]:
    """
    For a plumedid, returns the corresponding atom types and nrs.

    To convert from plumedid to atomnr, information from the plumed file is used.
    Then, the topology atoms section can be used to convert from atomnr to atomtype

    Parameters
    ----------
    plumedid:
        Identifier from a plumed input file (e.g d0).
    plumed:
        Parsed plumed input file
    top:
        Topology of the molecular system"""

    lookup_atomnr_plumedid = {
        entry["id"]: frozenset(entry["atoms"]) for entry in plumed["distances"]
    }
    atoms = get_protein_section(top, "atoms")
    lookup_atomtype_atomnr = {str(atom[0]): atom[1] for atom in atoms}
    atomnrs = sorted(lookup_atomnr_plumedid[plumedid], key=int)
    atomtypes = frozenset(
        [
            lookup_atomtype_atomnr[atomnrs[0]],
            lookup_atomtype_atomnr[atomnrs[1]],
        ]
    )
    logger.debug(f"Found atomtypes {atomtypes} for plumedid {plumedid}.")
    return atomtypes, atomnrs


def get_bondprm_from_atomtypes(
    atomtypes: frozenset, ffbonded: dict, lookup_edissoc_atomtype: dict
) -> tuple[float, float, float]:
    """Returns bond parameters (b0, kb, E_dis) for a set of two atomtypes.

    Parameters
    ----------
    atomtypes:
        Two atomtypes as defined in the respective force field
    ffbonded:
        Force field ffbonded.itp file parsed through the rtp parser
    lookup_edissoc_atomtype:
        Parsed file with dissociation energies per bond between two atomtypes or elements
    """
    atomtypes_list = list(atomtypes)
    lookup_ffbonded_atomtype = {
        frozenset(l[:2]): [float(l[3]), float(l[4])]
        for l in ffbonded["bondtypes"]["other"]
    }
    atomelements_list = [x[0] for x in atomtypes_list]

    # dissociation energy can be for bonds between atomtypes or elements or mixtures of both
    for comb in [
        atomtypes_list,
        [atomtypes_list[0], atomelements_list[1]],
        [atomelements_list[0], atomtypes_list[1]],
        atomelements_list,
    ]:
        if (comb_set := frozenset(comb)) in lookup_edissoc_atomtype.keys():
            E_dis = lookup_edissoc_atomtype[comb_set]
            break
    else:
        raise KeyError(
            f"Did not find dissociation energy for atomtypes {atomtypes} in edissoc file"
        )

    try:
        b0, kb = lookup_ffbonded_atomtype[atomtypes]
    except KeyError as e:
        raise KeyError(
            f"Did not find bond parameters for atomtypes {atomtypes} in ffbonded file"
        ) from e

    logger.debug(f"Found bondprm {[b0,kb,E_dis]} for atomtypes {atomtypes}.")
    return float(b0), float(kb), float(E_dis)


def morse_transition_rate(
    r_curr: list[float],
    r_0: float,
    dissociation_energy: float,
    k_f: float,
    k_0: float = 0.288,
    kT: float = 2.479,
) -> tuple[list[float], list[float]]:
    """Calculates reaction rate constant for a bond breaking event.

    Uses the Morse potential model for this calculation. For an array of bond distances of the same bond,
    first calculates the forces on the bond, then the minima and maxima of the shifted Morse potential
    to get an energy barrier and finally a reaction rate constant using the Arrhenius equation.
    For intramolecular reactions, the reaction rate constant is equal to the reaction rate.

    The calculation should be according to the derivation in the original KIMMDY paper: DOI: 10.1021/acs.jctc.9b00786

    Parameters
    ----------
    r_curr:
        Bond distances for a single bond, typically from a time series.
    r_0:
        Equilibrium bond length of the bond.
    dissociation energy:
        Dissociation energy of the bond.
    k_f:
        Spring constant of the bond.
    k_0:
        Prefactor of the Arrhenius equation in [1/ps]. Default value from fitting averaged C_a - N data to gromacs data, see original KIMMDY paper
        Alternatively 1/2pi sqrt(k/m).
    kT:
        Constant in the Arrhenius equation in GROMACS units [kJ mol-1], default for 310K.

    """
    rs = np.asarray(r_curr)
    beta = np.sqrt(k_f / (2 * dissociation_energy))

    # calculate forces on bond
    fs = (
        2
        * beta
        * dissociation_energy
        * np.exp(-beta * (rs - r_0))
        * (1 - np.exp(-beta * (rs - r_0)))
    )

    # if the bond is stretched beyond the inflection point, take the inflection point force because this force must have acted on the bond at some point
    r_inflection = (beta * r_0 + np.log(2)) / beta
    f_inflection = (
        2
        * beta
        * dissociation_energy
        * np.exp(-beta * (r_inflection - r_0))
        * (1 - np.exp(-beta * (r_inflection - r_0)))
    )
    fs_mask = rs > r_inflection
    fs[fs_mask] = f_inflection

    # calculate extrema of shifted potential i.e. get barrier height of V_eff = V_morse - F*X
    r_min = r_0 - 1 / beta * np.log(
        (
            beta * dissociation_energy
            + np.sqrt(
                beta**2 * dissociation_energy**2
                - 2 * dissociation_energy * beta * fs
            )
        )
        / (2 * beta * dissociation_energy)
    )
    r_max = r_0 - 1 / beta * np.log(
        (
            beta * dissociation_energy
            - np.sqrt(
                beta**2 * dissociation_energy**2
                - 2 * dissociation_energy * beta * fs
            )
        )
        / (2 * beta * dissociation_energy)
    )
    r_max = np.where(
        ~np.isfinite(r_max), 10 * r_0, r_max
    )  # set rmax to r0 * 10 where no rmax can be found

    v_max = dissociation_energy * (1 - np.exp(-beta * (r_max - r_0))) ** 2 - fs * (
        r_max - r_0
    )
    v_min = dissociation_energy * (1 - np.exp(-beta * (r_min - r_0))) ** 2 - fs * (
        r_min - r_0
    )
    # Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero
    delta_v = v_max - v_min

    # calculate reaction rate constant from barrier heigth
    k = k_0 * np.exp(-delta_v / kT)  # [1/ps]

    return k, fs


## GROMACS related functions


def get_gmx_dir() -> Path:
    """Returns the path to the gromacs installation

    This does not check if the installation is valid.
    It just returns the path to the gromacs data directory.
    If `gmx` is not executable it still returns the default
    gromacs data directory in `/usr/share/gromacs`.
    """

    # get the stder from calling `gmx` to search for the `Data prefix:`
    # line which contains the path to the gromacs installation
    r = sp.run(["gmx"], check=False, capture_output=True, text=True)

    gmx_prefix = "/usr"
    for l in r.stderr.splitlines():
        if l.startswith("Data prefix:"):
            gmx_prefix = l.split()[2]
            break

    gmx_dir = Path(gmx_prefix) / "share" / "gromacs"
    return gmx_dir


def check_gmx_version(config):
    """Check for an existing gromacs installation.

    If PLUMED is meant to be used it additionally checks for the keyword
    'MODIFIED' or 'plumed' in the version name.
    """
    # check for existing installation and get version
    try:
        version = [
            l
            for l in get_shell_stdout(
                f"{config.gromacs_alias} --quiet --version"
            ).split("\n")
            if "GROMACS version:" in l
        ][0]
    except Exception as e:
        m = "No system gromacs detected. With error: " + str(e)
        logger.error(m)
        raise SystemError(m)
    # check version for plumed patch if necessary
    if hasattr(config, "mds"):
        for md in config.mds.get_attributes():
            if config.mds.attr(md).use_plumed:
                if not ("MODIFIED" in version or "plumed" in version):
                    m = (
                        "GROMACS version does not contain 'MODIFIED' or "
                        "'plumed', aborting due to apparent lack of PLUMED patch."
                    )
                    logger.error(m)
                    logger.error("Version was: " + version)
                    if not config.dryrun:
                        raise SystemError(m)
    return version
