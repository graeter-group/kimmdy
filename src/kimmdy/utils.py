import subprocess as sp
from kimmdy.topology.utils import get_protein_section
import numpy as np
import logging
from typing import Optional
from pathlib import Path


def get_gmx_dir() -> Path:
    """returns the path to the gromacs installation"""

    # get the stder from calling `gmx` to search for the `Data prefix:`
    # line which contains the path to the gromacs installation
    r = sp.run(["gmx"], check=True, capture_output=True, text=True)

    gmx_prefix = None
    for l in r.stderr.splitlines():
        if l.startswith("Data prefix:"):
            gmx_prefix = l.split()[2]
            break

    if gmx_prefix is None:
        raise ValueError("Could not find gromacs installation")

    gmx_dir = Path(gmx_prefix) / "share" / "gromacs"
    return gmx_dir


def increment_logfile(f: Path) -> Path:
    backup_file_prefix = "#"
    backup_file_suffix = "#"
    logfile = f
    if logfile.exists():
        backup_count = 1
        backup_file = (
            f"{backup_file_prefix}{logfile}_{backup_count}{backup_file_suffix}"
        )
        while Path(backup_file).exists():
            backup_count += 1
            backup_file = (
                f"{backup_file_prefix}{logfile}_{backup_count}{backup_file_suffix}"
            )
        logfile.rename(backup_file)
    return logfile


def get_atominfo_from_plumedid(
    plumedid: str, plumed: dict, top: dict
) -> tuple[frozenset[str], list[str]]:
    """returns atomtypes for a plumedid with information from the plumed and topology file"""
    lookup_atomid_plumedid = {
        entry["id"]: frozenset(entry["atoms"]) for entry in plumed["distances"]
    }
    atoms = get_protein_section(top, "atoms")
    if not atoms:
        raise ValueError("Could not find atoms in topology file")
    lookup_atomtype_atomid = {str(atom[0]): atom[1] for atom in atoms}
    atomids = sorted(lookup_atomid_plumedid[plumedid], key=int)
    atomtypes_list = [
        lookup_atomtype_atomid[atomids[0]],
        lookup_atomtype_atomid[atomids[1]],
    ]
    atomtypes = frozenset(atomtypes_list)
    logging.debug(f"Found atomtypes {atomtypes} for plumedid {plumedid}.")
    return atomtypes, atomids


def get_bondprm_from_atomtypes(
    atomtypes: frozenset, ffbonded: dict, lookup_edissoc_atomtype: dict
) -> tuple[float, float, float]:
    """returns bond parameters (b0, kb, E_dis) for a set of atomtypes"""
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

    logging.debug(f"Found bondprm {[b0,kb,E_dis]} for atomtypes {atomtypes}.")
    return float(b0), float(kb), float(E_dis)


def morse_transition_rate(
    r_curr: list[float],
    r_0: float,
    dissociation_energies: float,
    k_f: float,
    k_0: float = 0.288,
    kT: float = 2.479,
) -> tuple[list[float], list[float]]:
    """calculates energy barrier crossing rate [in 1/ps]; barrier based on the model V = V_morse - F*X"""
    rs = np.asarray(r_curr)
    beta = np.sqrt(k_f / (2 * dissociation_energies))

    fs = (
        2
        * beta
        * dissociation_energies
        * np.exp(-beta * (rs - r_0))
        * (1 - np.exp(-beta * (rs - r_0)))
    )

    # inflection calculation
    r_inflection = (beta * r_0 + np.log(2)) / beta
    f_inflection = (
        2
        * beta
        * dissociation_energies
        * np.exp(-beta * (r_inflection - r_0))
        * (1 - np.exp(-beta * (r_inflection - r_0)))
    )
    fs_mask = rs > r_inflection
    fs[fs_mask] = f_inflection

    # calculate extrema of shifted potential i.o.t. get barrier hight
    r_min = r_0 - 1 / beta * np.log(
        (
            beta * dissociation_energies
            + np.sqrt(
                beta**2 * dissociation_energies**2
                - 2 * dissociation_energies * beta * fs
            )
        )
        / (2 * beta * dissociation_energies)
    )
    r_max = r_0 - 1 / beta * np.log(
        (
            beta * dissociation_energies
            - np.sqrt(
                beta**2 * dissociation_energies**2
                - 2 * dissociation_energies * beta * fs
            )
        )
        / (2 * beta * dissociation_energies)
    )
    logging.error(r_max)
    # set rmax to r0 * 10 where no rmax can be found
    r_max = np.where(~np.isfinite(r_max), 10 * r_0, r_max)
    logging.error(r_max)
    v_max = dissociation_energies * (1 - np.exp(-beta * (r_max - r_0))) ** 2 - fs * (
        r_max - r_0
    )
    v_min = dissociation_energies * (1 - np.exp(-beta * (r_min - r_0))) ** 2 - fs * (
        r_min - r_0
    )
    # Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero

    delta_v = v_max - v_min
    k = k_0 * np.exp(-delta_v / kT)  # [1/ps]

    return k, fs


def run_shell_cmd(s, cwd=None) -> sp.CompletedProcess:
    return sp.run(s, shell=True, cwd=cwd)


def run_gmx(s: str, cwd=None) -> Optional[sp.CalledProcessError]:
    result = run_shell_cmd(f"{s} -quiet", cwd)
    if result.returncode != 0:
        logging.error(f"Gromacs process failed with exit code {result.returncode}.")
        result.check_returncode()


def get_shell_stdout(s):
    process = sp.run(s, shell=True, capture_output=True, encoding="utf-8")
    return process.stdout


def check_gmx_version(config):
    """Check for an existing gromacs installation.

    If PLUMED is meant to be used it additionally checks for the keyword
    'MODIFIED' or 'plumed' in the version name.
    """
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
        logging.error(m)
        raise SystemError(m)

    if any(
        "plumed" in y
        for y in [
            config.mds.attr(x).get_attributes() for x in config.mds.get_attributes()
        ]
    ) and (not ("MODIFIED" in version or "plumed" in version)):
        m = "GROMACS version does not contain MODIFIED or plumed, aborting due to lack of PLUMED patch."
        logging.error(m)
        logging.error("Version was: " + version)
        if not config.dryrun:
            raise SystemError(m)
    return version
