import subprocess as sp
from kimmdy.topology.utils import get_protein_section
import numpy as np
import logging
from typing import Optional
from pathlib import Path



def get_gmx_dir() -> Path:
    """returns the path to the gromacs installation"""
    gmx_binary = Path(
        sp.run(["which", "gmx"], capture_output=True).stdout.decode().strip()
    )
    if not gmx_binary.exists():
        raise ValueError("Could not find gromacs installation")
    # resolve symlink if gmxbinary is a symlink
    gmx_binary = gmx_binary.resolve()
    gmx_dir = gmx_binary.parent.parent / "share" / "gromacs"
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
) -> tuple[frozenset, frozenset]:
    """returns atomtypes for a plumedid with information from the plumed and topology file"""
    lookup_atomid_plumedid = {
        entry["id"]: frozenset(entry["atoms"]) for entry in plumed["distances"]
    }
    atoms = get_protein_section(top, "atoms")
    if not atoms:
        raise ValueError("Could not find atoms in topology file")
    lookup_atomtype_atomid = {int(atom[0]): atom[1] for atom in atoms}
    atomids = lookup_atomid_plumedid[plumedid]
    atomids_list = list(atomids)
    atomtypes_list = [
        lookup_atomtype_atomid[atomids_list[0]],
        lookup_atomtype_atomid[atomids_list[1]],
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
    E_dis: float,
    k_f: float,
    k_0: float = 0.288,
    kT: float = 2.479,
) -> tuple[list[float], list[float]]:
    """calculates energy barrier crossing rate [in ps]; barrier based on the model V = V_morse - F*X"""
    r_curr = np.asarray(r_curr)
    beta = np.sqrt(k_f / (2 * E_dis))

    Fs = (
        2
        * beta
        * E_dis
        * np.exp(-beta * (r_curr - r_0))
        * (1 - np.exp(-beta * (r_curr - r_0)))
    )

    # inflection calculation
    r_infl = (beta * r_0 + np.log(2)) / beta
    F_infl = (
        2
        * beta
        * E_dis
        * np.exp(-beta * (r_infl - r_0))
        * (1 - np.exp(-beta * (r_infl - r_0)))
    )
    Fs_mask = r_curr > r_infl
    Fs[Fs_mask] = F_infl

    # calculate extrema of shifted potential i.o.t. get barrier hight
    rmin = r_0 - 1 / beta * np.log(
        (beta * E_dis + np.sqrt(beta**2 * E_dis**2 - 2 * E_dis * beta * Fs))
        / (2 * beta * E_dis)
    )
    rmax = r_0 - 1 / beta * np.log(
        (beta * E_dis - np.sqrt(beta**2 * E_dis**2 - 2 * E_dis * beta * Fs))
        / (2 * beta * E_dis)
    )
    # set rmax to r0 * 10 where no rmax can be found
    rmax = np.where(~np.isfinite(rmax), 10 * r_0, rmax)
    Vmax = E_dis * (1 - np.exp(-beta * (rmax - r_0))) ** 2 - Fs * (rmax - r_0)
    Vmin = E_dis * (1 - np.exp(-beta * (rmin - r_0))) ** 2 - Fs * (rmin - r_0)
    # Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero

    delta_V = Vmax - Vmin
    k = k_0 * np.exp(-delta_V / kT)  # [1/ps]

    return k, Fs


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
    'MODIFIED' in the version name.
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

    # i hate this
    if (
        any(
            "plumed" in y
            for y in [
                config.mds.attr(x).get_attributes() for x in config.mds.get_attributes()
            ]
        )
        and not "MODIFIED" in version
    ):
        m = "GROMACS version does not contain MODIFIED, aborting due to lack of PLUMED patch."
        logging.error(m)
        logging.error("Version was: " + version)
        if not config.dryrun:
            raise SystemError(m)
    return version


