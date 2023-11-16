"""
Utilities for building plugins, shell convenience functions and GROMACS related functions
"""
from __future__ import annotations
import subprocess as sp
import numpy as np
import re
import logging
from typing import Optional, TYPE_CHECKING
from pathlib import Path

if TYPE_CHECKING:
    from kimmdy.tasks import TaskFiles
    from kimmdy.topology.topology import Topology
    from kimmdy.parsing import Plumed_dict

logger = logging.getLogger(__name__)

TopologyAtomAddress = str
"""Address to an atom in the topology.

Corresponds to the 1-based id in the topology and coordinates file.
Note, that gromacs ids in the atomnr column of the gro file
can overflow due to the fixed width file format.
The line number - 2 (for the title and the number of atoms) is always
the correct the atom id.
"""


class longFormatter(logging.Formatter):
    def format(self, record):
        saved_name = record.name  # save and restore for other formatters if desired
        parts = saved_name.split(".")
        if len(parts) > 1:
            record.name = parts[0][0] + "." + ".".join(p[:10] for p in parts[1:])
        else:
            record.name = parts[0]
        result = super().format(record)
        record.name = saved_name
        return result


# input/output utility functions


def run_shell_cmd(s, cwd=None) -> sp.CompletedProcess:
    """Run command in shell."""
    return sp.run(s, shell=True, cwd=cwd, capture_output=True, text=True)


def run_gmx(s: str, cwd=None) -> Optional[sp.CalledProcessError]:
    """Run GROMACS command in shell.

    Adds a '-quiet' flag to the command and checks the return code.
    """
    result = run_shell_cmd(f"{s} -quiet", cwd)
    if result.returncode != 0:
        logger.error(f"Gromacs process with command {s} in {cwd} failed.")
        logger.error(f"Gromacs exit code {result.returncode}.")
        logger.error(f"Gromacs stdout:\n{result.stdout}.")
        logger.error(f"Gromacs stderr:\n{result.stderr}.")
        result.check_returncode()


def get_shell_stdout(s):
    """Run command in shell and capture stdout."""
    process = sp.run(s, shell=True, capture_output=True, encoding="utf-8")
    return process.stdout


def check_file_exists(p: Path):
    if not p.exists():
        m = f"File not found: {p}"
        raise LookupError(m)


# reaction plugin building blocks


def get_atomnrs_from_plumedid(
    plumedid: str,
    plumed: Plumed_dict,
) -> list[str]:
    """
    Convert from plumedid to atomnr, information from the plumed file is used.

    Parameters
    ----------
    plumedid:
        Identifier from a plumed input file (e.g d0).
    plumed:
        Parsed plumed input file
    """
    # lookup_atomnr_plumedid = {k: frozenset(v["atoms"])
    plumed_action = plumed["labeled_action"][plumedid]
    if a := plumed_action.get("atoms"):
        atomnrs = sorted(a, key=int)
        return atomnrs
    else:
        raise NotImplementedError(
            f"Can't get atomnrs for {plumedid}, is this for an unexpected plumed action?"
        )


def get_atominfo_from_atomnrs(
    atomnrs: list[str], top: Topology
) -> tuple[list[str], list[str]]:
    """Use topology atoms section to convert from atomnr to atomtype"""
    atomtypes = []
    atomnames = []
    for atomnr in atomnrs:
        atomtypes.append(top.atoms[atomnr].type)
        atomnames.append(top.atoms[atomnr].atom)
    return atomtypes, atomnames


def get_bondprm_from_atomtypes(
    atomtypes: list[str],
    ffbonded: dict,
) -> tuple[float, float]:
    """Returns bond parameters (b0, kb) for a set of two atomtypes.

    Parameters
    ----------
    atomtypes:
        Two atomtypes as defined in the respective force field
    ffbonded:
        Force field ffbonded.itp file parsed through the rtp parser
    """
    # search for b0 and kb for the given atomtypes in ffbonded bondtypes
    for bondtype in ffbonded["bondtypes"]["content"]:
        if set(atomtypes) == set(bondtype[:2]):
            b0, kb = [float(x) for x in bondtype[3:5]]
            break
    else:
        raise KeyError(
            f"Did not find bond parameters for atomtypes {atomtypes} in ffbonded file"
        )

    return b0, kb


def get_edissoc_from_atomnames(atomnames: list[str], edissoc: dict) -> float:
    """Returns dissociation energy E_dissoc for a set of two atomnames.

    Parameters
    ----------
    atomnames:
        Two atomnames as defined in the respective force field
    edissoc:
        Parsed file with dissociation energies per bond between two atomtypes or elements
    """
    # dissociation energy can be for bonds between atomtypes or elements or mixtures of both
    atomelements = [x[0] for x in atomnames]
    for comb in [
        atomnames,
        [atomnames[0], atomelements[1]],
        [atomelements[0], atomnames[1]],
        atomelements,
    ]:
        if E_dis := edissoc.get(tuple(comb)):
            break
        elif E_dis := edissoc.get(tuple(comb[::-1])):
            break
    else:
        raise KeyError(
            f"Did not find dissociation energy for atomtypes {atomnames} in edissoc file"
        )
    return E_dis


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
                (
                    beta**2 * dissociation_energy**2
                    - 2 * dissociation_energy * beta * fs
                )
                + 1e-7  # prevent rounding issue close to zero
            )
        )
        / (2 * beta * dissociation_energy)
    )
    r_max = r_0 - 1 / beta * np.log(
        (
            beta * dissociation_energy
            - np.sqrt(
                (
                    beta**2 * dissociation_energy**2
                    - 2 * dissociation_energy * beta * fs
                )
                + 1e-7  # prevent rounding issue close to zero
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


# GROMACS related functions


def get_gmx_dir(gromacs_alias: str = "gmx") -> Optional[Path]:
    """Returns the path to the gromacs installation"""

    # get the stder from calling `gmx` to search for the `Data prefix:`
    # line which contains the path to the gromacs installation
    try:
        r = sp.run([gromacs_alias], check=False, capture_output=True, text=True)
    except FileNotFoundError:
        logger.warning("GROMACS not found.")
        return None

    gmx_prefix = None
    for l in r.stderr.splitlines():
        if l.startswith("Data prefix:"):
            gmx_prefix = Path(l.split()[2])
            break

    if gmx_prefix is None:
        logger.warning("GROMACS data directory not found in gromacs message.")
        return None

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


def truncate_sim_files(files: TaskFiles, time: Optional[float], keep_tail: bool = True):
    """Truncates latest trr, xtc, edr, and gro to the time to a previous
    point in time.

    The files stay in place, the truncated tail is by default kept and renamed
    to '[...xtc].tail'

    Parameters
    ----------
    time
        Time in ps up to which the data should be truncated.
    files
        TaskFiles to get the latest files.
    """

    if time is None:
        logger.debug("time is None, nothing to truncate")
        return

    paths = {}
    paths["gro"] = files.input["gro"]
    for s in ["trr", "xtc", "edr"]:
        try:
            paths[s] = files.input[s]
        except FileNotFoundError:
            paths[s] = None

    # trr or xtc must be present
    if (traj := paths["trr"]) is None:
        if (traj := paths["xtc"]) is None:
            raise RuntimeError("No trajectory file!")

    # check time exists in traj
    p = sp.run(
        f"gmx -quiet -nocopyright check -f {traj}",
        text=True,
        capture_output=True,
        shell=True,
    )
    # FOR SOME REASON gmx check writes in stderr instead of stdout
    if m := re.search(r"Last frame.*time\s+(\d+\.\d+)", p.stderr):
        last_time = float(m.group(1))
        assert (
            last_time * 1.01 >= time
        ), "Requested to truncate trajectory after last frame"
    else:
        raise RuntimeError(f"gmx check failed:\n{p.stdout}\n{p.stderr}")
    logger.info(
        f"Truncating trajectories to {time:.4} ps. Trajectory time was {last_time:.4} ps"
    )

    # backup the tails of trajectories
    for trj in [paths["trr"], paths["xtc"]]:
        if trj is None:
            continue
        tmp = trj.rename(trj.with_name("tmp_backup_" + trj.name))
        if keep_tail:
            run_gmx(
                f"gmx trjconv -f {tmp} -b {time} -o {trj}",
            )
            trj.rename(str(trj) + ".tail")

        run_gmx(f"gmx trjconv -f {tmp} -e {time} -o {trj}")
        tmp.unlink()

    # backup the gro
    bck_gro = paths["gro"].rename(
        paths["gro"].with_name("tmp_backup_" + paths["gro"].name)
    )
    sp.run(
        f"gmx trjconv -f {traj} -s {bck_gro} -dump -1 -o {paths['gro']}",
        text=True,
        input="0",
        shell=True,
    )
    bck_gro.rename(str(paths["gro"]) + ".tail")
    if not keep_tail:
        bck_gro.unlink()

    # backup the edr
    if paths["edr"] is not None:
        bck_edr = paths["edr"].rename(
            paths["edr"].with_name("tmp_backup_" + paths["edr"].name)
        )
        run_shell_cmd(f"gmx eneconv -f {bck_edr} -e {time} -o {paths['edr']}")
        bck_edr.rename(str(paths["edr"]) + ".tail")
        if not keep_tail:
            bck_edr.unlink()
    return
