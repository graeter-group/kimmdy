"""
Utilities for building plugins, shell convenience functions and GROMACS related functions
"""

from __future__ import annotations

import logging
import re
import subprocess as sp
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import numpy as np

from kimmdy.constants import MARK_REACION_TIME, REACTION_EDR, REACTION_GRO, REACTION_TRR
from kimmdy.recipe import RecipeCollection
from kimmdy.constants import CONFIG_LOGS

if TYPE_CHECKING:
    from kimmdy.parsing import Plumed_dict
    from kimmdy.tasks import TaskFiles
    from kimmdy.topology.topology import Topology

logger = logging.getLogger(__name__)

TopologyAtomAddress = str
"""Address to an atom in the topology.

Corresponds to the 1-based id in the topology and coordinates file.
Note, that gromacs ids in the atomnr column of the gro file
can overflow due to the fixed width file format.
The line number - 2 (for the title and the number of atoms) is always
the correct the atom id.
"""


def flatten_recipe_collections(d: dict[str, RecipeCollection]) -> RecipeCollection:
    return RecipeCollection(recipes=[x for v in d.values() for x in v.recipes])


def field_or_none(l: list[str], i) -> Optional[str]:
    try:
        return l[i]
    except IndexError as _:
        return None


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


### IO utility functions ###


def run_shell_cmd(s, cwd=None) -> sp.CompletedProcess:
    """Run command in shell."""
    return sp.run(s, shell=True, cwd=cwd, capture_output=True, text=True)


def run_gmx(cmd: str, cwd=None) -> Optional[sp.CalledProcessError]:
    """Run GROMACS command in shell.

    Adds a '-quiet' flag to the command and checks the return code.
    """
    logger.debug(f"Starting Gromacs process with command {cmd} in {cwd}.")
    result = run_shell_cmd(f"{cmd} -quiet", cwd)
    if result.returncode != 0:
        logger.error(f"Gromacs process with command {cmd} in {cwd} failed.")
        logger.error(f"Gromacs exit code {result.returncode}.")
        logger.error(f"Gromacs stdout:\n{result.stdout}.")
        logger.error(f"Gromacs stderr:\n{result.stderr}.")
        result.check_returncode()


def get_shell_stdout(s):
    """Run command in shell and capture stdout."""
    process = sp.run(s, shell=True, capture_output=True, encoding="utf-8")
    return process.stdout


def check_file_exists(path: Path, option_name: str | None = None):
    if not path.exists():
        if option_name:
            m = f"File not found: {path} for {option_name}"
        else:
            m = f"File not found: {path}"
        raise LookupError(m)


### reaction plugin building blocks ###


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


def get_edissoc_from_atomnames(
    atomnames: list[str], edissoc: dict, residue: str = "_"
) -> float:
    """Returns dissociation energy E_dissoc for a set of two atomnames.

    Parameters
    ----------
    atomnames:
        Two atomnames as defined in the respective force field
    edissoc:
        Parsed file with dissociation energies per bond between two atomtypes or elements
    residue:
        Residue for which the atomnames are defined
    """
    if residue not in edissoc.keys():
        if "general" in edissoc.keys():
            logger.debug(
                f"residue {residue} not in edissoc keys: {edissoc.keys()}, using 'general' as residue."
            )
            residue = "general"
        else:
            raise KeyError(f"Did not find residue {residue} in edissoc file")

    try:
        interaction_key = tuple(sorted(atomnames))
        E_dis = edissoc[residue][interaction_key]
    except KeyError:
        # continue with guessed edissoc
        logger.warning(
            f"Did not find dissociation energy for atomtypes {atomnames}, residue {residue} in edissoc file, using standard value of 400.0"
        )
        E_dis = 400.0

    return E_dis


def morse_transition_rate(
    r_curr: list[float],
    r_0: float,
    dissociation_energy: float,
    k_f: float,
    frequency_factor: float = 0.288,
    temperature: float = 300,
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
    frequency_factor:
        Prefactor of the Arrhenius equation in [1/ps]. Default value from fitting averaged C_a - N data to gromacs data, see original KIMMDY paper
        Alternatively 1/2pi sqrt(k/m).
    temperature:
        Temperature for the Arrhenius equation in GROMACS units.

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
                (beta**2 * dissociation_energy**2 - 2 * dissociation_energy * beta * fs)
                + 1e-7  # prevent rounding issue close to zero
            )
        )
        / (2 * beta * dissociation_energy)
    )
    r_max = r_0 - 1 / beta * np.log(
        (
            beta * dissociation_energy
            - np.sqrt(
                (beta**2 * dissociation_energy**2 - 2 * dissociation_energy * beta * fs)
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
    R = 8.31446261815324e-3  # [kJ K-1 mol-1]
    k = frequency_factor * np.exp(-delta_v / (R * temperature))  # [1/ps]

    return k, fs


### GROMACS related functions ###


def get_gmx_dir(
    gromacs_alias: str = "gmx", grompp_prefix: Optional[str] = None
) -> Optional[Path]:
    """Returns the path to the gromacs installation"""

    # get the stder from calling `gmx` to search for the `Data prefix:`
    # line which contains the path to the gromacs installation

    # Add prefix if necesarry
    cmd = [gromacs_alias]
    if grompp_prefix:
        cmd.insert(0, grompp_prefix)
    try:
        r = sp.run(cmd, check=False, capture_output=True, text=True)
    except FileNotFoundError:
        logger.warning("GROMACS not found.")
        return None

    from_source = False
    gmx_prefix = None
    for l in r.stderr.splitlines():
        if l.startswith("Data prefix:"):
            gmx_prefix = Path(l.split()[2])
            if "(source tree)" in l:
                from_source = True
            break

    if gmx_prefix is None:
        logger.warning("GROMACS data directory not found in gromacs message.")
        return None

    if from_source:
        return Path(gmx_prefix) / "share"
    else:
        return Path(gmx_prefix) / "share" / "gromacs"


def check_gmx_version(config):
    """Check for an existing gromacs installation.

    If PLUMED is meant to be used it additionally checks for the keyword
    'MODIFIED' or 'plumed' in the version name.

    If slow growth pairs are used, it checks for gromacs version >= 2023.2
    """
    try:
        version = [
            l
            for l in get_shell_stdout(
                f"{config.grompp_prefix + ' ' if config.grompp_prefix else ''}{config.gromacs_alias} --quiet --version"
            ).split("\n")
            if "GROMACS version:" in l
        ][0]
    except Exception as e:
        m = "No system gromacs detected. With error: " + str(e)
        # NOTE: The logger is set up with information from the config
        # the the config can't use the logger.
        # Instead it collects the logmessages and displays them at the end.
        CONFIG_LOGS["errors"].append(m)
        raise SystemError(m)
    # check version for plumed patch if necessary
    if hasattr(config, "mds"):
        for md in config.mds.get_attributes():
            if config.mds.attr(md).use_plumed:
                if not ("MODIFIED" in version or "plumed" in version):
                    m = (
                        "GROMACS version does not contain 'MODIFIED' or "
                        "'plumed', aborting due to apparent lack of PLUMED patch."
                        f"Version is: {version}"
                    )
                    CONFIG_LOGS["errors"].append(m)
                    if not config.dryrun:
                        raise SystemError(m)
    if hasattr(config, "changer") and hasattr(config.changer, "coordinates"):
        if config.changer.coordinates.slow_growth not in ["", "morse_only"]:
            CONFIG_LOGS["debugs"].append(f"Gromacs version: {version}")
            major_minor = re.match(r".*(\d{4})\.(\d+).*", version)
            if major_minor is not None:
                year = int(major_minor.group(1))
                minor = int(major_minor.group(2))
                if year < 2023 or (year == 2023 and minor < 2):
                    m = "Note: slow growth of pairs is only supported by >= gromacs 2023.2. To disable morphing pairs in config: md.changer.coordinates.slow_growth: morse_only"
                    CONFIG_LOGS["errors"].append(m)
                    raise SystemError(m)
    return version


def write_reaction_time_marker(dir: Path, time: float):
    """Write out a file as marker for the reaction time."""
    logger.info(
        f"Writing reaction time marker {time} to {dir.name}/{MARK_REACION_TIME}"
    )
    with open(dir / MARK_REACION_TIME, "w") as f:
        f.write(str(time))


def read_reaction_time_marker(dir: Path) -> float | None:
    if not (dir / MARK_REACION_TIME).exists():
        return None
    with open(dir / MARK_REACION_TIME, "r") as f:
        return float(f.read())


def write_coordinate_files_at_reaction_time(files: TaskFiles, time: float):
    """Write out a gro file from the trajectory (xtc or trr) at the reaction time."""
    gro = files.input["gro"]
    if gro is None:
        m = "No gro file found from the previous md run."
        logger.error(m)
        raise FileNotFoundError(m)

    if gro.name == REACTION_GRO:
        m = f"The latest gro file registered already is a kimmdy reaction file. This state should not be possible unless multiple reactions where run in sequence without any MD in between (even relaxation)."
        logger.error(m)

    gro_reaction = gro.with_name(REACTION_GRO)
    trr_reaction = gro.with_name(REACTION_TRR)
    edr_reaction = gro.with_name(REACTION_EDR)

    if gro_reaction.exists() or trr_reaction.exists() or edr_reaction.exists():
        m = f"gro/trr/edr file at reaction time {time} already exists in {gro.parent.name}. Removing it. This may happen by restarting from a previous run."
        logger.error(m)
        gro_reaction.unlink(True)
        trr_reaction.unlink(True)
        edr_reaction.unlink(True)

    logger.info(
        f"Writing out gro/trr/edr files at reaction time {time} ps in {gro.parent.name}"
    )
    files.output["gro"] = gro_reaction
    files.output["trr"] = trr_reaction
    files.output["edr"] = edr_reaction

    # Prefer xtc over trr
    # It should have more frames and be smaller,
    # but sometimes the people only write a specific index group to the xtc,
    # in which case it fails and we try the trr
    # FIXME: fixme
    # this needs proper documetation for plugin authors
    # because one would want the trr file for the precision and velocities
    # at the raction time, but plugins may use the xtc file (smaller)
    # to determine the reaction and reaction time.
    # this can lead to a mismatch!
    wrote_file = False
    if files.input["xtc"] is not None:
        try:
            run_gmx(
                f"echo '0' | gmx trjconv -f {files.input['xtc']} -s {gro} -b {time} -dump {time} -o {gro_reaction}"
            )
            run_gmx(
                f"echo '0' | gmx trjconv -f {files.input['xtc']} -s {gro} -b {time} -dump {time} -o {trr_reaction}"
            )
            logger.info(
                f"Successfully wrote out gro/trr file {trr_reaction.name} at reaction time in {gro.parent.name} from xtc file."
            )
            wrote_file = True
        except sp.CalledProcessError:
            logger.error(
                f"Failed to write out gro/trr file {trr_reaction.name} at reaction time in {gro.parent.name} from xtc file because the xtc doesn't contain all atoms. Will try trr file."
            )

    if files.input["trr"] is not None:
        try:
            run_gmx(
                f"echo '0' | gmx trjconv -f {files.input['trr']} -s {gro} -b {time} -dump {time} -o {gro_reaction}"
            )
            run_gmx(
                f"echo '0' | gmx trjconv -f {files.input['trr']} -s {gro} -b {time} -dump {time} -o {trr_reaction}"
            )
            logger.info(
                f"Successfully wrote out gro/trr file at reaction time in {gro.parent.name} from trr file."
            )
            wrote_file = True
        except sp.CalledProcessError:
            logger.error(
                f"Failed to write out gro/trr file at reaction time in {gro.parent.name} from trr file."
            )

    if files.input["edr"]:
        try:
            run_gmx(
                f"gmx eneconv -f {files.input['edr']} -b {time} -e {time} -o {edr_reaction}"
            )
            logger.info(
                f"Successfully wrote out edr file at reaction time in {gro.parent.name} from edr file."
            )
        except sp.CalledProcessError:
            logger.error(
                f"Failed to write out edr file at reaction time in {gro.parent.name} from edr file."
            )

    if not wrote_file:
        m = f"No trajectory file found to write out gro/trr file at reaction time in {gro.parent.name}"
        logger.error(m)
        raise FileNotFoundError(m)


def get_task_directories(dir: Path, tasks: Union[list[str], str] = "all") -> list[Path]:
    """
    create list of subdirectories that match the tasks.
    If tasks is "all", all subdirectories are returned.

    Parameters
    ----------
    dir
        Directory to search for subdirectories
    tasks
        List of steps e.g. ["equilibrium", "production"]. Or a string "all" to return all subdirectories
    """
    directories = sorted(
        [
            p
            for p in dir.glob("*_*/")
            if p.is_dir() and "_" in p.name and p.name[0].isdigit()
        ],
        key=lambda p: int(p.name.split("_")[0]),
    )
    if tasks == "all":
        return directories
    else:
        return [d for d in directories if d.name.split("_")[1] in tasks]
