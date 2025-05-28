"""
Utilities for building plugins, shell convenience functions and GROMACS related functions
"""

from __future__ import annotations

import logging
import re
import subprocess as sp
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union


from kimmdy.constants import (
    CONFIG_LOGS,
    MARK_REACION_TIME,
    REACTION_GRO,
)
from kimmdy.recipe import RecipeCollection

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
        if config.changer.coordinates.slow_growth not in ["morse_only", False]:
            CONFIG_LOGS["debugs"].append(f"Gromacs version: {version}")
            major_minor = re.match(r".*(\d{4})\.(\d+).*", version)
            if major_minor is not None:
                year = int(major_minor.group(1))
                minor = int(major_minor.group(2))
                if year < 2023 or (year == 2023 and minor < 2):
                    m = "Note: slow growth of pairs is only supported by >= gromacs 2023.2. To disable morphing pairs in config: md.changer.coordinates.slow_growth: morse_only or remove it entirely to disable slow growth entirely."
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


def write_gro_at_reaction_time(files: TaskFiles, time: float):
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

    if gro_reaction.exists():
        m = f"gro file at reaction time {time} already exists in {gro.parent.name}. Removing it. This may happen by restarting from a previous run."
        logger.error(m)
        gro_reaction.unlink(True)

    logger.info(f"Writing out gro file at reaction time {time} ps in {gro.parent.name}")
    files.output["gro"] = gro_reaction

    wrote_file = False
    if files.input["xtc"] is not None:
        try:
            run_gmx(
                f"echo '0' | gmx trjconv -f {files.input['xtc']} -s {gro} -b {time} -dump {time} -o {gro_reaction}"
            )
            logger.info(
                f"Successfully wrote out gro file {gro_reaction.name} at reaction time in {gro.parent.name} from xtc file."
            )
            wrote_file = True
        except sp.CalledProcessError:
            logger.error(
                f"Failed to write out gro/trr file {gro_reaction.name} at reaction time in {gro.parent.name} from xtc file because the xtc doesn't contain all atoms. Will try trr file."
            )

    if files.input["trr"] is not None:
        try:
            run_gmx(
                f"echo '0' | gmx trjconv -f {files.input['trr']} -s {gro} -b {time} -dump {time} -o {gro_reaction}"
            )
            logger.info(
                f"Successfully wrote out gro file at reaction time in {gro.parent.name} from trr file."
            )
            wrote_file = True
        except sp.CalledProcessError:
            logger.error(
                f"Failed to write out gro file at reaction time in {gro.parent.name} from trr file."
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
