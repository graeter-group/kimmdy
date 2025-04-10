"""
The tasks module holds the TaskFiles class which organizes input and
output paths and the Task class for tasks in the runmanager queue.
"""

import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

from kimmdy.constants import MARK_DONE, MARK_STARTED
from kimmdy.parsing import read_plumed, write_time_marker
from kimmdy.utils import longFormatter

logger = logging.getLogger(__name__)


class AutoFillDict(dict):
    """Dictionary that gets populated by calling get_missing."""

    def __init__(self, get_missing: Callable[[str], Path | None]):
        self.get_missing = get_missing

    def __missing__(self, key: str) -> None | Path:
        v = self.get_missing(key)
        if v is not None:
            self[key] = v
            return v
        else:
            return None


@dataclass
class TaskFiles:
    """Class for Task input and output files and directories.

    Hosts the input and output file paths belonging to a task.
    A function or method that wants to be callable as a Task
    has to return a TaskFiles object.
    The input defaultdict is populated on the fly using
    get_latest of the runmanager to find newest files.
    Files which can not be found by get_latest must be added manually.

    Attributes
    ----------
    get_latest:
        Runmanager.get_latest function that returns paths to the latest file of given type.
    input:
        Input file paths for a Task. Is populated by get_latest or manually.
    output:
        Output file paths for a Task. Is populated by runmanager._discover_output_files or manually.
    outputdir:
        Output directory for a Task. Typically populated by create_task_directory called by Task.
    logger:
        Logger for a Task. Initialized in create_task_directory.

    Examples
    --------
    >>> class run():
    >>>     def get_latest(self, s):
    >>>         return f"latest {s}"
    >>> runmng = run()
    >>> files = TaskFiles(runmng.get_latest)
    >>> files.input
    >>> files.input["top"]
    {'top': 'latest top'}
    """

    get_latest: Callable[[str], Path | None]
    input: dict[str, Path | None] = field(default_factory=dict)
    output: dict[str, Path] = field(default_factory=dict)
    outputdir: Path = Path()
    logger: logging.Logger = logging.getLogger("kimmdy.basetask")

    def __post_init__(self):
        self.input = AutoFillDict(self.get_latest)


def create_task_directory(
    runmng, postfix: str, is_continuation: bool = False
) -> TaskFiles:
    """Creates TaskFiles object, output directory, logger and symlinks ff.

    Gets called when a Task is called (from the runmanager.tasks queue).
    """

    files = TaskFiles(runmng.get_latest)
    runmng.iteration += 1
    taskname = f"{runmng.iteration}_{postfix}"

    # create outputdir
    files.outputdir = runmng.config.out / taskname
    logger.debug(f"Creating Output directory: {files.outputdir}")
    if files.outputdir.exists() and not is_continuation:
        logger.warning(
            f"Output directory {files.outputdir} for the task already exists. Deleting."
        )
        shutil.rmtree(files.outputdir)
    if not is_continuation:
        files.outputdir.mkdir()

    # set up logger
    files.logger = logging.getLogger(f"kimmdy.{taskname}")
    hand = logging.FileHandler(files.outputdir / (taskname + ".log"))
    hand.setFormatter(
        longFormatter(
            "%(asctime)s %(name)-12s %(levelname)s: %(message)s", "%d-%m-%y %H:%M:%S"
        )
    )
    files.logger.addHandler(hand)
    # symlink force field
    if runmng.config.ff is not None:
        if not (files.outputdir / runmng.config.ff.name).exists():
            (files.outputdir / runmng.config.ff.name).symlink_to(runmng.config.ff)
        if not (files.outputdir / "residuetypes.dat").exists():
            (files.outputdir / "residuetypes.dat").symlink_to(
                runmng.config.ff / "residuetypes.dat"
            )

    return files


class Task:
    """A task to be performed as part of the RunManager Queue.

    A task consists of a function and its keyword arguments.
    Calling a taks calls the stored function.
    The function must return a TaskFiles object.

    Parameters
    ----------
    runmng
        Runmanager instance from which the task is called
    f
        Function that will be called when the task is called
    kwargs
        kwargs to be passed to f
    out
        If not None, an output dir will be created with this name
    """

    def __init__(
        self,
        runmng,
        f: Callable[..., Optional[TaskFiles]],
        kwargs: Optional[dict[str, Any]] = None,
        out: Optional[str] = None,
    ):
        self.runmng = runmng
        self.f = f
        if kwargs is None:
            kwargs = {}
        self.kwargs = kwargs
        self.name = self.f.__name__.strip("_")
        self.out = out

        logger.debug(f"Init task {self.name}\tkwargs: {self.kwargs}\tOut: {self.out}")

    def __call__(self) -> Optional[TaskFiles]:
        logger.debug(
            f"Starting task: {self.name} with args: {self.kwargs} in {self.runmng.iteration+1}_{self.out}"
        )
        logger.info(
            f"Starting task: {self.name} in {self.runmng.iteration+1}_{self.out}"
        )
        if self.out is not None:
            is_continuation = False
            if self.kwargs.get("continue_md") is True:
                logger.info(
                    f"Continuing task: {self.name} in {self.runmng.iteration}_{self.out}"
                )
                is_continuation = True
            self.kwargs.update(
                {
                    "files": create_task_directory(
                        self.runmng, self.out, is_continuation=is_continuation
                    )
                }
            )
            write_time_marker(self.kwargs["files"].outputdir / MARK_STARTED, self.name)
            logger.info(
                f"Wrote kimmdy start marker for task: {self.name} in {self.runmng.iteration}_{self.out}"
            )
        files = self.f(**self.kwargs)
        if self.out is not None:
            write_time_marker(self.kwargs["files"].outputdir / MARK_DONE, self.name)
            logger.info(
                f"Wrote kimmdy done marker for task: {self.name} in {self.runmng.iteration}_{self.out}"
            )
        logger.info(f"Finished task: {self.name}")
        if files is not None and files.logger:
            for h in files.logger.handlers:
                h.close()
        return files

    def __repr__(self) -> str:
        return str(self.name) + " args: " + str(self.kwargs)


# TODO: move this to appropriate place
def get_plumed_out(plumed: Path) -> Path:
    plumed_parsed = read_plumed(plumed)
    plumed_out = None
    for part in plumed_parsed["prints"]:
        if file := part.get("FILE"):
            if not plumed_out:
                plumed_out = file
                logger.debug(f"Found plumed print FILE {plumed_out}")
            else:
                raise RuntimeError("Multiple FILE sections found in plumed dat")
    if plumed_out is None:
        raise RuntimeError("No FILE section found in plumed dat")

    return plumed_out
