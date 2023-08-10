from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, TYPE_CHECKING, Union
import logging


class AutoFillDict(dict):
    def __init__(self, get_missing: Callable):
        self.get_missing = get_missing

    def __missing__(self, key):
        self[key] = self.get_missing(key)
        return self[key]


@dataclass
class TaskFiles:
    """Input and Output files and directories.

    Hosts the input and output files belonging to a task.
    A function or method that wants to be callable as a Task
    has to return a TaskFiles object.
    The input defaultdict is populated on the fly using
    get_latest of the runmanager to find newest files.
    Files which can not be found by get_latest must be added manually.

    Examples
    --------
    >>> class run():
    >>>     def get_latest(self, s):
    >>>         return f"latest {s}"
    >>> runmng = run()
    >>> files = TaskFiles(runmng)
    >>> files.input
    >>> files.input["top"]
    {'top': 'latest top'}
    """

    get_latest: Callable
    input: dict[str, Path] = field(default_factory=dict)
    output: dict[str, Path] = field(default_factory=dict)
    outputdir: Path = Path()

    def __post_init__(self):
        self.input = AutoFillDict(self.get_latest)


def create_task_directory(runmng, postfix: str) -> TaskFiles:
    """Creates TaskFiles object, output directory and symlinks ff."""
    files = TaskFiles(runmng.get_latest)
    runmng.iteration += 1
    files.outputdir = runmng.config.out / f"{runmng.iteration}_{postfix}"
    logging.debug(f"Creating Output directory: {files.outputdir}")
    files.outputdir.mkdir(exist_ok=runmng.from_checkpoint)
    if not (files.outputdir / runmng.config.ff.name).exists():
        (files.outputdir / runmng.config.ff.name).symlink_to(runmng.config.ff)
    return files


class Task:
    """A task to be performed as as a step in the RunManager.

    A task consists of a function and its keyword arguments.
    Calling a taks calls the stored function.
    The function must return a TaskFiles object.

    Parameters:
    -----------
    runmng : kimmdy.runmanager.Runmanager
        Runmanager instance
    f : Callable
        Will be called when the task is called
    kwargs : dict
        kwargs will be passed to f
    out : str, optional
        If not None, an output dir will be created with this name
    """

    def __init__(self, runmng, f: Callable[..., TaskFiles], kwargs=None, out=None):
        self.runmng = runmng
        self.f = f
        if kwargs is None:
            kwargs = {}
        self.kwargs = kwargs
        self.name = self.f.__name__
        self.out = out

        logging.info(
            f"Init task {self.name}\n\tkwargs: {self.kwargs}\n\tOut: {self.out}"
        )

    def __call__(self) -> TaskFiles:
        if self.out is not None:
            self.kwargs.update({"files": create_task_directory(self.runmng, self.out)})

        logging.debug(f"Calling task {self.name} with kwargs: {self.kwargs}")
        return self.f(**self.kwargs)

    def __repr__(self) -> str:
        return str(self.name) + " args: " + str(self.kwargs)


TaskMapping = dict[
    str,
    Union[list[Callable[..., Optional[TaskFiles]]], Callable[..., Optional[TaskFiles]]],
]
