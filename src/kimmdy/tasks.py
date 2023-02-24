from dataclasses import dataclass, field, InitVar
from pathlib import Path
from typing import Callable, Optional


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
    >>> files.input["tpr"]
    {'top': 'latest top'}

    """

    runmng: InitVar
    input: dict[str, Path] = field(default_factory=dict)
    output: dict[str, Path] = field(default_factory=dict)
    # default outputdir is current working directory
    outputdir: Path = Path()

    def __post_init__(self, runmng):
        self.input = AutoFillDict(runmng.get_latest)


class Task:
    """A task to be performed as as a step in the RunManager.

    A task consists of a function and its keyword arguments and is
    itself callable. The function must return a TaskFiles object.
    """

    def __init__(self, f: Callable[..., TaskFiles], kwargs={}):
        self.f = f
        self.kwargs = kwargs
        self.name = self.f.__name__

    def __call__(self) -> TaskFiles:
        return self.f(**self.kwargs)

    def __repr__(self) -> str:
        return str(self.f) + " args: " + str(self.kwargs)


TaskMapping = dict[str, list[Callable[..., Optional[TaskFiles]]]]
