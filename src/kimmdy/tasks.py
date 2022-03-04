from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable


@dataclass
class TaskFiles:
    """Input and Output files and directories.

    Hosts the input and output files belonging to a task.
    A function or method that wants to be callable as a Task
    has to return a TaskFiles object.
    """

    input: dict[str, Path] = field(default_factory=dict)
    output: dict[str, Path] = field(default_factory=dict)
    # default outputdir is current working directory
    outputdir: Path = Path()


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


TaskMapping = dict[str, Callable[..., TaskFiles]]
