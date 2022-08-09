from __future__ import annotations  # for 3.7 <= Python version < 3.10
from typing import TYPE_CHECKING  # fixes circular import issues for type hints

if TYPE_CHECKING:
    from kimmdy.runmanager import RunManager
    from kimmdy.config import Config
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from kimmdy.tasks import TaskFiles
import logging


class ConversionType(Enum):
    BREAK = auto()
    MOVE = auto()


@dataclass
class ConversionRecipe:
    """A ConversionReipe.
    encompasses a single transformation, e.g. moving one
    atom or braking one bond.

    Parameters
    ----------
    type : list[ConversionType.BREAK or .MOVE]
    atom_idx : list[(from, to)]
    """

    type: list[ConversionType] = field(default_factory=list)
    atom_idx: list[tuple[int, int]] = field(default_factory=list)


@dataclass
class ReactionResult:
    """A ReactionResult
    encompasses a list of transformations and their rates.

    Parameters
    ----------
    recipes : list of ConversionRecipes
    rates : list of rates
    """

    recipes: list[ConversionRecipe] = field(default_factory=list)
    rates: list[float] = field(default_factory=list)


class Reaction(ABC):
    """Reaction base class
    hast a type_scheme, which is a dict of types of possible entries in config.
    Used to read and check the input config.
    To not use this feature return empty dict.

    Example:
    ```python
    {"homolysis": {"edis": Path, "bonds": Path}}
    ```
    """

    type_scheme = dict()

    def __init__(self, name, runmng: RunManager):
        self.name = name
        self.runmng = runmng
        # sub config, settings of this specific reaction:
        self.config: Config = self.runmng.config.reactions.attr(self.name)

        logging.debug(f"Reaction {self.name} instatiated.")

    @abstractmethod
    def get_reaction_result(self, files: TaskFiles) -> ReactionResult:
        pass
