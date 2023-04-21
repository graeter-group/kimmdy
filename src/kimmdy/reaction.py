from __future__ import annotations
from collections.abc import Callable  # for 3.7 <= Python version < 3.10
from typing import TYPE_CHECKING

from kimmdy.topology.topology import (
    Topology,
)  # fixes circular import issues for type hints

if TYPE_CHECKING:
    from kimmdy.runmanager import RunManager
    from kimmdy.config import Config
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto
from kimmdy.tasks import TaskFiles
import logging


class ConversionType(Enum):
    BREAK = auto()
    BIND = auto()
    MOVE = auto()


@dataclass
class Conversion:
    """A Conversion.
    encompasses a single transformation, breaking of adding one bond.

    Parameters
    ----------
    type : ConversionType.BREAK or .BIND
    atom_idx : tuple(from, to)
    """

    type: ConversionType
    atom_idx: tuple[str, str]


ConversionRecipe = list[Conversion]
"""A ConversionRecipe.

is a list of Conversions to encompass one reaction outcome.
In the case of breaking a bond it is simply a list of length one,
but for e.g. moving an atom from one binding partner to another
it is a list with one BREAK and one BIND operation.
"""


@dataclass
class ReactionOutcome:
    """A ReactionOutcome
    encompasses a ConversionRecipe and its rate.
    """

    recipe: ConversionRecipe
    rate: float


ReactionResult = list[ReactionOutcome]
"""A ReactionResult
encompasses a list of ReactionOutcomes.
Each outcome has a ConversionRecipe for changing the topology and a rate
at which it is predicted to happen.
"""


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
        # TODO: does it even need it's own config explicitlye if if can always access the full config anyways?
        self.config: Config = self.runmng.config.reactions.attr(self.name)

        logging.debug(f"Reaction {self.name} instatiated.")

    @abstractmethod
    def get_reaction_result(self, files: TaskFiles) -> ReactionResult:
        pass
