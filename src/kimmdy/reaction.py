from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from kimmdy.tasks import TaskFiles

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
    @abstractmethod
    def get_reaction_result(self, files: TaskFiles) -> ReactionResult:
        pass

    @property
    def type_scheme(self) -> dict:
        """Dict of types of possible entries in config.
        Used to read and check the input config.
        To not use this feature return empty dict.

        Example:
        ```python
        {"homolysis": {"edis": Path, "bonds": Path}}
        ```
        """
        return dict()


