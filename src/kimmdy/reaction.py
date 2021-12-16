from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto


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
    type : ConversionType.BREAK or .MOVE
    atom_idx : (from, to)
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
    def get_reaction_result() -> ReactionResult:
        pass
