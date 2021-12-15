from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto


class ConversionType(Enum):
    BREAK = auto()
    MOVE = auto()

Topology = dict[str, list[list[str]]]

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
    type: ConversionType
    atom_idx: tuple[int, int]

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

