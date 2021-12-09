from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto


class ConversionType(Enum):
    BREAK = auto()
    MOVE = auto()


@dataclass
class ConversionRecipe:
    """
    atom_idx is a list with tuples indicating atoms
    being moved (from, to) a position for ConversionType.MOVE
    and and in the case of type == ConversionType.BREAK
    they indicate a bond (from, to) being broken.
    """
    type: ConversionType
    atom_type: list = field(default_factory=list)
    atom_idx: list[tuple[int, int]] = field(default_factory=list) 


@dataclass
class ReactionResult:
    recipe: ConversionRecipe
    rates: list = field(default_factory=list)


class Reaction(ABC):
    @abstractmethod
    def get_reaction_result() -> ReactionResult:
        pass

