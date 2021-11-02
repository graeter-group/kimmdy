from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto

class ConversionType(Enum):
    BREAK = auto()
    MOVE = auto()


@dataclass
class ConversionRecipe:
    type: ConversionType
    atom_idx: list[tuple[int, int]]
  
class Reaction(ABC):
    @abstractmethod
    def get_rates():
        pass

    @abstractmethod
    def get_recipe() -> ConversionRecipe:
        pass
