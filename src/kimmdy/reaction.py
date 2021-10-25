from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class ConversionRecipe:
    pass


class Reaction(ABC):
    @abstractmethod
    def get_rates():
        pass

    @abstractmethod
    def get_recipe() -> ConversionRecipe:
        pass
