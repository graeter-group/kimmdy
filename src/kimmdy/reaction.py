from abc import ABC, abstractmethod
from dataclasses import dataclass

class Reaction(ABC)
    @ abstractmethod
    def get_rates():
        pass

    @ abstractmethod
    def get_recipe() -> ConvertionRecipe:
        pass

@dataclass
class ConvertionRecipe():
    pass


