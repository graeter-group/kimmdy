from __future__ import annotations
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from kimmdy.runmanager import RunManager
    from kimmdy.config import Config
from abc import ABC, abstractmethod
from dataclasses import dataclass
from kimmdy.tasks import TaskFiles
import logging
from pathlib import Path
import dill
import csv

# Necessary before 3.11: https://peps.python.org/pep-0673/
from typing import TypeVar

TypeRecipe = TypeVar("TypeRecipe", bound="Recipe")


@dataclass
class RecipeStep(ABC):
    """ABC for all RecipeSteps."""

    pass


@dataclass
class Move(RecipeStep):
    """Change topology and/or coordinates to move an atom.

    Attributes
    ----------
    idx_to_move : int
    idx_to_bind : Union[int, None]
        Bonding partner to form bond with, default None.
    idx_to_break : Union[int, None]
        Bonding partner to break bond with, default None.
    new_coords : Union[list[float, float, float], None]
        Optional new xyz coordinates for atom to move to,
        default None.
    """

    idx_to_move: int
    idx_to_bind: Union[int, None] = None
    idx_to_break: Union[int, None] = None
    new_coords: Union[list[float, float, float], None] = None


@dataclass
class Break(RecipeStep):
    """Change topology to break a bond

    Attributes
    ----------
    atom_idxs : list[int, int]
        atom indices between which a bond should be removed
    """

    atom_idx_1: int
    atom_idx_2: int


@dataclass
class Bind(RecipeStep):
    """Change topology to form a bond

    Attributes
    ----------
    atom_idxs : list[int, int]
        atom indices between which a bond should be formed
    """

    atom_idx_1: int
    atom_idx_2: int


@dataclass
class Recipe:
    """A reaction path defined by one series of RecipeSteps.
    Defines everything necessart to build the
    product state from the educt state.

    Attributes
    ----------
    recipe_steps : list[RecipeStep]
        Single sequence of RecipeSteps to build product
    rates : list[float]
        Reaction rates corresponding 1:1 to times.
    times : list[float]
        List of times in ps, at which this reaction path applies.
        Must have same number of times as rates.
    avg_rates : Union[list[float], None]
        Optional, average rate corresponding to a time range, default None
    avg_timespans : Union[list[list[float,float]], None]
        Optional, per averaged rate the first and last time in ps of
        the averaged interval, default None

    """

    recipe_steps: list[RecipeStep]
    rates: list[float]
    times: list[float]
    avg_rates: Union[list[float], None] = None
    avg_timespans: Union[list[list[float, float]], None] = None

    def __post_init__(self):
        self.check_consistency()

    def calc_averages(self, window_size: int):
        """Calulate average rates over some window size

        Parameters
        ----------
        window_size : int
            Size of the window to average over,
            -1 to average over whole available range.
        """
        raise NotImplementedError

    def combine_with(self, other: TypeRecipe):
        """Combines this Recipe with another with the same RecipeSteps.

        Parameters
        ----------
        other : Recipe
        """

        if self.recipe_steps != other.recipe_steps:
            raise ValueError(
                "Error: Trying to combine reaction paths with "
                "different recipe_steps!\n"
                f"self: {self.recipe_steps}\n"
                f"other: {other.recipe_steps}"
            )

        self.check_consistency()
        other.check_consistency()

        self.rates += other.rates
        self.times += other.times

        if other.avg_rates is not None:
            if self.avg_rates is None:
                self.avg_rates = other.avg_rates
                self.avg_timespans = other.avg_timespans
            else:
                self.avg_rates += other.avg_rates
                self.avg_timespans += other.avg_timespans

    def check_consistency(self):
        """Run consistency checks for correct size of variables"""
        try:
            if len(self.rates) != len(self.times):
                raise ValueError(
                    "Times and rates are not of equal length\n"
                    f"\trates: {len(self.rates)}\n"
                    f"\times: {len(self.times)}"
                )

            if self.avg_rates is not None or self.avg_timespans is not None:
                if self.avg_rates is None or self.avg_timespans is None:
                    raise ValueError(
                        "Average times and average rates must be "
                        "of same type, but one is None\n"
                        f"\tavg_rates: {type(self.avg_rates)}\n"
                        f"\tavg_timespans: {type(self.avg_timespans)}"
                    )
                if len(self.avg_rates) != len(self.avg_timespans):
                    raise ValueError(
                        "Average times and average rates are not of equal length\n"
                        f"\tavg_rates: {len(self.avg_rates)}\n"
                        f"\tavg_timespans: {len(self.avg_timespans)}"
                    )

            double_counter = 0
            double_times = set()
            for i, time in enumerate(self.times):
                for time2 in self.times[i + 1 :]:
                    if abs(time - time2) < 1e-6:
                        double_counter += 1
                        double_times.add(time)
            if double_counter != 0:
                raise ValueError(
                    "Times are not unique! "
                    f"{double_counter} times found multiple times\n"
                    f"Times: {double_times}"
                )
        except ValueError as e:
            raise ValueError(
                f"Consistency error in Recipe {self.recipe_steps}" "" + e.args[0]
            )


@dataclass
class RecipeCollection:
    """A RecipeCollection encompasses a number of reaction paths.
    They can originate from multiple reaction plugins, but do not need to.
    """

    recipes: list[Recipe]

    def aggregate_reactions(self):
        """Combines reactions having the same sequence of RecipeSteps."""

        unique_recipes = []
        unique_recipes_idxs = []

        for i, recipe in enumerate(self.recipes):
            if recipe.recipe_steps not in unique_recipes:
                unique_recipes.append(recipe.recipe_steps)
                unique_recipes_idxs.append([i])
            else:
                for j, ur in enumerate(unique_recipes):
                    if recipe.recipe_steps == ur:
                        unique_recipes_idxs[j].append(i)

        # merge every dublicate into first reaction path
        for uri in unique_recipes_idxs:
            if len(uri) > 1:
                for uri_double in uri[1:]:
                    self.recipes[uri[0]].combine_with(self.recipes[uri_double])
        # only keep first of each reaction path
        urps = []
        for uri in unique_recipes_idxs:
            urps.append(self.recipes[uri[0]])
        self.recipes = urps

    def to_csv(self, path: Path):
        """Write a ReactionResult as defined in the reaction module to a csv file"""

        header = ["recipe_steps", "times", "rates", "avg_timespans", "avg_rates"]
        rows = []
        for i, rp in enumerate(self.recipes):
            rows.append([i] + [rp.__getattribute__(h) for h in header])
        header = ["index"] + header

        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)

    def to_dill(self, path: Path):
        with open(path, "wb") as f:
            dill.dump(self, f)

    @classmethod
    def from_dill(_, path: Path):
        with open(path, "rb") as f:
            return dill.load(f)


class ReactionPlugin(ABC):
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
    def get_recipe_collection(self, files: TaskFiles) -> RecipeCollection:
        pass
