"""Contains the Reaction Recipe, RecipeStep and RecipeCollection.
"""

from __future__ import annotations

import csv
import logging
from abc import ABC
from copy import copy
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Optional

import dill
import numpy as np

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from kimmdy.topology.topology import Topology


class RecipeStep(ABC):
    """Base class for all RecipeSteps.
    Indices can be accessed as 0-based or 1-based.
    ix: 0-based, int
    id: 1-based, str
    """


@dataclass
class Relax(RecipeStep):
    """Start a relaxation MD.

    The molecular system coordinates are far out of equilibrium after most topology changes.
    A relaxtion MD simulation using for example the slow growth method helps to reach the new equilibrium.
    """


class Place(RecipeStep):
    """Change topology and/or coordinates to place an atom.

    Either provide the index (ix_to_place) or the ID (id_to_place) of the atom to place.

    Parameters
    ----------
    new_coords :
        New xyz coordinates for atom to place to. Valid for the end point of the
        recipe timespan.

    ix_to_place : int
        Index of atom to place. 0-based.

    id_to_place : str
        Index of atom to place. 1-based
    """

    def __init__(
        self,
        new_coords: tuple[float, float, float],
        ix_to_place: Optional[int] = None,
        id_to_place: Optional[str] = None,
    ):
        """Create a new Place RecipeStep instance.

        of the first and second atom, either the id (1-based, str) or ix
        (0-based, int) can be given. If both are given, the ix is used.

        Parameters
        ----------
        new_coords :
            New xyz coordinates for atom to place to. Valid for the end point of the recipe timespan.
        ix_to_place :
            The index of the atom. zero-based, by default None
        id_to_place :
            The ID of the atom. one-based, by default None

        Raises
        ------
        ValueError
            If neither an index nor an ID is provided for any of the atoms.
        """
        self.new_coords = new_coords
        self._ix_to_place: int
        if id_to_place is not None:
            if type(id_to_place) is not str:
                raise ValueError(f"atom_id_1 is {type(id_to_place)}, should be str.")
            self._ix_to_place = int(id_to_place) - 1
        if ix_to_place is not None:
            if type(ix_to_place) is not int:
                raise ValueError(f"atom_ix_1 is {type(ix_to_place)}, should be int.")
            self._ix_to_place = ix_to_place
        if id_to_place is not None and ix_to_place is not None:
            logger.warning(f"Both atom_ix_1 and atom_id_1 are given, using atom_ix_1.")
        if self._ix_to_place is None:
            raise ValueError("id_ or ix_ to_place must be provided!")

    @property
    def ix_to_place(self) -> Optional[int]:
        return self._ix_to_place

    @ix_to_place.setter
    def ix_to_place(self, value: Optional[int]):
        if isinstance(value, property):
            return
        assert isinstance(value, int), f"ix_to_place is {type(value)}, should be int."
        self._ix_to_place = value

    @property
    def id_to_place(self) -> Optional[str]:
        return str(self._ix_to_place + 1)

    @id_to_place.setter
    def id_to_place(self, value: Optional[str]):
        if isinstance(value, property):
            return
        assert isinstance(value, str), f"id_to_place is {type(value)}, should be str."
        self._ix_to_place = int(value) - 1

    def __eq__(self, other):
        """Two Placements are equal if their atom indices are equal and they have the same new coordinates."""

        if not isinstance(other, Place) or type(self) != type(other):
            logger.warning(
                f"Comparing RecipeSteps with different types: {type(self)} and {type(other)}. Returning False."
            )
            return False
        return (
            self._ix_to_place == other._ix_to_place
            and self.new_coords == other.new_coords
        )

    def __hash__(self):
        return hash((self._ix_to_place, self.new_coords, type(self).__name__))

    def __repr__(self):
        return f"{type(self).__name__}(ix_to_place={self._ix_to_place}, new_coords={self.new_coords})"

    def __str__(self):
        return f"{type(self).__name__}({self._ix_to_place}, {self.new_coords})"


class BondOperation(RecipeStep):
    """Handle a bond operation on the recipe step.

    This class takes in either zero-based indices or one-base IDs for two atoms.

    Parameters
    ----------
    atom_ix_1 : int, optional
        The index of the first atom. zero-based, by default None
    atom_ix_2 : int, optional
        The index of the second atom. zero-based, by default None
    atom_id_1 : str, optional
        The ID of the first atom. one-based, by default None
    atom_id_2 : str, optional
        The ID of the second atom. one-based, by default None

    Raises
    ------
    ValueError
        If neither an index nor an ID is provided for any of the atoms.

    Notes
    -----
    Internally, this class stores the atom indices and converts IDs to indices as needed.

    """

    def __init__(
        self,
        atom_ix_1: Optional[int] = None,
        atom_ix_2: Optional[int] = None,
        atom_id_1: Optional[str] = None,
        atom_id_2: Optional[str] = None,
    ):
        """Create a new BondOperation instance.

        of the first and second atom, either the id (1-based, str) or ix
        (0-based, int) can be given. If both are given, the ix is used.

        Parameters
        ----------
        atom_ix_1 : int, optional
            The index of the first atom. zero-based, by default None
        atom_ix_2 : int, optional
            The index of the second atom. zero-based, by default None
        atom_id_1 : str, optional
            The ID of the first atom. one-based, by default None
        atom_id_2 : str, optional
            The ID of the second atom. one-based, by default None

        Raises
        ------
        ValueError
            If neither an index nor an ID is provided for any of the atoms.
        """
        self._atom_ix_1: int
        self._atom_ix_2: int
        if atom_id_1 is not None:
            if type(atom_id_1) is not str:
                raise ValueError(f"atom_id_1 is {type(atom_id_1)}, should be str.")
            self._atom_ix_1 = int(atom_id_1) - 1
        if atom_id_2 is not None:
            if type(atom_id_2) is not str:
                raise ValueError(f"atom_id_2 is {type(atom_id_2)}, should be str.")
            self._atom_ix_2 = int(atom_id_2) - 1
        if atom_ix_1 is not None:
            if type(atom_ix_1) is not int:
                raise ValueError(f"atom_ix_1 is {type(atom_ix_1)}, should be int.")
            self._atom_ix_1 = atom_ix_1
        if atom_ix_2 is not None:
            if type(atom_ix_2) is not int:
                raise ValueError(f"atom_ix_2 is {type(atom_ix_2)}, should be int.")
            self._atom_ix_2 = atom_ix_2
        if atom_ix_1 is not None and atom_id_1 is not None:
            logger.warning(f"Both atom_ix_1 and atom_id_1 are given, using atom_ix_1.")

    @property
    def atom_id_1(self) -> str:
        return str(self._atom_ix_1 + 1)

    @atom_id_1.setter
    def atom_id_1(self, value: str):
        if isinstance(value, property):
            return
        assert isinstance(value, str), f"atom_id_1 is {type(value)}, should be str."
        self._atom_ix_1 = int(value) - 1

    @property
    def atom_ix_1(self) -> int:
        return self._atom_ix_1

    @atom_ix_1.setter
    def atom_ix_1(self, value: int):
        if isinstance(value, property):
            return
        assert isinstance(value, int), f"atom_ix_1 is {type(value)}, should be int."
        self._atom_ix_1 = value

    @property
    def atom_id_2(self) -> str:
        return str(self._atom_ix_2 + 1)

    @atom_id_2.setter
    def atom_id_2(self, value: str):
        if isinstance(value, property):
            return
        assert isinstance(value, str), f"atom_id_2 is {type(value)}, should be str."
        self._atom_ix_2 = int(value) - 1

    @property
    def atom_ix_2(self) -> int:
        return self._atom_ix_2

    @atom_ix_2.setter
    def atom_ix_2(self, value: int):
        if isinstance(value, property):
            return
        assert isinstance(value, int), f"atom_ix_2 is {type(value)}, should be int."
        self._atom_ix_2 = value

    def __eq__(self, other):
        """Two BondOperations are equal if their atom indices are equal and they are of the same type."""

        if not isinstance(other, BondOperation) or type(self) != type(other):
            logger.warning(
                f"Comparing RecipeSteps with different types: {type(self)} and {type(other)}. Returning False."
            )
            return False
        return (
            self._atom_ix_1 == other._atom_ix_1 and self._atom_ix_2 == other._atom_ix_2
        )

    def __hash__(self):
        return hash((self._atom_ix_1, self._atom_ix_2, type(self).__name__))

    def __repr__(self):
        return f"{type(self).__name__}(atom_ix_1={self._atom_ix_1}, atom_ix_2={self._atom_ix_2})"

    def __str__(self):
        return f"{type(self).__name__}({self._atom_ix_1}, {self._atom_ix_2})"


class Break(BondOperation):
    """Change topology to break a bond

    Parameters
    ----------
    atom_ix_1 : int, optional
        The index of the first atom. zero-based, by default None
    atom_ix_2 : int, optional
        The index of the second atom. zero-based, by default None
    atom_id_1 : str, optional
        The ID of the first atom. one-based, by default None
    atom_id_2 : str, optional
        The ID of the second atom. one-based, by default None
    """


class Bind(BondOperation):
    """Change topology to form a bond

    Parameters
    ----------
    atom_ix_1 : int, optional
        The index of the first atom. zero-based, by default None
    atom_ix_2 : int, optional
        The index of the second atom. zero-based, by default None
    atom_id_1 : str, optional
        The ID of the first atom. one-based, by default None
    atom_id_2 : str, optional
        The ID of the second atom. one-based, by default None
    """


@dataclass
class CustomTopMod(RecipeStep):
    """A custom recipe step that can be used to define a custom topology modification.

    Parameters
    ----------
    f : Callable[[Topology], Topology]
        A function that takes a Topology object and modifies it in place.
    """

    f: Callable[[Topology], None]


@dataclass
class Recipe:
    """A reaction path defined by one series of RecipeSteps.
    Defines everything necessart to build the
    product state from the educt state.

    Parameters
    ----------
    recipe_steps : list[RecipeStep]
        Single sequence of RecipeSteps to build product
    rates : list[float]
        Reaction rates corresponding 1:1 to timespans.
    timespans : list[list[float, float]]
        List of half-open timespans (t1, t2] in ps, at which this rate is valid.
        Recipe steps which change the coordinates only need to be applicable
        at the first time in the interval.
        Must have same number of timespans as rates.
        t1 can equal t2 for the last frame.

    """

    recipe_steps: list[RecipeStep]
    rates: list[float]
    timespans: list[tuple[float, float]]

    def __post_init__(self):
        self.check_consistency()

    def copy(self):
        return Recipe(
            list(self.recipe_steps),
            list(self.rates),
            list(self.timespans),
        )

    def combine_with(self, other: Recipe):
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
        self.timespans += other.timespans

    def check_consistency(self):
        """Run consistency checks for correct size of variables"""
        try:
            if len(self.rates) != len(self.timespans):
                raise ValueError(
                    "Timespans and rates are not of equal length\n"
                    f"\trates: {len(self.rates)}\n"
                    f"\ttimespans: {len(self.timespans)}"
                )
            if not isinstance(self.recipe_steps, list):
                raise ValueError(
                    f"Recipe Steps must be a list, not {type(self.recipe_steps)}"
                )
            if not isinstance(self.timespans, list):
                raise ValueError(
                    f"timespans must be a list, not {type(self.timespans)}"
                )
            if not isinstance(self.rates, list):
                raise ValueError(f"rates must be a list, not {type(self.rates)}")

        except ValueError as e:
            raise ValueError(
                f"Consistency error in Recipe {self.recipe_steps}\n" + e.args[0]
            )

    def get_recipe_name(self):
        name = ""
        for rs in self.recipe_steps:
            name += " "
            if isinstance(rs, Place):
                if (ix := getattr(rs, "ix_to_place", None)) is not None:
                    name += str(ix)
                    name += "\u27A1"  # ➡
            elif isinstance(rs, Bind):
                if (ix := getattr(rs, "atom_ix_1", None)) is not None:
                    name += str(ix)
                    name += "\u27A1"  # ➡
                if (ix := getattr(rs, "atom_ix_2", None)) is not None:
                    name += str(ix)

            elif isinstance(rs, Break):
                if (ix := getattr(rs, "atom_ix_1", None)) is not None:
                    name += str(ix)
                    name += "\u26A1"  # ➡
                if (ix := getattr(rs, "atom_ix_2", None)) is not None:
                    name += str(ix)

            elif isinstance(rs, Relax):
                pass

            elif isinstance(rs, CustomTopMod):
                name += "🪛"

            else:
                logger.warning(f"get_recipe_name got unknown step type: {type(rs)}")
                name += "?".join(list(map(str, rs.__dict__.values())))
        return name


@dataclass
class RecipeCollection:
    """A RecipeCollection encompasses a number of reaction paths.
    They can originate from multiple reaction plugins, but do not need to.

    Returns
    -------
    unique_recipes_ixs
        List of lists binning the old recipe indices and maps to the new
        recipes list.
    """

    recipes: list[Recipe]

    def aggregate_reactions(self):
        """Combines reactions having the same sequence of RecipeSteps."""

        logger.debug(f"Aggregating {len(self.recipes)} recipes")
        unique_recipes = []
        unique_recipes_ixs = []

        for i, recipe in enumerate(self.recipes):
            if recipe.recipe_steps not in unique_recipes:
                unique_recipes.append(recipe.recipe_steps)
                unique_recipes_ixs.append([i])
            else:
                for j, ur in enumerate(unique_recipes):
                    if recipe.recipe_steps == ur:
                        unique_recipes_ixs[j].append(i)

        # merge every dublicate into first reaction path
        for uri in unique_recipes_ixs:
            # copy to not change outside references to this recip
            self.recipes[uri[0]] = self.recipes[uri[0]].copy()
            if len(uri) > 1:
                for uri_double in uri[1:]:
                    self.recipes[uri[0]].combine_with(self.recipes[uri_double])
        # only keep first of each reaction path
        urps = []
        for uri in unique_recipes_ixs:
            urps.append(self.recipes[uri[0]])
        self.recipes = urps
        return unique_recipes_ixs

    def to_csv(self, path: Path, picked_recipe=None):
        """Write a ReactionResult as defined in the reaction module to a csv file"""

        header = ["recipe_steps", "timespans", "rates"]
        rows = []
        for i, rp in enumerate(self.recipes):
            picked = rp == picked_recipe
            rows.append([i] + [picked] + [rp.__getattribute__(h) for h in header])
        header = ["index", "picked"] + header

        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)

    @classmethod
    def from_csv(cls, path: Path):
        """Create a RecipeCollection object from a CSV file
        Returns the recipe collection and a single recipe that was picked, otherwise None
        """
        recipes = []

        with open(path, "r") as f:
            reader = csv.reader(f)
            _ = next(reader)  # Skip the header row
            picked_rp = None
            for row in reader:
                _, picked_s, recipe_steps_s, timespans_s, rates_s = row

                picked = eval(picked_s)
                recipe_steps = eval(recipe_steps_s)
                timespans = eval(timespans_s)
                rates = eval(rates_s)

                recipe = Recipe(
                    recipe_steps=recipe_steps, timespans=timespans, rates=rates
                )
                recipes.append(recipe)
                if picked:
                    picked_rp = recipe

        return cls(recipes=recipes), picked_rp

    def to_dill(self, path: Path):
        with open(path, "wb") as f:
            dill.dump(self, f)

    @classmethod
    def from_dill(cls, path: Path):
        with open(path, "rb") as f:
            return dill.load(f)

    def calc_cumprob(self):
        """Calculate cumulative probability of all contained recipe steps.
        Sums up to 1 over all recipes. Assumes constant rate for given timespan
        and rate zero otherwise.
        """

        self.aggregate_reactions()

        cumprob = []

        for re in self.recipes:
            integral = 0
            for t, r in zip(re.timespans, re.rates):
                dt = t[1] - t[0]
                integral += r * dt
            cumprob.append(integral)

        return np.array(cumprob) / sum(cumprob)

    def calc_ratesum(self) -> tuple[list, list, list]:
        """Calculate the sum of rates over all timesteps

        Returns
        -------
        boarders
            flat list containing times of rate changes marking the
            boarders of the windows
        rate_windows
            flat list containing all rates in between the boarders.
            Each window is orderd as in recipe_windows
        recipe_windows
            flat list containing the recipes of the corresponding window.
            Each window is orderd as in rate_windows
        """

        self.aggregate_reactions()

        # non-overlapping window boaders
        boarders = set()

        for re in self.recipes:
            for ts in re.timespans:
                for t in ts:
                    boarders.add(t)

        boarders = sorted(boarders)
        rate_windows = [[] for _ in range(len(boarders) - 1)]
        recipe_windows = [[] for _ in range(len(boarders) - 1)]

        for re in self.recipes:
            for r, ts in zip(re.rates, re.timespans):
                left_idx = boarders.index(ts[0])
                right_idx = boarders.index(ts[1])
                # timespan <- boarder
                # the selected recipe must tell where it should be applied
                re_copy = copy(re)
                re_copy.timespans = [ts]
                re_copy.rates = [r]
                [l.append(r) for l in rate_windows[left_idx:right_idx]]
                [l.append(re_copy) for l in recipe_windows[left_idx:right_idx]]

        return boarders, rate_windows, recipe_windows

    def plot(
        self,
        outfile,
        highlight_r: Optional[Recipe] = None,
        highlight_t: Optional[float] = None,
    ):
        """Plot reaction rates over time

        Parameters
        ----------
        outfile : Path
            Where to save the plot, must have compatible suffix.
        highlight_r : Recipe, optional
            Recipe to highlight, by default None
        highlight_t : float, optional
            Time at which the reactions starts
        """

        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sns

        self.aggregate_reactions()
        cumprob = self.calc_cumprob()
        recipes = np.array(self.recipes)
        recipe_steps = np.empty(len(self.recipes), dtype=object)
        for i, r in enumerate(self.recipes):
            recipe_steps[i] = r.recipe_steps

        idxs = list(np.argsort(cumprob)[-8:])
        ref = np.empty((1,), dtype=object)
        if highlight_r is not None:
            ref[0] = highlight_r.recipe_steps
            i_to_highlight = np.nonzero(recipe_steps == ref)[0]
            idxs = list(set(np.concatenate([idxs, i_to_highlight])))

        cmap = sns.color_palette("husl", len(idxs))
        name_to_args = {}

        for r_i, re in enumerate(recipes):
            name = re.get_recipe_name()
            kwargs = {}
            kwargs["linestyle"] = "-"
            kwargs["linewidth"] = 1.0

            if name not in name_to_args:
                kwargs["color"] = (0.7, 0.7, 0.7)
                kwargs["label"] = "REMOVE"
                if r_i in idxs:
                    kwargs["color"] = cmap[idxs.index(r_i)]
                    kwargs["label"] = name
                name_to_args[name] = kwargs
            else:
                if r_i in idxs:
                    if name_to_args[name]["label"] == "REMOVE":
                        kwargs["color"] = cmap[idxs.index(r_i)]
                        kwargs["label"] = name
                        name_to_args[name] = kwargs

        plt.figure()
        if highlight_t is not None:
            plt.axvline(highlight_t, color="red")
        for r_i, re in enumerate(recipes):
            name = re.get_recipe_name()
            if highlight_r is not None and name == highlight_r.get_recipe_name():
                name_to_args[name]["linestyle"] = "-."
                name_to_args[name]["linewidth"] = 2.2

            for dt, r in zip(re.timespans, re.rates):
                marker = ""
                if dt[0] == dt[1]:
                    marker = "."
                plt.plot(
                    np.array(dt),
                    (r, r),
                    marker=marker,
                    **name_to_args[name],
                )
        plt.xlabel("time [ps]")
        plt.ylabel("reaction rate")
        plt.yscale("log")
        # removing duplicates in legend
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        by_label.pop("REMOVE", None)
        plt.legend(by_label.values(), by_label.keys())

        plt.savefig(outfile)
