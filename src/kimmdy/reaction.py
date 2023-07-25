"""
ReactionPlugin protocoll and reaction recipes.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from kimmdy.runmanager import RunManager
    from kimmdy.config import Config
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, InitVar
from kimmdy.tasks import TaskFiles
import logging
from pathlib import Path
import dill
import csv


class RecipeStep(ABC):
    """Base class for all RecipeSteps.
    Indices can be accessed as 0-based or 1-based.
    ix: 0-based, int
    id: 1-based, str
    """


@dataclass
class Move(RecipeStep):
    """Change topology and/or coordinates to move an atom.

    Parameters
    ----------
    ix_to_move : int
        Index of atom to move. 0-based.
    ix_to_bind : int
        Bonding partner to form bond with. 0-based.
    ix_to_break : int
        Bonding partner to break bond with, default None. 0-based.
    new_coords :
        Optional new xyz coordinates for atom to move to, and the associated
        time in ps default None.
    id_to_move : str
        Index of atom to move. 1-based
    id_to_bind : str
        Bonding partner to form bond with. 1-based
    id_to_break : str
        Bonding partner to break bond with, default None. 1-based
    """

    ix_to_move: Optional[int] = None
    ix_to_bind: Optional[int] = None
    ix_to_break: Optional[int] = None
    new_coords: Optional[tuple[tuple[float, float, float], float]] = None
    id_to_move: Optional[str] = None
    id_to_bind: Optional[str] = None
    id_to_break: Optional[str] = None

    _ix_to_move: Optional[int] = field(init=False, repr=False, default=None)
    _ix_to_bind: Optional[int] = field(init=False, repr=False, default=None)
    _ix_to_break: Optional[int] = field(init=False, repr=False, default=None)

    # During init without given parameters, setters recive **not** the default
    # value, but a property instance -> Error in type conversion

    @property
    def ix_to_move(self) -> Optional[int]:
        return self._ix_to_move

    @ix_to_move.setter
    def ix_to_move(self, value: Optional[int]):
        if isinstance(value, property):
            return
        self._ix_to_move = value

    @property
    def id_to_move(self) -> Optional[str]:
        if self._ix_to_move is None:
            return None
        return str(self._ix_to_move + 1)

    @id_to_move.setter
    def id_to_move(self, value: Optional[str]):
        if isinstance(value, property):
            return
        self._ix_to_move = int(value) - 1

    @property
    def ix_to_bind(self) -> Optional[int]:
        return self._ix_to_bind

    @ix_to_bind.setter
    def ix_to_bind(self, value: Optional[int]):
        if isinstance(value, property):
            return
        self._ix_to_bind = value

    @property
    def id_to_bind(self) -> Optional[str]:
        if self._ix_to_bind is None:
            return None
        return str(self._ix_to_bind + 1)

    @id_to_bind.setter
    def id_to_bind(self, value: Optional[str]):
        if isinstance(value, property):
            return
        self._ix_to_bind = int(value) - 1

    @property
    def ix_to_break(self) -> Optional[int]:
        return self._ix_to_break

    @ix_to_break.setter
    def ix_to_break(self, value: Optional[int]):
        if isinstance(value, property):
            return
        self._ix_to_break = value

    @property
    def id_to_break(self) -> Optional[str]:
        if self._ix_to_break is None:
            return None
        return str(self._ix_to_break + 1)

    @id_to_break.setter
    def id_to_break(self, value: Optional[str]):
        if isinstance(value, property):
            return
        self._ix_to_break = int(value) - 1


@dataclass
class SingleOperation(RecipeStep):
    """Handle a single operation on the recipe step.

    This class takes in either zero-based indices or one-base IDs for two atoms

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

    atom_ix_1: Optional[int] = None
    atom_ix_2: Optional[int] = None
    atom_id_1: Optional[str] = None
    atom_id_2: Optional[str] = None

    _atom_ix_1: int = field(init=False, repr=False, default=None)
    _atom_ix_2: int = field(init=False, repr=False, default=None)

    def __post_init__(self):
        e = ""
        if not (isinstance(self.atom_ix_1, int) or isinstance(self.atom_id_1, str)):
            e += "Exactly on of atom_ix_1 and atom_id_1 must be given!\n"
        if not (isinstance(self.atom_ix_2, int) or isinstance(self.atom_id_2, str)):
            e += "Exactly on of atom_ix_2 and atom_id_2 must be given!\n"
        if len(e) > 0:
            raise ValueError(e)

    @property
    def atom_id_1(self) -> str:
        if self._atom_ix_1 is None:
            return None
        return str(self._atom_ix_1 + 1)

    @atom_id_1.setter
    def atom_id_1(self, value: str):
        if isinstance(value, property):
            return
        self._atom_ix_1 = int(value) - 1

    @property
    def atom_ix_1(self) -> int:
        return self._atom_ix_1

    @atom_ix_1.setter
    def atom_ix_1(self, value: int):
        if isinstance(value, property):
            return
        self._atom_ix_1 = value

    @property
    def atom_id_2(self) -> str:
        if self._atom_ix_2 is None:
            return None
        return str(self._atom_ix_2 + 1)

    @atom_id_2.setter
    def atom_id_2(self, value: str):
        if isinstance(value, property):
            return
        self._atom_ix_2 = int(value) - 1

    @property
    def atom_ix_2(self) -> int:
        return self._atom_ix_2

    @atom_ix_2.setter
    def atom_ix_2(self, value: int):
        if isinstance(value, property):
            return
        self._atom_ix_2 = value


class Break(SingleOperation):
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


class Bind(SingleOperation):
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
        List of half-open timespans (t1, t2] in ps, at which this reaction
        path applies. Must have same number of timespans as rates.
        t1 can equal t2 for the first frame.

    """

    recipe_steps: list[RecipeStep]
    rates: list[float]
    timespans: list[tuple[float, float]]

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
        raise NotImplementedError("calc_averages not implemented yet")

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
                    f"\timespans: {len(self.timespans)}"
                )
            if len(self.rates) < 1:
                raise ValueError("Recipe empty! Use empty RecipeCollection instead!")

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
            if len(uri) > 1:
                for uri_double in uri[1:]:
                    self.recipes[uri[0]].combine_with(self.recipes[uri_double])
        # only keep first of each reaction path
        urps = []
        for uri in unique_recipes_ixs:
            urps.append(self.recipes[uri[0]])
        self.recipes = urps

    def to_csv(self, path: Path):
        """Write a ReactionResult as defined in the reaction module to a csv file"""

        header = ["recipe_steps", "timespans", "rates"]
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
    def from_dill(cls, path: Path):
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

    Parameters
    ----------
    name :
        Name of the reaction
    runmng :
        RunManager instance
    type_scheme : dict
        dict of types of possible entries in config
    """

    type_scheme: dict = dict()

    def __init__(self, name: str, runmng: RunManager):
        self.name = name
        self.runmng = runmng
        # sub config, settings of this specific reaction:
        self.config: Config = self.runmng.config.reactions.__getattribute__(self.name)

        logging.debug(f"Reaction {self.name} instatiated.")

    @abstractmethod
    def get_recipe_collection(self, files: TaskFiles) -> RecipeCollection:
        """Get a RecipeCollection as a result of the reaction.

        This is run as a [Task](`kimmdy.task.Task`) in the RunManager.
        How the RecipeCollection is built is up to the reaction.
        It has access to the current state of the system via the
        runmanager `self.runmng` and the files.

        Parameters
        ----------
        files :
            TaskFiles instance
        """
        pass
