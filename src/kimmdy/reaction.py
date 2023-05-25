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

TypeReactionPath = TypeVar("TypeReactionPath", bound="ReactionPath")


@dataclass
class Conversion(ABC):
    """ABC for all conversions.
    """    
    pass


@dataclass
class Move(Conversion):
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
class Break(Conversion):
    """Change topology to break a bond

    Attributes
    ----------
    atom_idxs : list[int, int]
        atom indices between which a bond should be removed
    """

    atom_idxs: list[int, int]


@dataclass
class Bind(Conversion):
    """Change topology to form a bond

    Attributes
    ----------
    atom_idxs : list[int, int]
        atom indices between which a bond should be formed
    """

    atom_idxs: list[int, int]


@dataclass
class ReactionPath:
    """A reaction path defined by one series of conversions.
    Defines everything necessart to build the
    product state from the educt state.

    Attributes
    ----------
    conversions : list[Conversion]
        Single sequence of conversions to build product
    rates : list[float]
        Reaction rates corresponding 1:1 to frames.
    frames : list[int]
        List of frame indices, in which this reaction path applies.
        Must have same number of frames as rates.
    avg_rates : Union[list[float], None]
        Optional, average rate corresponding to a frame range, default None
    avg_frames : Union[list[list[int,int]], None]
        Optional, per averaged rate the first and last frame index of
        the averaged interval, default None

    """

    conversions: list[Conversion]
    rates: list[float]
    frames: list[int]
    avg_rates: Union[list[float], None] = None
    avg_frames: Union[list[list[int, int]], None] = None

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

    def combine_with(self, other: TypeReactionPath):
        """Combines this ReactionPath with another with the same conversions.

        Parameters
        ----------
        other : ReactionPath
        """

        if self.conversions != other.conversions:
            raise ValueError(
                "Error: Trying to combine reaction paths with "
                "different conversions!\n"
                f"self: {self.conversions}\n"
                f"other: {other.conversions}"
            )

        self.check_consistency()
        other.check_consistency()

        self.rates += other.rates
        self.frames += other.frames

        if other.avg_rates is not None:
            if self.avg_rates is None:
                self.avg_rates = other.avg_rates
                self.avg_frames = other.avg_frames
            else:
                self.avg_rates += other.avg_rates
                self.avg_frames += other.avg_frames

    def check_consistency(self):
        """Run consistency checks for correct size of variables"""
        try:
            if len(self.rates) != len(self.frames):
                raise ValueError(
                    "Frames and rates are not of equal length\n"
                    f"\trates: {len(self.rates)}\n"
                    f"\tframes: {len(self.frames)}"
                )

            if self.avg_rates is not None or self.avg_frames is not None:
                if self.avg_rates is None or self.avg_frames is None:
                    raise ValueError(
                        "Average frames and average rates must be "
                        "of same type, but one is None\n"
                        f"\tavg_rates: {type(self.avg_rates)}\n"
                        f"\tavg_frames: {type(self.avg_frames)}"
                    )
                if len(self.avg_rates) != len(self.avg_frames):
                    raise ValueError(
                        "Average frames and average rates are not of equal length\n"
                        f"\tavg_rates: {len(self.avg_rates)}\n"
                        f"\tavg_frames: {len(self.avg_frames)}"
                    )

            double_counter = 0
            double_frames = set()
            for i, frame in enumerate(self.frames):
                for frame2 in self.frames[i + 1 :]:
                    if frame == frame2:
                        double_counter += 1
                        double_frames.add(frame)
            if double_counter != 0:
                raise ValueError(
                    "Frames are not unique! "
                    f"{double_counter} frames found multiple times\n"
                    f"Frames: {double_frames}"
                )
        except ValueError as e:
            raise ValueError(
                f"Consistency error in ReactionPath {self.conversions}" "" + e.args[0]
            )


@dataclass
class ReactionResults:
    """A ReactionResults encompasses a number of reaction paths.
    They can originate from multiple reaction plugins, but do not need to.
    """

    reaction_paths: list[ReactionPath]

    def aggregate_reactions(self):
        """Combines reactions having the same sequence of conversions."""

        unique_conversions = []
        unique_conversions_idxs = []

        for i, rp in enumerate(self.reaction_paths):
            if rp.conversions not in unique_conversions:
                unique_conversions.append(rp.conversions)
                unique_conversions_idxs.append([i])
            else:
                for j, uc in enumerate(unique_conversions):
                    if rp.conversions == uc:
                        unique_conversions_idxs[j].append(i)

        # merge every dublicate into first reaction path
        for uci in unique_conversions_idxs:
            if len(uci) > 1:
                for uci_double in uci[1:]:
                    self.reaction_paths[uci[0]].combine_with(
                        self.reaction_paths[uci_double]
                    )
        # only keep first of each reaction path
        urps = []
        for uci in unique_conversions_idxs:
            urps.append(self.reaction_paths[uci[0]])
        self.reaction_paths = urps

    def to_csv(self, path: Path):
        """Write a ReactionResult as defined in the reaction module to a csv file"""

        header = ["conversions", "frames", "rates", "avg_frames", "avg_rates"]
        rows = []
        for i, rp in enumerate(self.reaction_paths):
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
    def get_reaction_results(self, files: TaskFiles) -> ReactionResults:
        pass
