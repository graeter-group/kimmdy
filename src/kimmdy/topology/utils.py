from __future__ import annotations  # for 3.7 <= Python version < 3.10

import logging
from typing import TYPE_CHECKING, Callable, Optional

from gmx_top4py.topology.utils import (
    get_top_section,
    get_moleculetype_header,
    get_moleculetype_atomics,
    get_protein_section,
    get_selected_section, ######## necessary???
    set_top_section,
    attributes_to_list,
    is_not_none,
    get_by_permutations,
    match_atomic_item_to_atomic_type,
    increment_field,
    is_not_solvent_or_ion,
    get_is_selected_moleculetype_f, ####### necessary???
    get_residue_by_bonding,
)

from kimmdy.constants import ION_NAMES, REACTIVE_MOLECULEYPE, SOLVENT_NAMES
from kimmdy.topology.atomic import AtomId, MoleculeTypeHeader

if TYPE_CHECKING:
    from kimmdy.config import Config
    from kimmdy.topology.atomic import Atom, AtomicType, AtomicTypes
    from kimmdy.topology.topology import Topology

logger = logging.getLogger(__name__)


def get_reactive_section(top: dict, name: str) -> Optional[list[list]]:
    """Get content of a section in the Reactive moleculetype from a topology dict."""
    return get_top_section(top, name, moleculetype=REACTIVE_MOLECULEYPE)


def set_reactive_section(top: dict, name: str, value: list) -> Optional[list[list]]:
    """Set content of a section in the first moleculetype (protein) from a topology dict."""
    set_top_section(top, name, value, moleculetype=REACTIVE_MOLECULEYPE)

def get_is_reactive_predicate_from_config_f(cfg: Config) -> Callable[[str], bool]:
    """Returns whether a moleculetype name is configured to be recognized as reactive."""
    include = [x.lower() for x in cfg.include.split()]
    exclude = [x.lower() for x in cfg.exclude.split()]
    return get_is_reactive_predicate_f(include, exclude)


def get_is_reactive_predicate_f(
    include: list[str], exclude: list[str]
) -> Callable[[str], bool]:
    """Returns whether a moleculetype name is configured to be recognized as reactive."""
    default_excludes = [x.lower() for x in SOLVENT_NAMES + ION_NAMES]
    exclude = [s.lower() for s in exclude]
    include = [s.lower() for s in include]

    def f(name: str) -> bool:
        lower_name = name.lower()
        return lower_name not in exclude and (
            lower_name not in default_excludes or lower_name in include
        )

    return f


def get_residue_fragments(
    top: Topology,
    residue: list[Atom],
    start1: Atom,
    start2: Atom,
    iterations: int = 20,
) -> tuple[set, set]:
    """Splits a residue into fragments after a bond has been broken.

    Parameters
    ----------
    top
        Topology
    residue
        All atoms of current residue. Ok, when it contains more atoms
        as long as those are not connected to broken bond.
    start1
        First atom with broken bond
    start2
        Second atom with broken bond
    iterations
        Max number of bonds from start atoms to be included when
        building fragmets, by default 20

    Returns
    -------
        Two fragments, or one fragment and empty set in case the
        residue did not change its size.
    """
    residue_nrs = set([atom.nr for atom in residue])
    # could remove duplicate calculations
    fragments = [set([start1.nr]), set([start2.nr])]
    for fragment in fragments:
        for _ in range(iterations):
            neighbors = set()
            for nr in fragment:
                neighbors.update(set(top.atoms[nr].bound_to_nrs))
            fragment.update(neighbors.intersection(residue_nrs))
        if len(fragment) == len(residue):
            logger.debug(
                "Calculating residue fragments, but residue is whole! "
                "Could be intra-residue reaction!"
            )
            return fragment, set()
    return fragments[0], fragments[1]
