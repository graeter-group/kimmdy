import logging
from kimmdy.reaction import ConversionRecipe, ConversionType
from kimmdy.parsing import read_plumed, write_plumed, read_topol, write_topol, Topology
from pathlib import Path


def modify_top(recipe: ConversionRecipe, oldtop: Path, newtop: Path):
    logging.info(f"Reading: {oldtop} and writing modified topology to {newtop}.")
    topology = read_topol(oldtop)

    for type, pair in zip(recipe.type, recipe.atom_idx):
        if type == ConversionType.BREAK:
            topology = break_bond_top(topology, pair)
        elif type == ConversionType.MOVE:
            topology = move_bond_top(topology, pair)

    write_topol(topology, newtop)


def break_bond_top(topology: Topology, breakpair: tuple[int, int]) -> Topology:
    """Break bonds in topology.
    removes bond, angles and dihedrals where breakpair was involved
    """
    topology["bonds"] = [
        bond
        for bond in topology["bonds"]
        if not bond[0] in breakpair and not bond[1] in breakpair
    ]
    topology["pairs"] = [
        pair
        for pair in topology["pairs"]
        if not pair[0] in breakpair and not pair[1] in breakpair
    ]
    topology["angles"] = [
        angle
        for angle in topology["angles"]
        if not angle[0] in breakpair and not angle[1] in breakpair
    ]
    topology["dihedrals"] = [
        dihedral
        for dihedral in topology["dihedrals"]
        if not dihedral[1] in breakpair and not dihedral[2] in breakpair
    ]
    return topology


def move_bond_top(topology: Topology, movepair: tuple[int, int]) -> Topology:
    raise NotImplementedError(
        "Topology Changer for moving Atoms is not implemented yet."
    )


def break_bond_plumed(plumeddat, breakpair, newplumeddist):
    new_distances = []
    broken_distances = []
    for line in plumeddat["distances"]:
        if breakpair[0] in line["atoms"] or breakpair[1] in line["atoms"]:
            broken_distances.append(line["id"])
        else:
            new_distances.append(line)

    plumeddat["distances"] = new_distances

    for line in plumeddat["prints"]:
        line["ARG"] = [id for id in line["ARG"] if not id in broken_distances]
        line["FIlE"] = newplumeddist

    return plumeddat


def move_bond_plumed(plumeddat, movepair, newplumeddist):
    raise NotImplementedError(
        "Plumeddat Changer for moving Atoms is not implemented yet."
    )


def modify_plumed(
    recipe: ConversionRecipe,
    oldplumeddat: Path,
    newplumeddat: Path,
    newplumeddist: Path,
):

    logging.info(
        f"Reading: {oldplumeddat} and writing modified plumed input to {newplumeddat}. Also writing {newplumeddist}."
    )
    plumeddat = read_plumed(oldplumeddat)

    for type, pair in zip(recipe.type, recipe.atom_idx):
        if type == ConversionType.BREAK:
            plumeddat = break_bond_plumed(plumeddat, pair, newplumeddist)
        elif type == ConversionType.MOVE:
            plumeddat = move_bond_plumed(plumeddat, pair, newplumeddist)

    write_plumed(plumeddat, newplumeddat)
