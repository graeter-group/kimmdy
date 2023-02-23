from __future__ import annotations
import logging
from typing import Optional
from kimmdy.reaction import ConversionRecipe, ConversionType
from kimmdy.parsing import (
    read_topol,
    write_topol,
    write_plumed,
    read_plumed
)
from kimmdy.topology.topology import Topology
from pathlib import Path

def modify_top(
        recipe: ConversionRecipe, oldtop: Path, newtop: Path, ffdir: Path, ffpatch: Optional[Path], topology: Optional[Topology]
):
    logging.info(f"Reading: {oldtop} and writing modified topology to {newtop}.")
    if topology is None:
        topologyDict = read_topol(oldtop)
        topology = Topology(topologyDict, ffdir, ffpatch)

    for conversion in recipe:
        if conversion.type == ConversionType.BREAK:
            topology.break_bond(conversion.atom_idx)
        elif conversion.type == ConversionType.BIND:
            topology.bind_bond(conversion.atom_idx)
    topology._update_dict()
    write_topol(topology.top, newtop)


def modify_plumed(
    recipe: ConversionRecipe,
    oldplumeddat: Path,
    newplumeddat: Path,
    plumeddist: Path,
):
    logging.info(
        f"Reading: {oldplumeddat} and writing modified plumed input to {newplumeddat}."
    )
    plumeddat = read_plumed(oldplumeddat)

    for conversion in recipe:
        if conversion.type == ConversionType.BREAK:
            plumeddat = break_bond_plumed(plumeddat, conversion.atom_idx, plumeddist)

    # TODO: handle BIND / MOVE
    write_plumed(plumeddat, newplumeddat)


def break_bond_plumed(plumeddat, breakpair, plumeddist):
    new_distances = []
    broken_distances = []
    breakpair = [str(x) for x in breakpair]
    for line in plumeddat["distances"]:
        if all(x in line["atoms"] for x in breakpair):
            broken_distances.append(line["id"])
        else:
            new_distances.append(line)

    plumeddat["distances"] = new_distances

    for line in plumeddat["prints"]:
        line["ARG"] = [id for id in line["ARG"] if not id in broken_distances]
        line["FILE"] = plumeddist

    return plumeddat

