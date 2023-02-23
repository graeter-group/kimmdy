from __future__ import annotations
import logging
from typing import Optional
from kimmdy.reaction import ConversionRecipe, ConversionType
from kimmdy.parsing import (
    read_topol,
    write_topol,
)
from kimmdy.topology.topology import Topology
from pathlib import Path

def modify_top(
        recipe: ConversionRecipe, oldtop: Path, newtop: Path, ffdir: Path, ffpatch: Path, topology: Optional[Topology]
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

