from __future__ import annotations
import logging
from typing import Optional
from kimmdy.reaction import Recipe, Bind, Break, Move, RecipeStep
from kimmdy.parsing import read_top, write_top, write_plumed, read_plumed
from kimmdy.topology.topology import Topology
from pathlib import Path


def modify_top(
    recipe_steps: list[RecipeStep],
    oldtop: Path,
    newtop: Path,
    ffdir: Path,
    ffpatch: Optional[Path],
    topology: Optional[Topology],
) -> Topology:
    logging.info(f"Reading: {oldtop} and writing modified topology to {newtop}.")
    if topology is None:
        topologyDict = read_top(oldtop)
        topology = Topology(topologyDict, ffdir, ffpatch)

    focus = set()
    for step in recipe_steps:
        if isinstance(step, Break):
            topology.break_bond([str(step.atom_idx_1), str(step.atom_idx_2)])
            focus.add(str(step.atom_idx_1))
            focus.add(str(step.atom_idx_2))
        elif isinstance(step, Bind):
            topology.bind_bond([step.atom_idx_1, step.atom_idx_2])
            focus.add(str(step.atom_idx_1))
            focus.add(str(step.atom_idx_2))
        elif isinstance(step, Move):
            top_done = False
            if step.idx_to_bind is not None and step.idx_to_break is None:
                # implicit H-bond breaking
                topology.move_hydrogen([step.idx_to_move, step.idx_to_bind])
                focus.add(str(step.idx_to_move))
                focus.add(str(step.idx_to_bind))
                top_done = True
            if step.idx_to_bind is not None and not top_done:
                topology.bind_bond([step.idx_to_move, step.idx_to_bind])
                focus.add(str(step.idx_to_move))
                focus.add(str(step.idx_to_bind))
            if step.idx_to_break is not None and not top_done:
                topology.break_bond([step.idx_to_move, step.idx_to_break])
                focus.add(str(step.idx_to_move))
                focus.add(str(step.idx_to_break))
            if step.new_coords is not None:
                raise NotImplementedError("Changing coordinates not implemented!")

        else:
            raise NotImplementedError(f"RecipeStep {step} not implemented!")
    topology._update_dict()
    write_top(topology.top, newtop)

    topology.patch_parameters(list(focus))

    return topology


def modify_plumed(
    recipe_steps: list[RecipeStep],
    oldplumeddat: Path,
    newplumeddat: Path,
    plumeddist: Path,
):
    logging.info(
        f"Reading: {oldplumeddat} and writing modified plumed input to {newplumeddat}."
    )
    plumeddat = read_plumed(oldplumeddat)

    for step in recipe_steps:
        if isinstance(step, Break):
            plumeddat = break_bond_plumed(
                plumeddat, (step.atom_idx_1, step.atom_idx_2), plumeddist
            )
        else:
            # TODO: handle BIND / MOVE
            logging.WARNING(f"Plumed changes for {step} not implemented!")

    write_plumed(plumeddat, newplumeddat)


def break_bond_plumed(plumeddat, breakpair: list[int, int], plumeddist):
    new_distances = []
    broken_distances = []
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
