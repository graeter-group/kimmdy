from __future__ import annotations
import logging
import MDAnalysis as mda
from typing import Optional
from kimmdy.reaction import Recipe, Bind, Break, Move, RecipeStep
from kimmdy.parsing import read_top, write_top, write_plumed, read_plumed
import numpy as np
from kimmdy.tasks import TaskFiles
from kimmdy.topology.topology import Topology
from pathlib import Path


def modify_coords(
    recipe_steps: list[RecipeStep],
    files: TaskFiles,
):
    logging.debug(f"Entering modify_coords with recipe_steps {recipe_steps}")

    trr = files.input["trr"]
    tpr = files.input["tpr"]
    trr_out = files.outputdir / "coord_mod.trr"
    gro_out = files.outputdir / "coord_mod.gro"
    files.output["trr"] = trr_out
    files.output["gro"] = gro_out

    assert not trr_out.exists(), f"Error: Output trr exists! {trr_out}"
    assert not gro_out.exists(), f"Error: Output gro exists! {gro_out}"

    u = mda.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

    ttime = None
    to_move = []
    for step in recipe_steps:
        if isinstance(step, Move) and step.new_coords is not None:
            assert not step.is_one_based, "RecipeSteps should be zero based!"
            if ttime is None:
                ttime = step.new_coords[1]

            elif abs(ttime - step.new_coords[1]) > 1e-5:  # 0.01 fs
                m = (
                    f"Multiple coordinate changes requested at different times!"
                    "\nRecipeSteps:{recipe_steps}"
                )
                logging.error(m)
                raise ValueError(m)
            to_move.append(step.idx_to_move)

    np.unique(to_move, return_counts=True)

    for ts in u.trajectory:
        if abs(ts.time - ttime) > 1e-5:  # 0.01 fs
            continue
        for step in recipe_steps:
            if isinstance(step, Move) and step.new_coords is not None:
                atm_move = u.select_atoms(f"index {step.idx_to_move}")
                atm_move.position = step.new_coords[0]

        break
    else:
        raise LookupError(
            f"Did not find time {ttime} in trajectory "
            f"with length {u.trajectory[-1].time}"
        )

    u.atoms.write(trr_out)
    u.atoms.write(gro_out)

    files.output = {"trr": trr_out}
    files.output = {"gro": gro_out}

    logging.debug(
        f"Exit modify_coords, final coordinates written to {trr_out.parts[-2:]}"
    )


def modify_top(
    recipe_steps: list[RecipeStep],
    files: TaskFiles,
    ffpatch: Optional[Path],
    topology: Optional[Topology],
):
    files.output = {"top": files.outputdir / "topol_mod.top"}
    oldtop = (files.input["top"],)
    newtop = (files.output["top"],)
    ffdir = (files.input["ff"],)
    # why is this a tuple??

    logging.info(f"Reading: {oldtop} and writing modified topology to {newtop}.")
    if topology is None:
        topologyDict = read_top(oldtop)
        topology = Topology(topologyDict, ffdir, ffpatch)

    focus = set()
    for step in recipe_steps:
        step = step.one_based()
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

        else:
            raise NotImplementedError(f"RecipeStep {step} not implemented!")
    topology._update_dict()
    write_top(topology.top, newtop)

    topology.patch_parameters(list(focus))

    return


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
