import logging
from kimmdy.reaction import (
    Break,
    Recipe,
    RecipeCollection,
    ReactionPlugin,
)
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    morse_transition_rate,
    get_atominfo_from_plumedid,
    get_bondprm_from_atomtypes,
)
from kimmdy.parsing import (
    read_topol,
    read_rtp,
    read_plumed,
    read_distances_dat,
    read_edissoc,
)
import MDAnalysis as mda
from pathlib import Path


class Homolysis(ReactionPlugin):
    """Homolytic bond breaking leading to 2 radicals.
    Implementation for time-varying rates
    """

    type_scheme = {"homolysis": {"dat": Path, "itp": Path}}

    def get_recipe_collection(self, files: TaskFiles):
        logging.debug("Getting recipe for reaction: homolysis")

        # Initialization of filepaths
        plumed_dat = files.input["plumed.dat"]
        distances_dat = files.input["distances.dat"]
        topol_top = files.input["top"]
        ffbonded_itp = files.input["ffbonded.itp"]
        edissoc_dat = files.input["edissoc.dat"]

        # Initialization of objects from files
        distances = read_distances_dat(distances_dat)
        plumed = read_plumed(plumed_dat)
        top = read_topol(topol_top)
        ffbonded = read_rtp(ffbonded_itp)
        edissoc = read_edissoc(edissoc_dat)

        #
        recipes = []
        for plumedid, dists in distances.items():
            if plumedid == "time":
                continue
            # get from plumedid to b0 and kb of the bond via atomtypes
            atomtypes, atomids = get_atominfo_from_plumedid(plumedid, plumed, top)
            b0, kb, E_dis = get_bondprm_from_atomtypes(atomtypes, ffbonded, edissoc)

            logging.debug(
                f"plumedid: {plumedid}, atomids: {atomids}, atomtypes: {atomtypes}, b0: {b0}, kb: {kb}, E_dis: {E_dis}"
            )

            k_avg, _ = morse_transition_rate([sum(dists) / len(dists)], b0, E_dis, kb)
            # averaging distances works here because we typically have
            # one conformational state per calculation

            # converto to zero-base
            atomids = [i - 1 for i in list(atomids)]

            recipes.append(
                Recipe(
                    recipe_steps=[Break(*list(atomids))],
                    rates=[*k_avg],
                    timespans=[[distances["time"][0], distances["time"][-1]]],
                )
            )

        return RecipeCollection(recipes)
