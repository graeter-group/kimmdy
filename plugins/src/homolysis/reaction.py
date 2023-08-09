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
    read_top,
    read_rtp,
    read_plumed,
    read_distances_dat,
    read_edissoc,
)


class Homolysis(ReactionPlugin):
    """Homolytic bond breaking leading to 2 radicals.
    Implementation for time-varying rates
    """

    def get_recipe_collection(self, files: TaskFiles):
        logging.debug("Getting recipe for reaction: homolysis")

        # Initialization of filepaths
        plumed_dat = files.input["plumed.dat"]
        distances_dat = files.input["distances.dat"]
        topol_top = files.input["top"]
        ffbonded_itp = self.config.itp
        edissoc_dat = self.config.edis

        # Initialization of objects from files
        distances = read_distances_dat(distances_dat)
        plumed = read_plumed(plumed_dat)
        top = read_top(topol_top)
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

            recipes.append(
                Recipe(
                    recipe_steps=[Break(atom_id_1 = atomids[0], atom_id_2 = atomids[1])],
                    rates=[*k_avg],
                    timespans=[(distances["time"][0], distances["time"][-1])],
                )
            )

        return RecipeCollection(recipes)
