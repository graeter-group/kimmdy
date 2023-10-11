from kimmdy.recipe import (
    Break,
    Relax,
    Recipe,
    RecipeCollection,
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    morse_transition_rate,
    get_atomnrs_from_plumedid,
    get_atominfo_from_atomnrs,
    get_bondprm_from_atomtypes,
    get_edissoc_from_atomnames,
)
from kimmdy.parsing import (
    read_top,
    read_plumed,
    read_distances_dat,
    read_edissoc,
)


class Homolysis(ReactionPlugin):
    """Homolytic bond breaking leading to 2 radicals.
    Implementation for time-varying rates
    """

    def get_recipe_collection(self, files: TaskFiles):
        logger = files.logger
        logger.debug("Getting recipe for reaction: homolysis")

        # Initialization of filepaths
        files.input["itp"] = self.config.itp
        files.input["edis"] = self.config.edis

        # Initialization of objects from files
        distances = read_distances_dat(files.input["plumed_out"])
        plumed = read_plumed(files.input["plumed"])
        top = self.runmng.top
        ffbonded = read_top(files.input["itp"])
        edissoc = read_edissoc(files.input["edis"])

        #
        recipes = []
        for plumedid, dists in distances.items():
            if plumedid == "time":
                continue
            # get from plumedid to b0 and kb of the bond via atomtypes
            atomnrs = get_atomnrs_from_plumedid(plumedid, plumed)
            atomtypes, atomnames = get_atominfo_from_atomnrs(atomnrs, top)
            b0, kb = get_bondprm_from_atomtypes(atomtypes, ffbonded)
            E_dis = get_edissoc_from_atomnames(atomnames, edissoc)

            # logger.debug(
            #     f"plumedid: {plumedid}, atomids: {atomnrs}, atomtypes: {atomtypes}, b0: {b0}, kb: {kb}, E_dis: {E_dis}"
            # )

            k_avg, _ = morse_transition_rate([sum(dists) / len(dists)], b0, E_dis, kb)
            # averaging distances works here because we typically have
            # one conformational state per calculation

            recipes.append(
                Recipe(
                    recipe_steps=[
                        Break(atom_id_1=atomnrs[0], atom_id_2=atomnrs[1]),
                        Relax(),
                    ],
                    rates=[*k_avg],
                    timespans=[(distances["time"][0], distances["time"][-1])],
                )
            )

        return RecipeCollection(recipes)
