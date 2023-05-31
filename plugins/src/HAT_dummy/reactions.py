from kimmdy.reaction import (
    Break,
    Bind,
    Recipe,
    RecipeCollection,
    ReactionPlugin,
)
from .HAT_utils import cap_single_rad, find_radicals
import logging
import MDAnalysis as mda
from numpy.random import default_rng

rng = default_rng()


class HAT_reaction(ReactionPlugin):
    """HAT reaction"""

    def get_reaction_result(self, files) -> RecipeCollection:
        def get_reaction_rate():
            return rng.random()

        logging.info("Starting HAT reaction, will do cool things")

        tpr = files.input["tpr"]
        trr = files.input["trr"]
        u = mda.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

        logging.info(u.atoms[:40].types)
        rads = find_radicals(u)
        logging.debug(f"{rads} for {tpr}")
        logging.debug([u.atoms[:20].elements, u.atoms[:20].types])

        recipes = []
        for rad in rads:
            if len(rad.atoms) == 0:
                logging.debug("no radical found, returning empty RecipeCollection")
                return RecipeCollection([])

            bonded_rad = rad[0].bonded_atoms
            logging.info(f"radical: {rad}")
            logging.info(f"bonded_rad: {bonded_rad}")

            subsystems = cap_single_rad(
                u, u.trajectory[-1], rad, bonded_rad, h_cutoff=3.5
            )
            logging.info("made subsystem")
            for subsystem in subsystems:
                from_H = subsystem["meta"]["indices"][0]
                h_partner = u.atoms[from_H].bonded_atoms[0]
                from_H_nr = str(u.atoms[from_H].index + 1)
                h_partner_nr = str(h_partner.index + 1)
                logging.info(u.atoms[from_H].resname)
                if u.atoms[from_H].resname not in [
                    "NME",
                    "ACE",
                ]:  # doesn't work with capping groups at the moment
                    rad_nr = str(rad.atoms[0].index + 1)

                    rate = get_reaction_rate()
                    recipes.append(
                        Recipe(
                            [
                                Break(from_H_nr, h_partner_nr),
                                Bind(from_H_nr, rad_nr),
                            ],
                            rates=[rate],
                            timespans=[u.trajectory[0].time, u.trajectory[-1].time],
                        )
                    )
                    logging.info(f"Recipe: {recipes[-1]}, rate: {rate}")

        logging.info(f"Returning these recipes to runmanager: {recipes}")
        return RecipeCollection(recipes)
