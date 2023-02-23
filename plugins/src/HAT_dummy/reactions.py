from kimmdy.reaction import Conversion, Reaction, ReactionOutcome, ReactionResult, ConversionRecipe, ConversionType
from .HAT_utils import cap_single_rad, find_radicals
import logging
import MDAnalysis as mda
from numpy.random import default_rng

rng = default_rng()


class HAT_reaction(Reaction):
    """HAT reaction"""

    def get_reaction_result(self, files) -> ReactionResult:
        def get_reaction_rates():
            return rng.random()

        logging.info("Starting HAT reaction, will do cool things")

        tpr = files.input["tpr"]
        trr = files.input["trr"]
        u = mda.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

        logging.info(u.atoms[:40].types)
        rads = find_radicals(u)
        logging.debug(f"{rads} for {tpr}")
        logging.debug([u.atoms[:20].elements, u.atoms[:20].types])

        outcomes = []
        for rad in rads:
            if len(rad.atoms) == 0:
                logging.info("no radical found, returning zero rate recipe")
                return [ReactionOutcome([], 0)]
            bonded_rad = rad[0].bonded_atoms
            logging.info(f'radical: {rad}')
            logging.info(f'bonded_rad: {bonded_rad}')

            subsystems = cap_single_rad(
                u, u.trajectory[-2], rad, bonded_rad, h_cutoff=3.5
            )
            logging.info('made subsystem')
            for subsystem in subsystems:
                from_H = subsystem["meta"]["indices"][0]
                from_H_nr = str(u.atoms[from_H].index + 1)
                logging.info(u.atoms[from_H].resname)
                if u.atoms[from_H].resname not in [
                    "NME",
                    "ACE",
                ]:  # doesn't work with capping groups at the moment
                    rad_nr = str(rad.atoms[0].index + 1)
                    recipe = [
                        # FIXME: which is the previous H binding partner?
                        # to break it's bond so that the H can move.
                        Conversion(ConversionType.BREAK, (from_H_nr, rad_nr)),
                        Conversion(ConversionType.BIND, (from_H_nr, rad_nr)),
                    ]
                    rate = get_reaction_rates()
                    outcomes.append(ReactionOutcome(recipe, rate))
                    logging.info(f'Made outcome with recipe: {recipe} and rate: {rate}')

        logging.info(f"Returning exactly these recipes to runmanager: {outcomes}")
        return outcomes

