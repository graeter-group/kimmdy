from kimmdy.reaction import Reaction, ReactionResult, ConversionRecipe, ConversionType
from .HAT_utils import cap_single_rad, find_radicals
import logging
from pathlib import Path
import MDAnalysis as MDA
import numpy as np
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
        u = MDA.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

        logging.debug(u.atoms[:40].types)
        rads = find_radicals(u)
        logging.debug(f"{rads} for {tpr}")
        logging.debug([u.atoms[:20].elements, u.atoms[:20].types])

        RR = ReactionResult(recipes=[], rates=[])
        for rad in rads:
            if len(rad.atoms) == 0:
                logging.info("no radical found, returning zero rate recipe")
                return ReactionResult(recipes=[], rates=[0])
            bonded_rad = rad[0].bonded_atoms
            logging.warning([rad, bonded_rad])

            subsystems = cap_single_rad(
                u, u.trajectory[-2], rad, bonded_rad, h_cutoff=3.5
            )

            for subsystem in subsystems:
                from_H = subsystem["meta"]["indices"][0]
                from_H_nr = str(u.atoms[from_H].index + 1)
                logging.warning(u.atoms[from_H].resname)
                if u.atoms[from_H].resname not in [
                    "NME",
                    "ACE",
                ]:  # doesn't work with capping groups at the moment
                    rad_nr = str(rad.atoms[0].index + 1)
                    CR = ConversionRecipe(
                        type=[ConversionType.MOVE], atom_idx=[[from_H_nr, rad_nr]]
                    )
                    RR.recipes.append(CR)
                    RR.rates.append(get_reaction_rates())

        logging.warning(f"Returning exactly these recipes to runmanager: {RR}")
        return RR
