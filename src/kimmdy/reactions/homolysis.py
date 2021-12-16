import logging
from kimmdy.reaction import Reaction, ConversionRecipe, ConversionType, ReactionResult
from kimmdy.utils import (
    identify_atomtypes,
    find_distances,
    find_Edis,
    find_bond_param,
    calc_av_rate,
)


class Homolysis(Reaction):
    """Homolytic bond breaking leading to 2 radicals.
    Implemented according to kimmdy 1.0
    """

    def get_reaction_result(
        self, plumed_dat, distances_dat, top, ffbonded_itp, edissoc_dat
    ):
        logging.info("Getting recipe for reaction: homolysis")

        # read out bond distances and atomtypes
        list_of_breakpairs_and_distances = find_distances(plumed_dat, distances_dat)
        dic_of_nbrs_to_atomtypes = identify_atomtypes(top)

        result = ReactionResult(recipes=[], rates=[])

        logging.info("Parameters for calc_av_rate:")
        logging.info(list_of_breakpairs_and_distances[0][0:10])

        # go through all possible breakpairs, calculate their rupture rates
        for j in range(len(list_of_breakpairs_and_distances)):

            # get parameters (distances, atomtypes etc) for current potential breakpair
            breakpair = (
                list_of_breakpairs_and_distances[j][0],
                list_of_breakpairs_and_distances[j][1],
            )
            distances = list_of_breakpairs_and_distances[j][2:]

            atomtypes = []
            atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[0]])
            atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[1]])

            r_0, k_f = find_bond_param(
                atomtypes, ffbonded_itp
            )  # read out from gromacs force field
            E_dis = find_Edis(
                atomtypes, edissoc_dat
            )  # read out from gromacs force field

            # calculate rupture probabilties
            k = calc_av_rate(distances, float(r_0), float(E_dis), float(k_f))
            if j == 0:
                logging.info(E_dis)
                logging.info(r_0)
                logging.info(k_f)
                logging.info(k)

            result.rates.append(k)
            result.recipes.append(
                ConversionRecipe(type=[ConversionType.BREAK], atom_idx=[breakpair])
            )

        return result
