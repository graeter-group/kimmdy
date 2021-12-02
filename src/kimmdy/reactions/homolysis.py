import logging
from kimmdy.reaction import Reaction, ConversionRecipe, ConversionType
from kimmdy.utils import identify_atomtypes


class Homolysis(Reaction):
    """Homolytic bond breaking leading to 2 radicals"""

    def get_rates(self):
        logging.info("Generating rates for reaction: homolysis")
        rates = [0, 1, 42]
        return rates

    def get_recipe(self):
        logging.info("Generating conversion recipe for reaction: homolysis")
        recipe = ConversionRecipe(ConversionType.BREAK, [(1, 2)])
        return recipe

    def rupturerates(self, plumed_dat, distances_dat, top, ffbonded_itp, edissoc_dat):
        # read out bond distances and atomtypes
        list_of_breakpairs_and_distances = self.find_distances(
            plumed_dat, distances_dat
        )
        dic_of_nbrs_to_atomtypes = identify_atomtypes(top)

        list_of_recipes = []
        list_of_rates = []

        logging.info("Parameters for calc_av_rate:")
        logging.info(list_of_breakpairs_and_distances[0][0:10])

        ##go through all possible breakpairs, calculate their rupture rates
        for j in range(len(list_of_breakpairs_and_distances)):

            # get parameters (distances, atomtypes etc) for current potential breakpair
            breakpair = [
                list_of_breakpairs_and_distances[j][0],
                list_of_breakpairs_and_distances[j][1],
            ]
            distances = list_of_breakpairs_and_distances[j][2:]

            atomtypes = []
            atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[0]])
            atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[1]])

            r_0, k_f = self.find_bond_param(
                atomtypes, ffbonded_itp
            )  # read out from gromacs force field
            E_dis = self.find_Edis(
                atomtypes, edissoc_dat
            )  # read out from gromacs force field

            # calculate rupture probabilties
            k = self.calc_av_rate(distances, float(r_0), float(E_dis), float(k_f))
            if j == 0:
                logging.info(E_dis)
                logging.info(r_0)
                logging.info(k_f)
                logging.info(k)

            list_of_rates.append(k)
            list_of_recipes.append(
                ConversionRecipe(ConversionType.BREAK, breakpair, atomtypes)
            )
        return (list_of_rates, list_of_recipes)
