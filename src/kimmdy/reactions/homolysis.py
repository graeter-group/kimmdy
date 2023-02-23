import logging
from kimmdy.reaction import Conversion, Reaction, ConversionRecipe, ConversionType, ReactionOutcome, ReactionResult
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    identify_atomtypes,
    find_distances,
    find_Edis,
    find_bond_param,
    calc_av_rate,
)
from pathlib import Path


class Homolysis(Reaction):
    """Homolytic bond breaking leading to 2 radicals.
    Implemented according to kimmdy 1.0
    """

    type_scheme = {"homolysis": {"edis": Path, "bonds": Path}}

    def get_reaction_result(self, files: TaskFiles):
        logging.debug("Getting recipe for reaction: homolysis")

        plumed_dat = files.input["plumed.dat"]
        distances_dat = files.input["distances.dat"]
        top = files.input["top"]
        files.input["ffbonded.itp"] = self.config.bonds
        files.input["edissoc.dat"] = self.config.edis
        ffbonded_itp = files.input["ffbonded.itp"]
        edissoc_dat = files.input["edissoc.dat"]

        # read out bond distances and atomtypes
        list_of_breakpairs_and_distances = find_distances(plumed_dat, distances_dat)
        dic_of_nbrs_to_atomtypes = identify_atomtypes(top)

        result = ReactionResult()

        # go through all possible breakpairs, calculate their rupture rates
        for i, _ in enumerate(list_of_breakpairs_and_distances):

            # get parameters (distances, atomtypes etc) for current potential breakpair
            breakpair = (
                list_of_breakpairs_and_distances[i][0],
                list_of_breakpairs_and_distances[i][1],
            )
            distances = list_of_breakpairs_and_distances[i][2:]

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

            outcome = ReactionOutcome(
                recipe=[Conversion(ConversionType.BREAK, breakpair)],
                rate=k
            )
            result.append(outcome)

        return result
