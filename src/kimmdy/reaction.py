from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto
import numpy as np
import logging

from kimmdy.utils import get_data_from_file


class ConversionType(Enum):
    BREAK = auto()
    MOVE = auto()


@dataclass
class ConversionRecipe:
    type: ConversionType
    atom_idx: list  # old definition made some problems
    atom_type: list


class Reaction(ABC):
    @abstractmethod
    def get_rates():
        pass

    @abstractmethod
    def get_recipe() -> ConversionRecipe:  # no clue what this is doing
        pass

    def find_bond_param(self, atomtypes, filepath):
        # reads bond parameters vom gromacs forcefield file #based on ff.bonded.itp from amber99sb-ildn.ff
        data_all, data_array = get_data_from_file(filepath)
        comb1 = [atomtypes[0], atomtypes[1]]
        comb2 = [atomtypes[1], atomtypes[0]]

        found_first_entry = (
            False  # only use first entry correspondnig to [ bondtypes ] section
        )
        for i in range(len(data_array)):
            if (
                list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2
            ):  # or data_array[i][:2] == (atomtype2, atomtype1)) :
                r_0 = data_array[i][3]
                k_f = data_array[i][4]
                found_first_entry = True
                break  # stop when first entry found

        if not found_first_entry:
            # print("Warning: No bond parameters found")
            # logging.warning('No bond parameters found')
            # r_0 = 0.145   #default value in same order of magnitude
            # k_f = 280000   #default value in same order of magnitude
            print("Found no entry for " + str(atomtypes))
            logging.error("Found no entry for " + str(atomtypes))

        # print ('Info: Considering ' + str(atomtypes) + ': r_0 = ' + str(r_0))
        # logging.info ('Considering' + str(atomtypes) + ': r_0 = ' + str(r_0))

        return r_0, k_f

    def find_Edis(self, atomtypes, filepath):
        # reads dissociation energy from edissoc.dat (gromacs)

        data_all, data_array = get_data_from_file(filepath)
        comb1 = [atomtypes[0], atomtypes[1]]
        comb2 = [atomtypes[1], atomtypes[0]]
        Edis = 0
        for i in range(len(data_array)):
            if (
                list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2
            ):  # or data_array[i][:2] == (atomtype2, atomtype1)) :
                Edis = data_array[i][2]
        if not Edis:
            # print("Info: Used simplified atom types to determine dissociaton energy for morse potential")
            # logging.debug("Used simplified atom types to determine dissociaton energy for morse potential")
            comb1 = [atomtypes[0][0], atomtypes[1][0]]
            comb2 = [atomtypes[1][0], atomtypes[0][0]]
            for i in range(len(data_array)):
                if (
                    list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2
                ):  # or data_array[i][:2] == (atomtype2, atomtype1)) :
                    Edis = data_array[i][2]
        # print ('E_dis = ' + str(Edis))
        # logging.debug('E_dis = ' + str(Edis))
        if not Edis:
            # print("Warning: No morse dissociation energy found. Used 350 as default value")
            # logging.warning("No morse dissociation energy found. Used 350 as default value")
            Edis = 350  # default value in same order of magnitude
        return Edis

    def find_distances(self, plumedfile, datafile):
        data_all, data_array = get_data_from_file(plumedfile)

        # get all pairs from plumedfile
        list_of_pairs_and_distances = []
        nbr_of_pairs = 0
        for i in range(len(data_all)):
            if "DISTANCE" in data_all[i]:
                if "broken" in data_all[i]:  # leave out already broken distances
                    continue
                split1 = data_all[i].split(",")
                atom2 = split1[-1][:-2]
                split2 = split1[0].split("=")
                atom1 = split2[-1]
                nbr_of_pairs += 1
                list_of_pairs_and_distances.append([str(atom1), str(atom2)])

        # get all distances from datafile
        # Note: Open file line by line i.o.t. to avoid memory error
        list_of_distances = []
        header = True  # used to skip header
        ctr = 0
        with open(datafile) as f:
            for line in f:
                line_split = []
                ctr += 1
                if header:
                    header = False
                    continue
                line_split = np.asarray(line.split())
                for k in range(
                    1, nbr_of_pairs + 1
                ):  # shift by +1 i.o.t. leave out time (first entry)
                    distance = float(line_split[k])
                    list_of_pairs_and_distances[k - 1].append(distance)
        nbr_of_data_points = ctr - 1

        # print ('Collected distances from ' + str(datafile) + ' for ' + str(nbr_of_pairs) + ' pairs with ' + str(nbr_of_data_points) + ' distances per pair.')
        # logging.info('Collected distances from ' + str(datafile) + ' for ' + str(nbr_of_pairs) + ' pairs with ' + str(nbr_of_data_points) + ' distances per pair.')

        return list_of_pairs_and_distances

    def calc_transition_rate(self, r_curr, r_0, E_dis, k_f):
        # parameters
        kT = 2.479  # k_B T at 310K #in Gromacs units kj *mol^-1
        # tau =  0.16    #unfitted / theoretical pre-exponential factor #h/kT = 0.16 ps from transition state theory
        k_0 = 0.288  # pre-exponential factor #from fitting averaged C_a - N data to gromacs data, see paper  #or: 1/2pi sqrt(k/m)

        # calculates energy barrier crossing rate [in ps]; barrier based on the model V = V_morse - F*X

        beta = np.sqrt(k_f / (2 * E_dis))  # [beta] =1/nm since beta = sqrt(k/2D)

        # calculate inflection point corresponding to point with maximal force
        r_infl = (beta * r_0 + np.log(2)) / beta

        # calculate current force in bond F = -del V / del x
        if r_curr > r_infl:
            # logging.debug('Used maximum force for bond Evans model since position behind inflection point found')
            F = (
                2
                * beta
                * E_dis
                * np.exp(-beta * (r_infl - r_0))
                * (1 - np.exp(-beta * (r_infl - r_0)))
            )
        else:
            F = (
                2
                * beta
                * E_dis
                * np.exp(-beta * (r_curr - r_0))
                * (1 - np.exp(-beta * (r_curr - r_0)))
            )
        # logging.debug('Calculated force in bond F = ' + str(F))

        # calculate extrema of shifted potential i.o.t. get barrier hight
        rmin = r_0 - 1 / beta * np.log(
            (beta * E_dis + np.sqrt(beta ** 2 * E_dis ** 2 - 2 * E_dis * beta * F))
            / (2 * beta * E_dis)
        )
        rmax = r_0 - 1 / beta * np.log(
            (beta * E_dis - np.sqrt(beta ** 2 * E_dis ** 2 - 2 * E_dis * beta * F))
            / (2 * beta * E_dis)
        )

        Vmax = E_dis * (1 - np.exp(-beta * (rmax - r_0))) ** 2 - F * (
            rmax - r_0
        )  # Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero
        Vmin = E_dis * (1 - np.exp(-beta * (rmin - r_0))) ** 2 - F * (rmin - r_0)

        delta_V = Vmax - Vmin

        k = k_0 * np.exp(-delta_V / kT)  # [1/ps]

        if (
            float(r_curr) > rmax
        ):  # already jumped over the barrier? Even if not "open" in gromacs morse potential?
            pass
        if F <= 0.0:  # negative force: Vmax -> infinity impliying k -> 0
            k = 0.0
            # logging.info('Found negative force, most likely due to compression. Rate replaced with zero.')

        return k, F  # [0,1]

    def calc_av_rate(self, distances, r_0, E_dis, k_f):
        # average distances first, if necessary
        dist = []
        if len(distances) > 1:
            r_av = sum(distances[2:]) / len(distances[2:])
        else:
            r_av = float(distances[0])
            print(r_av)
        k, F = self.calc_transition_rate(r_av, r_0, E_dis, k_f)
        print(r_av, k, F)

        return k
