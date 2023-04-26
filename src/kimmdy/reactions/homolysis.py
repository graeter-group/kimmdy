import logging
from kimmdy.reaction import (
    Conversion,
    Reaction,
    ConversionRecipe,
    ConversionType,
    ReactionOutcome,
    ReactionResult,
)
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    identify_atomtypes,
    find_distances,
    find_Edis,
    find_bond_param,
    calc_av_rate,
    morse_transition_rate
)
from kimmdy.parsing import read_topol, read_rtp, read_plumed, read_distances_dat, read_edissoc
from pathlib import Path


class Homolysis(Reaction):
    """Homolytic bond breaking leading to 2 radicals.
    New implementation for time-varying rates
    """

    type_scheme = {"homolysis": {"edis": Path, "bonds": Path}}

    def get_reaction_result(self, files: TaskFiles):
        logging.debug("Getting recipe for reaction: homolysis")

        plumed_dat = files.input["plumed.dat"]
        distances_dat = files.input["distances.dat"]
        topol_top = files.input["top"]
        files.input["ffbonded.itp"] = self.config.bonds
        files.input["edissoc.dat"] = self.config.edis
        # shouldn't config files automatically be in the TaskFiles?
        # no, only does self.latest_files as defined in runmgr
        ffbonded_itp = files.input["ffbonded.itp"]
        edissoc_dat = files.input["edissoc.dat"]
       
        distances = read_distances_dat(distances_dat)
        plumed = read_plumed(plumed_dat)
        lookup_atomid_plumedid = {entry['id']:frozenset(entry['atoms']) for entry in plumed['distances']}
        top = read_topol(topol_top)
        lookup_atomtype_atomid = {atom[0]:atom[1] for atom in top['atoms']}


        ffbonded_dict = read_rtp(ffbonded_itp)
        lookup_ffbonded_atomtype = {frozenset(l[:2]):[float(l[3]),float(l[4])] for l in ffbonded_dict['bondtypes']['other']}
        lookup_edissoc_atomtype = read_edissoc(edissoc_dat)


        ts = distances['time']
        result = ReactionResult()
        for plumedid, dists in distances.items():
            if plumedid == 'time':
                continue
            atomids = lookup_atomid_plumedid[plumedid]
            atomids_list = list(atomids)
            atomtypes_list = [lookup_atomtype_atomid[atomids_list[0]],lookup_atomtype_atomid[atomids_list[1]]]
            atomtypes = frozenset(atomtypes_list)
            atomelements_list = [x[0] for x in atomtypes_list]
            for comb in [atomtypes_list,[atomtypes_list[0],atomelements_list[1]],[atomelements_list[0],atomtypes_list[1]],atomelements_list]:
                if (comb_set := frozenset(comb)) in lookup_edissoc_atomtype.keys():
                    E_dis = lookup_edissoc_atomtype[comb_set]
                    break
            b0, kb = lookup_ffbonded_atomtype[atomtypes]
            
            print(plumedid,atomids,atomtypes,b0,kb,E_dis)

            k_reaction = morse_transition_rate(dists,b0,E_dis,kb)
            k_avg = calc_av_rate(dists, b0, E_dis, kb)

            outcome = ReactionOutcome(
                recipe=[Conversion(ConversionType.BREAK, atomids_list)], rate=k_avg, r_ts=k_reaction, ts=ts
            )
            result.append(outcome)   
            if len(result) == 1:
                return result    

        return result











        # # read out bond distances and atomtypes
        # list_of_breakpairs_and_distances = find_distances(plumed_dat, distances_dat)
        # dic_of_nbrs_to_atomtypes = identify_atomtypes(top)

        # result = ReactionResult()

        # # go through all possible breakpairs, calculate their rupture rates
        # for i, _ in enumerate(list_of_breakpairs_and_distances):
        #     # get parameters (distances, atomtypes etc) for current potential breakpair
        #     breakpair = (
        #         list_of_breakpairs_and_distances[i][0],
        #         list_of_breakpairs_and_distances[i][1],
        #     )
        #     distances = list_of_breakpairs_and_distances[i][2:]

        #     atomtypes = []
        #     atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[0]])
        #     atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[1]])

        #     r_0, k_f = find_bond_param(
        #         atomtypes, ffbonded_itp
        #     )  # read out from gromacs force field
        #     E_dis = find_Edis(
        #         atomtypes, edissoc_dat
        #     )  # read out from gromacs force field

        #     # calculate rupture probabilties
        #     k = calc_av_rate(distances, float(r_0), float(E_dis), float(k_f))

        #     outcome = ReactionOutcome(
        #         recipe=[Conversion(ConversionType.BREAK, breakpair)], rate=k, r_ts=k, ts=0
        #     )
        #     result.append(outcome)

        #return result
