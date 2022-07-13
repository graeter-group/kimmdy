from kimmdy.reaction import Reaction, ReactionResult, ConversionRecipe, ConversionType
from .HAT_utils import cap_single_rad
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


        # tpr = files['HAT']['tprpath']
        # trr = files['HAT']['trrpath']
        tpr = files.input["tpr"]
        trr = files.input["trr"]
        logging.info([str(tpr),str(trr)])
        u = MDA.Universe(str(tpr),str(trr),topology_format='tpr',format='trr')


        rad = u.select_atoms('resname ALA and name CA')
        bonded_rad = rad[0].bonded_atoms
        #print(rad)
        #print(bonded_rad)

        #rad[0].type = 'C'
        #print(rad.elements)
        #print(rad[0].element)

        subsystems = cap_single_rad(u,u.trajectory[10],rad,bonded_rad,h_cutoff=3)

        #print(subsystems)

        #print(rad.atoms[0].index +1,rad.atoms)
        RR = ReactionResult(recipes=[],rates=[])
        for subsystem in subsystems:
            from_H = subsystem['meta']['indices'][0]
            #print(u.atoms[from_H].index + 1,u.atoms[from_H])
            from_H_nr = u.atoms[from_H].index + 1
            rad_nr = rad.atoms[0].index +1
            CR = ConversionRecipe(type=ConversionType.MOVE,atom_idx=[from_H_nr,rad_nr])
            RR.recipes.append(CR)
            RR.rates.append(get_reaction_rates())

        logging.info(RR)
        return RR

    @property
    def type_scheme(self):
        """Dict of types of possible entries in config.
        Used to read and check the input config.
        To not use this feature return empty dict
        {"HAT":{'tprpath': Path, 'trrpath': Path}}
        """
        return dict()


#tprpath = Path("/hits/fast/mbm/hartmaec/workdir/HAT_reaction/md_Ala_delHA/md.tpr")
#trrpath = Path("/hits/fast/mbm/hartmaec/workdir/HAT_reaction/md_Ala_delHA/trjout.pdb")
# tprpath = Path('/hits/fast/mbm/hartmaec/kimmdys/kimmdy_topology/example/example_ala/test_out_012/production_2/prod.tpr')
# trrpath = Path('/hits/fast/mbm/hartmaec/kimmdys/kimmdy_topology/example/example_ala/test_out_012/production_2/prod.trr')

# files = {"HAT":{'tprpath':tprpath,'trrpath':trrpath}}

# HAT_reaction.get_reaction_result(None,files)
