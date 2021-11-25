from __future__ import annotations
import logging
import queue
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from typing import Callable
from kimmdy.utils import run_shell_cmd
from kimmdy.config import Config
from kimmdy.reaction import ConversionType, ConversionRecipe
#from kimmdy.runmanager import RunManager

import subprocess as sp

class MDManager:
    def __init__(self,top,struct,idx,mdppath,iteration,basename,dryrun):
        self.topfile = top
        self.struct = struct
        self.indexfile = idx
        self.mdpfile = mdppath
        self.it = str(iteration)

        self.outgro = basename + self.it + ".gro"
        self.outtpr = basename + self.it + ".tpr"
        self.outtrr= basename + self.it + ".trr"

        self.dryrun = dryrun
        pass

    def write_mdp(self):
        pass

    def minimzation(self):
        if self.dryrun:
            logging.info("Pretending to run minimization")
            return

        run_shell_cmd(
            f"gmx grompp -p {self.topfile} -c {self.grofile} -r {self.grofile} -f {self.mdpfile} -o {self.tprfile}"
        )
        run_shell_cmd(f"gmx mdrun -s {self.tprfile} -c {self.outgro}")

    def equilibration(self,ensemble):

        if self.dryrun:
            logging.info("Pretending to run {} equilibration".format(ensemble))
            return

        run_shell_cmd(f"gmx grompp -p {self.topfile} -c {self.grofile} -r {self.grofile} -f {self.mdpfile} -o {self.tprfile}")
        run_shell_cmd(f"gmx mdrun -v -s {self.tprfile} -c {self.outgro}")

    def equilibrium(self):
        # equilibrium before pulling md
         run_shell_cmd(f"gmx grompp -p {self.topfile} -c {self.struct} -f {self.mdpfile} -n {self.indexfile} -o {self.outtpr} -maxwarn 5")
         run_shell_cmd(f"gmx mdrun -v -s {self.outtpr} -c {self.outgro} -o {self.outtrr}") #use mpirun mdrun_mpi on cluster / several nodes are used
         return self.outgro


    def production(self,cpt,plumedfile):
        # normal pulling md
        run_shell_cmd(f" gmx grompp -p {self.topfile} -c {self.struct} -f {self.mdpfile} -n {self.indexfile} -t {cpt} -o {self.outtpr} -maxwarn 5")
        run_shell_cmd(f"gmx mdrun -v -s {self.outtpr} -c {self.outgro} -plumed {plumedfile} -o {self.outtrr}")  # use mpirun mdrun_mpi on cluster
        return self.outgro


    def relaxation(self,cpt):
        ## equil after break -->> called from changer
        run_shell_cmd(f" gmx grompp -p {self.topfile} -c {self.struct} -f {self.mdpfile} -n {self.indexfile} -t {cpt} -o {self.outtpr} -maxwarn 5")
        run_shell_cmd(f"gmx mdrun -v -s {self.outtpr} -c {self.outgro} -o {self.outtrr}") #use mpirun mdrun_mpi on cluster / several nodes are used
        return self.outgro

