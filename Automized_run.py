"""
Part of reactive Kinetic Monte Carlo / Molecular Dynamics Simulations (rKMC/MD)
This modules contains all functions that are directly executing task in the terminal or need to communicate with the terminal
"""

import numpy as np
import subprocess as sp
import os


def do_production_run(grofile, topfile, mdpfile, indexfile, tprfile, trrfile):

    command = sp.Popen('gmx grompp -p '+topfile+' -c '+grofile+' -r '+grofile+' -f ' +
                       mdpfile + ' -n ' + indexfile + ' -o '+tprfile + ' -maxwarn 5', shell=True)
    command.wait()
    # use gmx mdrun if not on cluster  #use mpirun mdrun_mpi on cluster
    command = sp.Popen('mpirun mdrun_mpi -v -s '+tprfile +
                       ' -o ' + trrfile, shell=True)
    command.wait()


def do_equilibration(grofile, topfile, mdpfile, tprfile, outgro):

    command = sp.Popen('gmx grompp -p '+topfile+' -c '+grofile +
                       ' -r '+grofile+' -f '+mdpfile+' -o '+tprfile, shell=True)
    command.wait()
    command = sp.Popen('gmx mdrun -v -s '+tprfile+' -c '+outgro, shell=True)
    command.wait()


def do_energy_minimisation(grofile, topfile, mdpfile, tprfile, outgro):

    command = sp.Popen('gmx grompp -p '+topfile+' -c '+grofile +
                       ' -r '+grofile+' -f '+mdpfile+' -o '+tprfile, shell=True)
    command.wait()
    command = sp.Popen('gmx mdrun -s '+tprfile+' -c '+outgro, shell=True)
    command.wait()


def energy_min_after_break(oldcpt, newmdp, oldtpr, newtpr, newtop, indexfile, outgro, newtrr):
    command = sp.Popen(" gmx grompp -f " + newmdp + " -p " + newtop + ' -n ' + indexfile + " -c " +
                       oldtpr + " -r "+oldtpr+" -o " + newtpr + " -t " + oldcpt + ' -maxwarn 5', shell=True)
    command.wait()
    command = sp.Popen('mpirun mdrun_mpi -v -s '+newtpr + ' -c ' + outgro +
                       ' -o ' + newtrr, shell=True)  # use mpirun mdrun_mpi on cluster
    command.wait()


def do_plumed_run(grofile, topfile, mdpfile, indexfile, tprfile, trrfile, plumedfile):

    command = sp.Popen('gmx grompp -p '+topfile+' -c '+grofile+' -r '+grofile+' -f ' +
                       mdpfile + ' -n ' + indexfile + ' -o '+tprfile + ' -maxwarn 70', shell=True)
    command.wait()
    # use gmx mdrun if not on cluster  #use mpirun mdrun_mpi on cluster
    command = sp.Popen('mpirun mdrun_mpi  -v -s '+tprfile +
                       ' -o ' + trrfile + ' -plumed ' + plumedfile, shell=True)
    command.wait()


def continue_run(oldcpt, newmdp, oldtpr, newtpr, newtop, indexfile, newtrr, plumedfile):
    command = sp.Popen(" gmx grompp -f " + newmdp + " -p " + newtop + ' -n ' + indexfile + " -c " +
                       oldtpr + " -r "+oldtpr+" -o " + newtpr + " -t " + oldcpt + ' -maxwarn 75', shell=True)
    command.wait()
    command = sp.Popen('mpirun mdrun_mpi -maxh 22.0 -v -s '+newtpr+' -o ' + newtrr +
                       ' -plumed ' + plumedfile, shell=True)  # use mpirun mdrun_mpi on cluster
    command.wait()


def create_dir(dir_name):
    command = sp.Popen('mkdir -p ' + dir_name, shell=True)
    command.wait()


def change_to_dir(dir_name):
    os.chdir(dir_name)


def rerun_trajectory(trrfile, tprfile, edrfile, topfile, grofile, mdpfile):

    command = sp.Popen('gmx grompp -p '+topfile+' -c '+grofile +
                       ' -f '+mdpfile+' -o '+tprfile, shell=True)  # Warum nochmal?
    command.wait()
    command = sp.Popen('gmx mdrun -rerun '+trrfile+' -s ' +
                       tprfile+' -v -e '+edrfile, shell=True)
    command.wait()


if __name__ == "__main__":
    pass
