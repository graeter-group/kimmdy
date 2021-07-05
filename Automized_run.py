"""
Part of reactive Kinetic Monte Carlo / Molecular Dynamics Simulations (rKMC/MD)
This modules contains all functions that are directly executing task in the terminal or need to communicate with the terminal
"""

import numpy as np
import subprocess as sp
import os


def do_production_run(grofile, topfile, mdpfile, indexfile, tprfile, trrfile):

    run_shell_cmd(f'gmx grompp -p {topfile} -c {grofile} -r {grofile} -f {mdpfile} -n {indexfile} -o {tprfile} -maxwarn 5')
    run_shell_cmd(f'gmx mdrun -v -s {tprfile} -o {trrfile}') # use mpirun mdrun_mpi on cluster



def do_equilibration(grofile, topfile, mdpfile, tprfile, outgro):

    run_shell_cmd(f'gmx grompp -p {topfile} -c {grofile} -r {grofile} -f {mdpfile} -o {tprfile}')
    run_shell_cmd(f'gmx mdrun -v -s {tprfile} -c {outgro}')


def do_energy_minimisation(grofile, topfile, mdpfile, tprfile, outgro):

    run_shell_cmd(f'gmx grompp -p {topfile} -c {grofile} -r {grofile} -f {mdpfile} -o {tprfile}')
    run_shell_cmd(f'gmx mdrun -s {tprfile} -c {outgro}')


def energy_min_after_break(oldcpt, newmdp, oldtpr, newtpr, newtop, indexfile, outgro, newtrr):
    run_shell_cmd(f' gmx grompp -f {newmdp} -p {newtop} -n {indexfile} -c {oldtpr} -r {oldtpr} -o {newtpr} -t {oldcpt} -maxwarn 5')
    run_shell_cmd(f'gmx mdrun -v -s {newtpr} -c {outgro}  -o {newtrr}')  # use mpirun mdrun_mpi on cluster


def do_plumed_run(grofile, topfile, mdpfile, indexfile, tprfile, trrfile, plumedfile):

    run_shell_cmd(f'gmx grompp -p {topfile} -c {grofile} -r {grofile} -f {mdpfile} -n {indexfile} -o {tprfile} -maxwarn 70')
    # use gmx mdrun if not on cluster  #use mpirun mdrun_mpi on cluster
    run_shell_cmd(f'mpirun mdrun_mpi  -v -s {tprfile} -o {trrfile} -plumed {plumedfile}')


def continue_run(oldcpt, newmdp, oldtpr, newtpr, newtop, indexfile, newtrr, plumedfile):
    run_shell_cmd(f' gmx grompp -f {newmdp} -p {newtop} -n {indexfile} -c {oldtpr} -r {oldtpr} -o {newtpr} -t {oldcpt} -maxwarn 75')
    run_shell_cmd(f'gmx mdrun -maxh 22.0 -v -s {newtpr} -o {newtrr} -plumed {plumedfile}')  # use mpirun mdrun_mpi on cluster


def create_dir(dir_name):
    run_shell_cmd(f'mkdir -p {dir_name}')


def change_to_dir(dir_name):
    os.chdir(dir_name)


def rerun_trajectory(trrfile, tprfile, edrfile, topfile, grofile, mdpfile):

    run_shell_cmd(f'gmx grompp -p {topfile} -c {grofile} -f {mdpfile} -o {tprfile}')  # Warum nochmal?
    run_shell_cmd(f'gmx mdrun -rerun {trrfile} -s {tprfile} -v -e {edrfile}')


def run_shell_cmd(s):
    command = sp.Popen(s, shell=True)
    command.wait()
