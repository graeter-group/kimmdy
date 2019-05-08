"""
Part of reactive Kinetic Monte Carlo / Molecular Dynamics Simulations (rKMC/MD)
This modules contains all functions that are directly executing task in the terminal or need to communicate with the terminal
"""


import numpy as np
import subprocess as sp
import os

def do_production_run(grofile,topfile,mdpfile, indexfile, tprfile,trrfile):

    command=sp.Popen('gmx grompp -p '+topfile+' -c '+grofile+' -r '+grofile+' -f '+mdpfile + ' -n ' +indexfile +  ' -o '+tprfile + ' -maxwarn 5' ,shell=True)
    command.wait()
    command=sp.Popen('gmx mdrun -v -s '+tprfile+ ' -o ' +trrfile,shell=True) #use mpirun mdrun_mpi on cluster / several nodes are used
    command.wait()
    
    

def do_equilibration(grofile,topfile,mdpfile,tprfile,outgro):

    command=sp.Popen('gmx grompp -p '+topfile+' -c '+grofile+' -r '+grofile+' -f '+mdpfile+' -o '+tprfile,shell=True)
    command.wait()
    command=sp.Popen('gmx mdrun -v -s '+tprfile+' -c '+outgro,shell=True)
    command.wait()


def do_energy_minimisation(grofile,topfile,mdpfile,tprfile,outgro):

    command=sp.Popen('gmx grompp -p '+topfile+' -c '+grofile+' -r '+grofile+' -f '+mdpfile+' -o '+tprfile,shell=True)
    command.wait()
    command=sp.Popen('gmx mdrun -s '+tprfile+' -c '+outgro,shell=True)
    command.wait()

def equilibration_after_break(oldcpt, newmdp, oldtpr, newtpr, newtop, indexfile, outgro, newtrr):
    command = sp.Popen(" gmx grompp -f " + newmdp + " -p " + newtop + ' -n ' + indexfile + " -c "+ oldtpr + " -r "+oldtpr+" -o " + newtpr + " -t " + oldcpt + ' -maxwarn 5', shell = True)
    command.wait()
    command=sp.Popen('gmx mdrun -v -s '+newtpr+ ' -c '+ outgro + ' -o ' + newtrr,shell=True)  # use mpirun mdrun_mpi on cluster
    command.wait()

def continue_run (oldcpt, newmdp, oldtpr, newtpr, newtop, indexfile, newtrr, plumedfile):
    command = sp.Popen(" gmx grompp -f " + newmdp + " -p " + newtop + ' -n ' + indexfile + " -c "+ oldtpr + " -r "+oldtpr+" -o " + newtpr + " -t " + oldcpt + ' -maxwarn 5', shell = True)
    command.wait()
    command=sp.Popen('gmx mdrun -v -s '+newtpr+' -o '+ newtrr + ' -plumed ' + plumedfile,shell=True)  # use mpirun mdrun_mpi on cluster
    command.wait()
 

def concat_all_trr(nbr_of_files):
    #beachte Reihenfolge! Wird so einfach hintereinander gehaengt. t.b.d.
    command=sp.Popen('gmx trjcat -f *.trr -o combined.trr -cat  -settime' ,shell=True)
    command.wait()
    inputstr = "0"
    command=command.communicate(input=inputstr.encode()) # setze erstes Frame auf 0 und dann alle auf "C" 
    command.wait()
    for i in range(nbr_of_files -1):
        inputstr = "c"
        command=command.communicate(input=inputstr.encode())
        command.wait()

def create_dir(dir_name):
    command=sp.Popen('mkdir -p ' + dir_name ,shell=True)
    command.wait()

def change_to_dir(dir_name):
    os.chdir(dir_name)

            

if __name__ == "__main__":
    pass
    
