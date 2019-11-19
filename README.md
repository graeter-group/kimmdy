# kimmdy
reactive MD pipeline for GROMACS using Kinetic Monte Carlo / Molecular Dynamics (KIMMDY)

This is a short readme / manual to KIMMDY 

### Python Pipeline for hybrid reactive Kinetic Monte Carlo / Molecular Dynamics Simulations (KIMMDY)
### by Benedikt Rennekamp. Heidelberg Institute of Theoretical Studies (HITS)
### Current version: 1.0, lastly updated on 31/03/2019
### Please cite paper accordingly if KIMMDY is used (https://doi.org/10.1021/acs.jctc.9b00786):  
Hybrid Kinetic Monte Carlo / Molecular Dynamics Simulations of Bond Scissions in Proteins

Benedikt Rennekamp, Fabian Kutzki, Agnieszka Obarska-Kosinska, Christopher Zapp, and Frauke Gräter

Journal of Chemical Theory and Computation

DOI: 10.1021/acs.jctc.9b00786 


Overall function: Enabling (Homolytic) Bond Scission in Tensed Proteins in all-atom MD simulations.

Software requirments:
- To be used with a Gromacs version that is patched with Plumed. Tested for Gromacs 2018.1 patched with Plumed 2.4.2 (see www.plumed.org)
- Tested with Python 2.7. Uses the modules: logging, numpy, random, subprocess, os

Scheme:

a) Do all-atom MD of your tensed protein (constant force) for a first, shorter sampling run. Default: 1ns simulation time (see paper for details)

During this time, the bond elongation of possbile bond rupture rates will be monitored.

b) Rate Calculation: Using these distances, for each bond a rupture rate based on Transition State Theory will be calculated.

c) Kinetic Monte Carlo: These rates will be used as input for a KMC step, determining
    i) which bond breaks (current implementation: rejection-free) and ii) the corresponding time step of that transition.

d) Adjustment of the topology accorinding to the break (removal of bonds, angles, pairs, dihedrals in topology)

e) Continuation of Simulation at that point (be aware of the system's time jump due to the Monte Carlo Step)


Short manual, how to use KIMMDY as provided in this folder:

1) Set-up regularly your pulling MD simulation (constant force) in Gromacs. There are example files provided for an individual collagen triple helix, that can be used for a quick exemplary test of KIMMDY. However, it is highly recommended to adjust these to your system at hand and consider general MD best practices before using KIMMDY.
2) Before the actual programm can be used, a "plumedfile" has to be generated that contains the distances that should be monitored during the main run, i.e. the bond elongation of all possible bond rupture candidates. An exemplary function to generate this file is provided in the "functions.py" module under the name "write_conditions_in_plumedfile(topfile, indexfile, indexgroup, parameterfile)". Based on the indexfile and the topology, it takes all C_a - N and C_a - C bonds in the backbone as candidates. As stated in the paper, this is a reasonable choice for collagen-based material.
The exemplary output of this function is provided in the file "plumed.dat". In this example, the bond distances will be printed every 50 steps ('stride frequency') to the file 'distances.dat'
3) Before starting the program, check the following points: 
The evalutation of the bond rupture rates is based on a (tilted) Morse potential for each bond. For the parameterization, the parameters from the force field will be read out. The current implementation is build on the amber99sb-star-ildnp force field and the file 'ffbonded.itp' from that force field is used as input. In principle, you can replace this with any other force field and input file, but than you might need to adjust the respective read-out function "find_bond_param(atomtypes, filepath)".
In a similar way, the dissociation energies are provided in an extra file 'edissoc.dat' that it provided by Gromacs and read out by 'func.find_Edis(atomtypes, filepath_edis)'. Please note that we adjusted the C_a - N bond energy to 348 as described in the paper. 
4) In the main.py-file, you can provide a few parameters (input file names, choosing if you like an automated equilibration, how many runs you like to conduct). See the commented code for more details. After choosing these parameters, you can start the program by executing main.py. 
For each run, a subfolder in the current working directory will be created and named by the number i of the run ('run_i'). 
5) The current implementation provides a basic logging functionality, storing essential steps and results in 'log.log' for each run. More information will be printed if you change the log-level to 'debug' in the main function.
6) As an output, KIMMDY will yield the bond that breaks after the cycle, the corresponding time evolution due to the Kinetic Monte Carlo Step and the adjusted topology 'broken_topol.top', from which the simulation can be continued (together with the state-file 'run_x.cpt').
Note that this new topology will be written to the original directory and not the newly created sub-directory of the individual run.
7) An exemplary automated continuation of the simulation is also provided in the code.

 
Feel free to message Benedikt Rennekamp (benedikt.rennekamp@h-its.org) or Prof. Dr. Frauke Graeter (frauke.graeter@h-its.org) for questions and suggestions!



