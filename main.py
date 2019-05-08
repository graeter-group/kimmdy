### Python Pipeline for hybrid reactive Kinetic Monte Carlo / Molecular Dynamics Simulations (rKMC/MD) ###
### by Benedikt Rennekamp. Heidelberg Institute of Theoretical Studies (HITS)
### Current version: 1.0, lastly updated on 30/04/2019
### Please cite paper accordingly if used. 

"""
Overall function:
Enabling (Homolytic) Bond Scission in Tensed Proteins in all-atom MD simulations.

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


(Main) input parameter / files:

Name            Variable type       What                                            Default vale
dir_name        string              name of subdirectory that will be created       run_n with n: number of runs of rKMC/MD
                                    where files will be stored
topfile         string              (path to) topology input file                   ../topol.top
grofile         string              (path to) gromacs (.gro) input file             ../npt.gro
mdpfile         string              (path to) simulation parameter (.mdp) file      depends on the respective step that is conducted
filepath_bonds  string              (path to) the force field parameter file from   ../ffbonded.itp
                                    which the bond parameters will be taken 
filepath_edis   string              (path to) the file with dissociation energies   ../edissoc.dat
indexfile       string              (path to) the index file of the system          ../index_backbone.ndx
plumedfile      string              (path to) the plumed input file specifying      ../plumed.dat
                                    which distanes will be monitored

(Main) output files:

datafile        string              (path to) the output file where plumed writes   ../distances.dat 
                                    the monitored distances into
new_top         string              new topology after the bond rupture             ../broken_topol.top
logfile         string              log file with intermediate steps and results    ../log.log
                                    Note: loglevel can be adjusted (Default: INFO)
                                    to DEBUG to get more output
trrfile         string              trajectory of the simulation, cut into parts    str(counter) + 'run' + str(n) + '.trr'
                                    for the different steps of the scheme


"""



import functions as func      #Own Module: all functions except interaction with terminal / outside world
import Automized_run as auto  #Own Module: all functions interacting with the terminal / outside world
import logging


def run(n, max_nbr_of_stops, equil_and_min, kinetic = True, rejection = False):
    
    #start logging and create directory to save files
    print("#####Info: Start run_" + str(n)+'#####')
    dir_name = "run_" + str(n)
    auto.create_dir(dir_name)    #creates directory to store files
    auto.change_to_dir(dir_name)
    
    log = logging.getLogger()
    for hdlr in log.handlers[:]:  # remove old handlers i.o.t. create a new logfile when starting a new run
        log.removeHandler(hdlr)
    logging.basicConfig(filename='log.log', level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    logging.info("####---Started run" + str(n)+" of Reactive MD ---#####")

    #initialize some counters
    nbr_of_stops = 0
    counter_length = len(str(max_nbr_of_stops))
    counter = "{0:0={counter_length}d}".format(nbr_of_stops, counter_length=counter_length)

    ### energy minimization, nvt and npt equilibration, if needed. 
    if equil_and_min:
        topfile = '../topol.top'
        grofile = '../npt.gro'
        outgro = 'min' + str(n) + '.gro'
        mdpfile = '../min.mdp'
        tprfile = "em" + str(n) + ".tpr"    
        auto.do_energy_minimisation(grofile,topfile,mdpfile,tprfile,outgro) 

        grofile = outgro
        outgro = 'nvt' + str(n) + '.gro'
        tprfile = 'nvt' + str(n) + '.tpr'
        mdpfile = '../temp_equil.mdp'
        auto.do_equilibration(grofile,topfile,mdpfile,tprfile,outgro)

        grofile = outgro
        outgro = 'npt' + str(n) + '.gro'
        tprfile = 'npt' + str(n) + '.tpr'
        mdpfile = '../p_equil.mdp'
        auto.do_equilibration(grofile,topfile,mdpfile,tprfile,outgro)

    #Choose input options for (first) production run
    filepath_bonds = "../ffbonded.itp"
    filepath_edis = "../edissoc.dat"
        
    trrfile = str(counter) + 'run' + str(n) + '.trr'
    tprfile = str(counter) + 'run' + str(n) + '.tpr'
    edrfile = 'run' + str(n) + '.edr'
    topfile = '../topol.top'
    indexfile = '../index_backbone.ndx'
    grofile = '../npt.gro'  #use npt.gro for centered version
    if equil_and_min:
        grofile = outgro

    #Short run to equilbrate beforehand under the (lower or) same force such that this phase won't be used for the distance sampling
    mdpfile = "../pullf1500_equil.mdp" 
    auto.do_production_run(grofile,topfile,mdpfile, indexfile, tprfile,trrfile)

    #More input files for main run
    mdpfile = "../pullf1500.mdp"
    oldcpt = "state.cpt"
    plumedfile = "../plumed.dat"
    datafile = "distances.dat"
    oldtpr = tprfile
    nbr_of_stops += 1
    counter = "{0:0={counter_length}d}".format(nbr_of_stops, counter_length=counter_length)
    newtrr = str(counter) + 'run' + str(n) + '.trr'
    newtpr = str(counter) + 'run' + str(n) + '.tpr'
    
    #First main run
    auto.continue_run(oldcpt, mdpfile, oldtpr, newtpr, topfile, indexfile, newtrr, plumedfile)


    ### let it run until it stops then run again
    continuation = True
    ruptured_bonds = []
    first_break = True
    while continuation: 
        logging.info('------ Simulation paused --------')
        #func.del_backup_files_and_step_files ()  #clean up folder, enable if needed
        nbr_of_stops += 1
        counter = "{0:0={counter_length}d}".format(nbr_of_stops, counter_length=counter_length)

        #read out bond distances and atomtypes
        list_of_breakpairs_and_distances = func.find_distances(plumedfile, datafile)
        dic_of_nbrs_to_atomtypes = func.identify_atomtypes(topfile)

        list_of_nbrs_and_atomtypes = []
        list_of_rates = []
        
        ##go through all possible breakpairs, calculate their rupture rates
        for j in range(len(list_of_breakpairs_and_distances)):
            if rejection == True: #don't need all rates if rejection MC. Not used yet. 
                break
            
            #get parameters (distances, atomtypes etc) for current potential breakpair
            breakpair = [list_of_breakpairs_and_distances[j][0], list_of_breakpairs_and_distances[j][1]]
            distances = list_of_breakpairs_and_distances[j][2:]
            
            atomtypes = []
            atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[0]])
            atomtypes.append(dic_of_nbrs_to_atomtypes[breakpair[1]])
            
            r_0, k_f = func.find_bond_param(atomtypes, filepath_bonds) #read out from gromacs force field
            E_dis = func.find_Edis(atomtypes, filepath_edis)  #read out from gromacs force field
            
            #calculate rupture probabilties
            k = func.calc_av_rate(distances, float(r_0), float(E_dis), float(k_f))
            list_of_rates.append(k)
            list_of_nbrs_and_atomtypes.append([breakpair, atomtypes])
            
            print ('Info: For ' + str(breakpair) +' ' + str(atomtypes) + ' calculated rupture rate using ' + str((len(distances))) + ' distances per bond: ' + str(k) + ' per ps.')
            logging.info('For ' + str(breakpair) + ' ' +  str(atomtypes) + ' calculated rupture rate using ' + str((len(distances))) + ' distances per bond: ' + str(k) + ' per ps.')


        #Rejection-free kinetic MC scheme
        if kinetic == True and len(list_of_rates) > 0:
            print ("DO kinetic MC")
            logging.info("DO kinetic MC")
            breakpair, atomtypes, delta_t = func.do_kinetic_mc (list_of_nbrs_and_atomtypes, list_of_rates)
            rupture_time = delta_t / 1000.0 #convert ps to ns        
            if first_break:
                first_rupture = nbr_of_stops
                first_rupture_time = rupture_time
                first_break = False                     
            ruptured_bonds.append([breakpair, atomtypes, rupture_time])
            newtop = '../broken_topol.top'
            func.modify_top (topfile, newtop, breakpair)
            topfile = newtop
                
        #Rejection Kinetic Monte carlo Scheme. #Try-out implemenation, not used yet but left here for exploration.
        if rejection == True:
            #find one reaction candidate, calculate rate only for that one and do MC to see if accepted
            r_0 = 0.1 #[1/ps] #Estimated upper rate limit
            rupture, breakpair, atomtypes, delta_t = func.do_rejection_KMC(list_of_breakpairs_and_distances, r_0, topfile, filepath_bonds, filepath_edis)
            print ("---Info: RKMC Transition accepted? " + str(rupture) + "---")
            logging.info("--- RKMC Transition accepted? = " + str(rupture) + "---")
            if rupture:
                rupture_time = rupture_time = delta_t / 1000.0 #convert ps to ns 
                if first_break:
                    first_rupture = nbr_of_stops
                    first_rupture_time = rupture_time
                    first_break = False                     
                ruptured_bonds.append([breakpair, atomtypes, rupture_time])
                newtop = '../broken_topol.top'
                func.modify_top (topfile, newtop, breakpair)
                topfile = newtop
                
  
        logging.info('###End of sub-run evaluation: In this run the following bonds ruptured (or were already ruptured): ' + str(ruptured_bonds))
        print('###End of sub-run evaluation: Info: In this run the following bonds ruptured (or were already ruptured): ' + str(ruptured_bonds))


        ### Note: In principle, it is sufficient for one cycle of rKMC/MD to stop here.
        #In the following, a sample implemenation of the continued simulation is provided. However, this has to be adjusted to the system at hand.
        #In particular, the simulation time of the equilibration (here in '../broken_equil.mdp') and the second major MD run (here in '../pullf1000_broken.mdp') need to be long enough to reach a new equilibrium before a secondary rupture can reasonably be done.
        #If you like to continue manually your simulation, simply use the new topology ('broken_topol.top') that was created and your desired other input files.
        #See paper for more details!


        ###choose input parameters for next run
            
        #Decide what to do next
        equilrun_after_break = True #run with shorter timestep after the breackage to stabilize transition.
        modify_plumed = True   #use to modify your plumed-file for a secondary rupture. Excludes the already ruptured bond from the candidates.

        oldcpt = "state.cpt"
        oldtpr = tprfile
        
        #short equilibration after rupture. Uses a shorter timestep and harmonic bonds instead of Morse to stabilize transition   
        if equilrun_after_break:
            mdpfile = '../broken_equil_f1000.mdp'
            newtpr = str(counter) +'run' + str(n) + '_equil' + '.tpr'
            newtrr = str(counter) +'run' + str(n) + '_equil' + '.trr'
            outgro = str(counter) + 'run' + str(n) + '_equil' + '.gro'
            auto.equilibration_after_break(oldcpt, mdpfile, oldtpr, newtpr, topfile, indexfile, outgro, newtrr)
            oldtpr = outgro #continue afterwards with new gro from min (i.e. new ensemble!) instead of old tpr
            nbr_of_stops += 1
            counter = "{0:0={counter_length}d}".format(nbr_of_stops, counter_length=counter_length)

        #rupture happened: change plumed-inputfile etc.
        if modify_plumed:    
            newtpr = str(counter) +'run' + str(n) + '_rupture' + '.tpr'
            newtrr = str(counter) + 'run' + str(n) + '_ruptrue' + '.trr'
            plumedfile_new= str(counter) + 'plumed'+ str(n) + '.dat'
            distancefile_new = str(counter) + 'distances'+ str(n) + '.dat'
            func.modify_plumedfile(plumedfile, plumedfile_new, distancefile_new, ruptured_bonds)
            mdpfile = '../pullf1000_broken.mdp'
            logging.info("Since rupture event happend, will use the following mdp-file for next run: " + mdpfile)
            print("Info: Since rupture event happend, will use the following mdp-file for next run: " + mdpfile)
            newtpr = str(counter) +'run' + str(n) + '_broken' + '.tpr'
            newtrr = str(counter) + 'run' + str(n) + '_broken' + '.trr'
        else:        
        #continue with old mdp if no (further) rupture
            newtpr = str(counter) + 'run'+ str(n) + '_continue' + '.tpr'
            newtrr = str(counter) + 'run' + str(n) + '_continue'+ '.trr'


        #next subrun
        auto.continue_run(oldcpt, mdpfile, oldtpr, newtpr, topfile, indexfile, newtrr, plumedfile_new)

        #conditions to stop the loop: error occured or specified max number of runs reached
        print('Info: ### Sub-run with current number of stops = ' + str(nbr_of_stops) + ' finished. Starting again evalution of possible bond ruptures.###')
        logging.info('Sub-run with current number of stops = ' + str(nbr_of_stops) +' finished. Starting again evalution of possible bond ruptures.###')
        logfile = 'md.log'
        error = func.check_if_error_occured(logfile)   #bool #true if error occured
        if nbr_of_stops >= max_nbr_of_stops:
            continuation = False
        if error:
            continuation = False
        
    auto.change_to_dir('../')
       
    print ("Info: End")
    logging.info("###---Run of Reactive MD finished---###")



if __name__ == "__main__":
    max_nbr_of_stops = 1 #how many times the simulation may stop and continue again (mamximum)
    equil_and_min = False #Do energy min, NVT and NPT equilibration if true
    for i in range (4, 5):    #choose how many runs you like to conduct. Will create subfolders for each run, labeled by number: run_i 
        run(i, max_nbr_of_stops, equil_and_min)
