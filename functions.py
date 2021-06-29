"""
Part of reactive Kinetic Monte Carlo / Molecular Dynamics Simulations (KIMMDY)
This modules contains all python-internal functions, i.e. all except the ones that are directly executing tasks in the terminal or need to communicate with the outside world. 
"""

import numpy as np
import logging
import random
import pprint as pp
import MDAnalysis

def get_data_from_file(filepath):
    file = open(filepath, 'r')
    data_all = []  #array with each entry corresponding to one (string) line
    data_array = []  #array of all lines with each subarray containing one value
    settings = []
    for line in file:
        data_all.append (line)
        line_array = np.asarray(line.split())   
        data_array.append(line_array)
    file.close()

    return data_all, data_array

    
def store_linelist_to_file(data, filepath):
    file = open (filepath, "w")
    for line in data:
        file.write(line)
    file.close()


def identify_atomtypes(filepath):

    dic_of_atoms_to_groups = {}
    file = open(filepath, 'r')
    atoms = False
    
    for line in file:
        #start collecting data when [ atoms ] reached
        if 'atoms' in line:
            atoms = True
            continue
        if atoms == False: 
            continue
        
        #reached [ bonds ], stop
        if 'bonds' in line: 
            break
        
        #leave out comments, includes, and empty lines
        if len(line) <= 1:
            continue
        elif line[0] == ";": 
           continue
        elif line[0] == "#":
            continue
        
        #hopefully, these are all and exclusively the lines now with atomnbrs, types etc. 
        line_array = np.asarray(line.split())
        nbr = line_array[0]
        atomtype = line_array[1]
        dic_of_atoms_to_groups.update({nbr:atomtype})    
    file.close()

    return dic_of_atoms_to_groups          
 

def find_bond_param(atomtypes, filepath):
    #reads bond parameters vom gromacs forcefield file #based on ff.bonded.itp from amber99sb*-ildnp.ff
    data_all, data_array = get_data_from_file(filepath)
    comb1 = [atomtypes[0], atomtypes[1]]
    comb2 = [atomtypes[1], atomtypes[0]]     
    
    found_first_entry = False  #only use first entry correspondnig to [ bondtypes ] section
    for i in range(len(data_array)):
        if ( list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2) : #or data_array[i][:2] == (atomtype2, atomtype1)) :
            r_0 = data_array[i][3]
            k_f = data_array[i][4]
            found_first_entry = True
            break #stop when first entry found

    if not found_first_entry:
        print("Warning: No bond parameters found")
        logging.warning('No bond parameters found')
        r_0 = 0.145   #default value in same order of magnitude
        k_f = 280000   #default value in same order of magnitude

        if comb2_equil in dic_of_equil_r0 :
            r_0 = dic_of_equil_r0[comb2_equil]
        elif comb1_equil in dic_of_equil_r0:
            r_0 = dic_of_equil_r0[comb1_equil]
        else:
            print ('Found no entry for ' + str(atomtypes))
            logging.info('Found no entry for ' + str(atomtypes))

    print ('Info: Considering ' + str(atomtypes) + ': r_0 = ' + str(r_0))
    logging.info ('Considering' + str(atomtypes) + ': r_0 = ' + str(r_0))
    
    return r_0, k_f

def find_Edis(atomtypes, filepath):
    #reads dissociation energy from edissoc.dat (gromacs)

    data_all, data_array = get_data_from_file(filepath)
    comb1 = [atomtypes[0], atomtypes[1]]
    comb2 = [atomtypes[1], atomtypes[0]]
    Edis = 0
    for i in range(len(data_array)):
        if list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2 : #or data_array[i][:2] == (atomtype2, atomtype1)) :
            Edis = data_array[i][2]
    if not Edis:
        #print("Info: Used simplified atom types to determine dissociaton energy for morse potential")
        logging.debug("Used simplified atom types to determine dissociaton energy for morse potential")
        comb1 = [atomtypes[0][0], atomtypes[1][0]]
        comb2 = [atomtypes[1][0], atomtypes[0][0]]
        for i in range(len(data_array)):
            if list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2 : #or data_array[i][:2] == (atomtype2, atomtype1)) :
                Edis = data_array[i][2]
    #print ('E_dis = ' + str(Edis))
    logging.debug('E_dis = ' + str(Edis))
    if not Edis:
        print("Warning: No morse dissociation energy found. Used 350 as default value")
        logging.warning("No morse dissociation energy found. Used 350 as default value")
        Edis = 350   #default value in same order of magnitude
    return Edis



def find_distances (plumedfile, datafile):
    data_all, data_array = get_data_from_file(plumedfile)
    
    #get all pairs from plumedfile
    list_of_pairs_and_distances = []
    nbr_of_pairs = 0
    for i in range(len(data_all)):
        if 'DISTANCE' in data_all[i]:
            if 'broken' in data_all[i]: #leave out already broken distances
                continue
            split1 = data_all[i].split(',')
            atom2 = split1[-1][:-2]
            split2 = split1[0].split('=')
            atom1 = split2[-1]
            nbr_of_pairs += 1 
            list_of_pairs_and_distances.append([str(atom1), str(atom2)])
   
    #get all distances from datafile
    #Note: Open file line by line i.o.t. to avoid memory error  
    list_of_distances = []
    header = True #used to skip header
    ctr = 0
    with open (datafile) as f:
        for line in f:
            line_split = []
            ctr += 1
            if header:
                header = False
                continue
            line_split = np.asarray(line.split())
            for k in range(1, nbr_of_pairs + 1): #shift by +1 i.o.t. leave out time (first entry)
                distance = float(line_split[k])
                list_of_pairs_and_distances[k-1].append(distance)
    nbr_of_data_points = ctr - 1

    print ('Collected distances from ' + str(datafile) + ' for ' + str(nbr_of_pairs) + ' pairs with ' + str(nbr_of_data_points) + ' distances per pair.')
    logging.info('Collected distances from ' + str(datafile) + ' for ' + str(nbr_of_pairs) + ' pairs with ' + str(nbr_of_data_points) + ' distances per pair.')
    
    return list_of_pairs_and_distances
              

def calc_transition_rate(r_curr, r_0, E_dis, k_f):
    #parameters
    kT = 2.479      #k_B T at 310K #in Gromacs units kj *mol^-1
    #tau =  0.16    #unfitted / theoretical pre-exponential factor #h/kT = 0.16 ps from transition state theory 
    k_0 =  0.288    #pre-exponential factor #from fitting averaged C_a - N data to gromacs data, see paper  #or: 1/2pi sqrt(k/m)
    
    #calculates energy barrier crossing rate [in ps]; barrier based on the model V = V_morse - F*X
    
    beta = np.sqrt(k_f / (2*E_dis))  #[beta] =1/nm since beta = sqrt(k/2D)

    #calculate inflection point corresponding to point with maximal force
    r_infl = (beta*r_0 + np.log(2)) / beta
        
    #calculate current force in bond F = -del V / del x
    if r_curr > r_infl:
        logging.debug('Used maximum force for bond Evans model since position behind inflection point found')
        F = 2*beta*E_dis*np.exp(-beta*(r_infl-r_0))*(1-np.exp(-beta*(r_infl-r_0)))
    else:
        F = 2*beta*E_dis*np.exp(-beta*(r_curr-r_0))*(1-np.exp(-beta*(r_curr-r_0)))
    logging.debug('Calculated force in bond F = ' + str(F))

    
    #calculate extrema of shifted potential i.o.t. get barrier hight
    rmin = r_0 - 1/beta * np.log((beta * E_dis + np.sqrt(beta**2 * E_dis **2 - 2*E_dis*beta*F))/(2*beta*E_dis))
    rmax = r_0 - 1/beta * np.log((beta * E_dis - np.sqrt(beta**2 * E_dis **2 - 2*E_dis*beta*F))/(2*beta*E_dis))

    Vmax = E_dis*(1-np.exp(-beta*(rmax-r_0)))**2 - F * (rmax - r_0)   #Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero
    Vmin = E_dis*(1-np.exp(-beta*(rmin-r_0)))**2 - F * (rmin - r_0)

    delta_V = Vmax - Vmin

    k = k_0 * np.exp(- delta_V/kT)     #[1/ps]

    if float(r_curr) > rmax:   #already jumped over the barrier? Even if not "open" in gromacs morse potential?
        pass
    if F <= 0.0:  #negative force: Vmax -> infinity impliying k -> 0
        k = 0.0
        logging.info('Found negative force, most likely due to compression. Rate replaced with zero.')

    return k, F   

    
def calc_av_rate(distances, r_0, E_dis, k_f):
    #average distances first, if necessary
    dist = []
    if len(distances) > 1:    
        r_av = sum(distances[2:]) / len(distances[2:])
    else:
        r_av = float(distances[0])
        print (r_av)
    k, F = calc_transition_rate(r_av, r_0, E_dis, k_f)
    print(r_av, k,F)

    return k


def delete_bonded_interactions(oldtop, newtop, breakpair):
    #function that cuts topology into the parts at the breakpair deleting all bonded interactions

    data_all, data_array = get_data_from_file(oldtop)

    print('Info: Start deleting bond ' + str(breakpair) + ' in Topology')
    logging.info('Info: Start deleting bond ' + str(breakpair) +' in Topology')

    proper_set = False
    list_of_pairs = []
    list_of_dihedrals = []
    deleted = 0
    reached_pairs = False
    data_new = data_all[:]   #get an independent copy with slice

    #remove bond, angles and dihedrals where breakpair was involved
    #save pairs, dihedrals and their respective positions in new data set for next step
    for i in range(len(data_all)):
        if breakpair[0] in data_array[i][0:4] and breakpair[1] in data_array[i][0:4]:    #[0:4] since value "func" should not be considered (leads e.g. to problem if '1' in breakpair)
            #print ("Deleted: " + data_new[i - deleted])  #shift -
            logging.debug ('Deleted in .top: '+  data_new[i - deleted])
            data_new.remove(data_new[i-deleted])
            deleted = deleted + 1
        if "[ bonds ]" in data_all[i]:
            bonds = i
        if "[ pairs ]" in data_all[i]:
            pairs = i   #Note: This is the position in the NEW data since removel of all above interaction has already happend
            reached_pairs = True
        if reached_pairs and i > pairs +1:
            list_of_pairs.append(list(data_array[i][0:2]))
        if "[ angles ]" in data_all[i]:
            angles = i
            reached_pairs = False
            list_of_pairs = list_of_pairs[:-2]
        if "[ dihedrals ]" in data_all[i]:
            if not proper_set:
                dihedrals = i
                proper_set = True
            else:
                improper = i
                proper_set = False
        if proper_set and i > dihedrals +1:
            list_of_dihedrals.append(data_array[i])
                

    list_of_dihedrals = list_of_dihedrals[:-1]
    deleted_pairs = 0

    #go through all dihedrals and thus find pairs to be deleted if breakpair is in the dihedral
    logging.debug('Deleted bonds, angles and dihedrals. Now starting deleting pairs: ')
    pairs_to_be_deleted = []
    for j in range(len(list_of_dihedrals)):
        if (breakpair[0] in list_of_dihedrals[j][0:4]) and (breakpair[1] in list_of_dihedrals[j][0:4]):
            if float(list_of_dihedrals[j][0]) < float(list_of_dihedrals[j][3]): #pairs are sorted 
                pair = [list_of_dihedrals[j][0], list_of_dihedrals[j][3]]
            else:
                pair = [list_of_dihedrals[j][3], list_of_dihedrals[j][0]]
            pairs_to_be_deleted.append(pair)
    for k in range(len(list_of_pairs)):
        if list_of_pairs[k] in pairs_to_be_deleted:
            #print ("Deleted: " + data_new[k + pairs - deleted_pairs + 1])  #shift -
            logging.debug('Deleted in .top:' + data_new[k + pairs - deleted_pairs + 1])
            data_new.remove(data_new[k + pairs - deleted_pairs +1 ])
            deleted_pairs += 1
                    
    store_linelist_to_file(data_new, newtop)


def modify_radical_interactions(oldtop, newtop, old_bond, reference_bond, old_to_new):    
    #function that changes topology: deletes old bond, creates new bonded interactions for reaction_pair by comparing to reference bond and changes the reference from old to new

    #delete old bond from new topology
    delete_bonded_interactions(oldtop, newtop, old_bond)
    data_all, data_array = get_data_from_file(newtop)
    reference_bond_str = str(reference_bond)
    newbond = reference_bond_str.replace(old_to_new[0],old_to_new[1],1)

    print('Info: Start creating new bond ' + newbond +' in Topology')
    logging.info('Info: Start reating new bond ' + newbond +' in Topology')

    proper_set = False
    list_of_pairs = []
    list_of_dihedrals = []
    reached_pairs = False
    data_new = data_all[:]   #get an independent copy with slice
    modified = 0

    #copy bond, angles and dihedrals of the old breakpair and change old binding partner to the new one
    #save pairs, dihedrals and their respective positions in new data set for next step
    for i in range(len(data_all)):
        if reference_bond[0] in data_array[i][0:4] and reference_bond[1] in data_array[i][0:4]:    #modify the ones with old breakpair #[0:4] since value "func" should not be considered (leads e.g. to problem if '1' in breakpair)
            logging.debug ('Copied in .top: ' + data_all[i])
            new_interaction_to_be_inserted = data_all[i]
            new_interaction_to_be_inserted = new_interaction_to_be_inserted.replace(old_to_new[0],old_to_new[1],1) #replace (first) occurance of old bond with new bonding partner 
            data_new = np.append(data_new, new_interaction_to_be_inserted)
            #old_interaction = '; former interaction: ' + str(data_all[i])
            #data_new = np.append(data_new, old_interaction)
            #modified = modified + 1
        else:
            data_new = np.append(data_new, data_all[i])
            
        if "[ bonds ]" in data_all[i]:
            bonds = i
        if "[ pairs ]" in data_all[i]:
            pairs = i   #Note: This is the position in the NEW data since removel of all above interaction has already happend
            reached_pairs = True
        if reached_pairs and i > pairs +1:
            list_of_pairs.append(list(data_array[i][0:2]))
        if "[ angles ]" in data_all[i]:
            angles = i
            reached_pairs = False
            list_of_pairs = list_of_pairs[:-2]
        if "[ dihedrals ]" in data_all[i]:
            if not proper_set:
                dihedrals = i
                proper_set = True
            else:
                improper = i
                proper_set = False
        if proper_set and i > dihedrals +1:
            list_of_dihedrals.append(data_array[i])                

    list_of_dihedrals = list_of_dihedrals[:-1]

    #go through all dihedrals and thus find pairs to be modified if reference bond is in the dihedral. Only take the ones where new bonding partner is directly involved and delete the others
    logging.debug('Added bonds, angles and dihedrals. Now start finding pairs: ')
    pairs_to_be_modified = []
    deleted_pairs = 0 
    for j in range(len(list_of_dihedrals)):
        if (reference_bond[0] in list_of_dihedrals[j][0:4]) and (reference_bond[1] in list_of_dihedrals[j][0:4]):
            if float(list_of_dihedrals[j][0]) < float(list_of_dihedrals[j][3]): #pairs are sorted 
                pair = [list_of_dihedrals[j][0], list_of_dihedrals[j][3]]
            else:
                pair = [list_of_dihedrals[j][3], list_of_dihedrals[j][0]]
            pairs_to_be_modified.append(pair)
    for k in range(len(list_of_pairs)):
        if list_of_pairs[k] in pairs_to_be_modified:
            if old_to_new[0] in list_of_pairs[k]:
                old_interaction = data_all[k + pairs +  2]
                new_interaction_to_be_inserted = old_interaction.replace(old_to_new[0],old_to_new[1],1) #replace (first) occurance of old bond with new bonding partner, if in there
                print ('will modify pair ' + old_interaction + 'to '  +  new_interaction_to_be_inserted )
                logging.debug('will modify pair ' + old_interaction + 'to '  +  new_interaction_to_be_inserted )
                data_new = np.insert(data_new, k + pairs - deleted_pairs + 2 ,  new_interaction_to_be_inserted)
                data_new = np.delete(data_new, k + pairs - deleted_pairs + 3) # delete old pair 
            else: #delete other pairs
                logging.debug('will delete pair ' + data_new[k + pairs - deleted_pairs + 2])
                print ('will delete pair ' + data_new[k + pairs - deleted_pairs + 2])
                data_new= np.delete(data_new, k + pairs - deleted_pairs +2)
                deleted_pairs += 1
                    
    store_linelist_to_file(data_new, newtop)

def modify_hbond(referencetop, oldtop, newtop, old_bond, reference_bond, old_to_new):    
    #function that changes topology: deletes old bond, creates new bonded interactions by comparing to reference bond and changes the reference from old to new

    #delete old bond from new topology
    delete_bonded_interactions(oldtop, newtop, old_bond)
    data_all, data_array = get_data_from_file(newtop)
    intermediatetop = 'broken_topol' + str(old_bond) + '.top'
    delete_bonded_interactions(referencetop, intermediatetop, old_bond)
    reference_all, reference_array = get_data_from_file(intermediatetop)  #original topology that includes the new h-bond but not the old one anymore
    reference_bond_str = str(reference_bond)
    newbond = reference_bond_str.replace(old_to_new[0],old_to_new[1],1)

    print('Info: Start creating new bond ' + newbond +' in Topology')
    logging.info('Info: Start reating new bond ' + newbond +' in Topology')

    proper_set = False

    #store positions in new topology so that we now where to insert the interactions
    for j in range(len(data_all)):
        if "[ bonds ]" in data_all[j]:
            bonds = j
        if "[ pairs ]" in data_all[j]:
            pairs = j
        if "[ angles ]" in data_all[j]:
            angles = j
        if "[ dihedrals ]" in data_all[j]:
            if not proper_set:
                dihedrals = j
                proper_set = True
            else:
                improper = j

    proper_set = False
    list_of_pairs = []
    list_of_dihedrals = []
    reached_pairs = False
    reached_dihedrals = False
    reached_angles = False
    data_new = data_all[:]   #get an independent copy with slice
    inserted = 0

    #copy bond, angles and dihedrals of the old breakpair and change old binding partner to the new one
    #save pairs, dihedrals and their respective positions in new data set for next step
    for i in range(len(reference_all)):
        if reference_bond[0] in reference_array[i][0:4] and reference_bond[1] in reference_array[i][0:4]:    #modify the ones with old breakpair #[0:4] since value "func" should not be considered (leads e.g. to problem if '1' in breakpair)
            logging.debug ('Copied in .top: ' + reference_all[i])
            new_interaction_to_be_inserted = reference_all[i]
            new_interaction_to_be_inserted = new_interaction_to_be_inserted.replace(old_to_new[0],old_to_new[1],1) #replace (first) occurance of old bond with new bonding partner 
            if reached_angles == False:
                #print (' add  bond' + str(new_interaction_to_be_inserted))
                logging.debug(' add  bond' + str(new_interaction_to_be_inserted))
                data_new = np.insert(data_new, bonds + 2 ,new_interaction_to_be_inserted)
                inserted = inserted +1
            elif reached_dihedrals == False:
                #print (' add  angle' + str(new_interaction_to_be_inserted))
                logging.debug(' add  angle' + str(new_interaction_to_be_inserted))
                data_new = np.insert(data_new, angles +2 + inserted ,new_interaction_to_be_inserted)
                inserted = inserted +1
            elif reached_dihedrals == True and proper_set == True: 
                data_new = np.insert(data_new, dihedrals +2 + inserted,new_interaction_to_be_inserted)
                print (' add  dihedral' + str(new_interaction_to_be_inserted))
                logging.debug(' add  dihedral' + str(new_interaction_to_be_inserted))
                inserted = inserted +1
            elif reached_dihedrals == True and proper_set == False:
                data_new = np.insert(data_new, improper +2 + inserted ,new_interaction_to_be_inserted)
                #print (' add  improper' + str(new_interaction_to_be_inserted))
                logging.debug (' add  improper' + str(new_interaction_to_be_inserted))
                inserted = inserted +1
            
        if "[ pairs ]" in reference_all[i]:
            pairs_ref = i
            reached_pairs = True
        if reached_pairs and i > pairs_ref +1:
            list_of_pairs.append(list(reference_all[i][0:2]))
        if "[ angles ]" in reference_all[i]:
            reached_angles = True
            reached_pairs = False
            list_of_pairs = list_of_pairs[:-2]
        if "[ dihedrals ]" in reference_all[i]:
            reached_dihedrals = True
            if not proper_set:
                dihedrals_ref = i
                proper_set = True
            else:
                improper_ref = i
                proper_set = False
        if proper_set and i > dihedrals +1:
            list_of_dihedrals.append(reference_array[i])                

    list_of_dihedrals = list_of_dihedrals[:-1]

    #go through all dihedrals and thus find pairs to be modified if reference bond is in the dihedral. 
    logging.debug('Added bonds, angles and dihedrals. Now start finding pairs: ')
    pairs_to_be_modified = []
    deleted_pairs = 0
    for j in range(len(list_of_dihedrals)):
        if (reference_bond[0] in list_of_dihedrals[j][0:4]) and (reference_bond[1] in list_of_dihedrals[j][0:4]):
            if float(list_of_dihedrals[j][0]) < float(list_of_dihedrals[j][3]): #pairs are sorted 
                pair = [list_of_dihedrals[j][0], list_of_dihedrals[j][3]]
            else:
                pair = [list_of_dihedrals[j][3], list_of_dihedrals[j][0]]
            pairs_to_be_modified.append(pair)
            
    inserted = 0
    for pair in pairs_to_be_modified:
        pair_str = str(pair[0]) + ' ' + str(pair[1])
        if old_to_new[0] in pair:
            new_interaction_to_be_inserted = pair_str.replace(old_to_new[0],old_to_new[1],1) #replace (first) occurance of old bond with new bonding partner, if in there
            new_interaction_to_be_inserted = new_interaction_to_be_inserted + '     1 \n'  #add function only after replacement to avoid mistakes e.g. by confusing atom 1 with function 1
            #print ('add pair '+  new_interaction_to_be_inserted + ' which was formerly ' + str(pair) )
            logging.debug('add pair '+  new_interaction_to_be_inserted + ' which was formerly ' + str(pair) )
            data_new = np.insert(data_new, pairs + inserted + 3,  new_interaction_to_be_inserted)
            inserted = inserted + 1
        else:
            pair_str = pair_str + '     1'
            logging.warning('Warning: will add pair that might have not been there before: ' + pair_str)
            print ('Warning: will add pair that might have not been there before: ' + pair_str)
            data_new= np.insert(data_new, pairs + inserted + 3, pair_str)
            inserted = inserted + 1            
                    
    store_linelist_to_file(data_new, newtop)
    
    modify_atomtype(newtop, referencetop, old_to_new)

    
def modify_atomtype (topology, reference_topology, old_to_new):
    # modifies atomtype in topology of the new atom to the atomtype of the old atom (that was at that position before) read out from the original reference topology

    reference_atom = old_to_new[0]
    atom = old_to_new[1]

    reference_all, reference_array = get_data_from_file(reference_topology)
    atoms = False
    for j in range(len(reference_all)):
        if 'atoms' in reference_all[j]:
            atoms = True
        if atoms and reference_atom in reference_array[j][0]:
            reference_type = reference_array[j][1]
            reference_charge = reference_array[j][6]
            print (' will change atomtype of ' + str(atom) + ' to type ' + str (reference_type) + ' with partial charge ' + str(reference_charge))
            logging.debug (' will change atomtype of ' + str(atom) + ' to type ' + str (reference_type) + ' with partial charge ' + str(reference_charge))
            break                       

    data_all, data_array = get_data_from_file(topology)
    atoms = False
    for i in range(len(data_all)):
        if 'atoms' in data_all[i]:
            atoms = True      
        if atoms and atom in data_array[i][0]:
            atom_type = data_array[i][1]
            atom_charge = data_array[i][6]
            line = data_all[i]
            print (line)
            line = line.replace(atom_type, reference_type)
            line = line.replace(atom_charge, reference_charge)
            print (line)
            data_all = np.delete(data_all, i)
            data_all = np.insert(data_all, i, line)            
            break

    store_linelist_to_file(data_all, topology)

    #maybe improve memory usage by reading in file line by line only, as indicated below
    '''
    #file = open(topology, 'wr')
    found_atom = False
    found_reference = False
    with open (topology) as f:
    for line in f:

        
        if '[ bonds ] ' in line:
            break
        
    file.close()
    '''

    
def modify_plumedfile(plumedfile, plumedfile_new, distancefile_new, ruptured_bonds):
    nbr_of_ruptures = len(ruptured_bonds)
    rupture_pairs = []
    cond_groups = []
    broken_distnbr = []
    for i in range (nbr_of_ruptures):
        rupture_pairs.append(ruptured_bonds[i][0])

    file = open(plumedfile_new, "wr")   #open in  append mode
    with open (plumedfile) as f:
        for line in f:
            if 'ATOMS=' in line:
                split1 = line.split(',')
                atom2 = split1[-1][:-2]
                split2 = split1[0].split('=')
                atom1 = split2[-1]
                split3 = split2[0].split(':')
                dist_nbr = split3[0]
                if [atom1, atom2] in rupture_pairs:
                    line = '# --already broken-- ' + line
                    broken_distnbr.append(dist_nbr)
            if 'PRINT' in line:
                line.find(dist_nbr)
                for dist_nbr in broken_distnbr:
                    line = line.replace(dist_nbr + ',', '')
                split = line.split()
                file_old = split[-1]
                line = line.replace(file_old, 'FILE=' + str(distancefile_new))
                
            file.write(line)
    file.close()
    print ('Adjusted Plumedfile: Commented out distances '+ str(broken_distnbr) + ' corresponding to breakpairs: ' + str(rupture_pairs) + '. \n Will now write distances to: ' + str(distancefile_new))
    logging.info('Adjusted Plumedfile: Commented out distances '+ str(broken_distnbr) + ' corresponding to breakpairs: ' + str(rupture_pairs) + '. \n Will now write distances to: ' + str(distancefile_new))
    

def find_candidates(grofile, trrfile, cutoff_distance, reaction_time, radical):
    topology = grofile
    trajectory = trrfile
    radical_atom = u.select_atoms('bynum ' + str(radical))
    
    u = MDAnalysis.Universe(topology, trajectory)
    print (u, ' with frames in trajectory = ' + str(len(u.trajectory)))
    logging.info (u, ' with frames in trajectory = ' + str(len(u.trajectory)))
    candidates = np.array([])
    time = np.array([])
    frame_ctr = 0
    distances_rad = np.array([])

    # at t = 0
    nearR = u.select_atoms('around ' + str(distance) + ' bynum '+ str(radical))
    H_nearR = nearR.select_atoms('type H')
    print ('H-Atoms nearby '+ str(radical) +' (max. distance = '+str(distance) + ' Angstroem) at point of breakage: ')
    for atom in H_nearR:
        print (atom)

    u_curr = u.select_atoms('around ' + str(distance) + ' bynum '+str(radical), updating = True)
    u_currH = u_curr.select_atoms('type H', updating = True)

    for ts in u.trajectory:
        if u.trajectory.time > reaction_time:
            break
        H_nearR += u_currH
        time = np.append(time, u.trajectory.time)
        candidates = np.append(candidates, len(u_curr))
        frame_ctr += 1
        radical_pos = np.split(radical_atom.position,3)
    
        
            
    H_nearR = H_nearR.select_atoms('type H') #deleting duplicates
    print ('H-Atoms nearby '+ str(radical) +' (at least once within max. distance = '+str(distance) + ' Angstroem) between break and '
           + str(reaction_time) +'ps after breakage (number of frames = ' + str(frame_ctr) + ') :')
    for atom in H_nearR:
        print (atom)

    print ('Sorted by resnames: ' + str(H_nearR.groupby(['resnames'])))
    print ('Sorted by atom names: ' + str(H_nearR.groupby(['names'])))
    print ('Now analyzing the distances of these candidates to the radical in the aforementioned time window after the breakage:')

    #average distance for candidates
    dist_dic = {}
    list_av_dist = np.array([])
    for atom in H_nearR:
        distances = np.array([])          
        for ts in u.trajectory:
            if u.trajectory.time > reaction_time:
                break
            atom_pos = np.split(atom.position,3)
            radical_pos = np.split(radical_atom.position,3)
            dist = np.sqrt ((radical_pos[0] - atom_pos[0])**2 + (radical_pos[1] - atom_pos[1])**2 +  (radical_pos[2] - atom_pos[2])**2)
            distances = np.append(distances, dist)
        av_distance = np.average(distances)
        if av_distance <= 1.4: #skip H-atoms that are already bound to the radical
            continue
        else:
            list_av_dist = np.append(list_av_dist, av_distance)
            dist_dic.update({av_distance : atom})
            #print (atom.index, 'has average distance to radical [Angstroem]: ' + str(av_distance))

    min_dist = np.amin(list_av_dist)
    min_candidate = dist_dic[min_dist]
    
    print ('Average distances and their respective atoms: ', dist_dic)
    print (' Candidate with smallest average distance ' + str(min_dist) + ' = ' + str (min_candidate))

    return min_candidate

    

def do_kinetic_mc(list_of_nbrs_and_atomtypes, list_of_rates):
    #rejection-free kinetic Monte Carlo 
    #compare e.g. https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo#Rejection-free_KMC

    total_rate = sum(list_of_rates)  #sum all rates
    random.seed()
    t = random.random() #t in [0.0,1.0)
    rate_running_sum = 0
    for i in range(len(list_of_rates)): 
        rate_running_sum += list_of_rates[i]
        if (t*total_rate) <= rate_running_sum:
            breakpair = list_of_nbrs_and_atomtypes[i][0]
            atomtypes = list_of_nbrs_and_atomtypes[i][1]
            rate = list_of_rates[i]                 
            break
    u = random.random()
    delta_t = np.log(1/u)/total_rate

    print ('Rate = ' + str(rate) + ', breakpair = ' + str (breakpair) + ', atomtypes = ' + str(atomtypes) + 'jump [ps] = ' + str(delta_t))
    logging.info(('Rate = ' + str(rate) + ', breakpair = ' + str (breakpair) + ', atomtypes = ' + str(atomtypes) + 'jump [ps] = ' + str(delta_t)))

    return breakpair, atomtypes, delta_t

def do_rejection_KMC (list_of_breakpairs, k_0, topfile, filepath_bonds, filepath_edis):
    #find candidate and calculate respective rupture rate
    nbr_of_breakpairs = len(list_of_breakpairs)
    random.seed()
    t = random.randint(0, nbr_of_breakpairs -1)  #choose index randomly. Note: Indices go from 0 to length-1 since computer people think differently
    breakpair = [list_of_breakpairs[t][0], list_of_breakpairs[t][1]]
    distances = list_of_breakpairs[t][2:]
    atomtypes = identify_atomtypes(breakpair,topfile)
    r_0, k_f = find_bond_param(atomtypes, filepath_bonds) #read out from gromacs force field
    E_dis = find_Edis(atomtypes, filepath_edis)  #read out from gromacs force field
    k = calc_av_rate(distances, float(r_0), float(E_dis), float(k_f)) 
    print ('Info: For ' + str(breakpair) + str(atomtypes) + ' as candidate in r-kMC calculated rupture rate using ' + str((len(distances))) + ' distances per bond: ' + str(k) + ' per ps.')
    logging.info('For ' + str(breakpair) + str(atomtypes) + ' as candidate in r-kMC calculated rupture rate using ' + str((len(distances))) + ' distances per bond: ' + str(k) + ' per ps.')

    u = random.random() # [0,1)
    acceptance_probability = float(k) / k_0
    if acceptance_probability > u:
        rupture = True
        v = random.random()
        delta_t = np.log(1/v)/(nbr_of_breakpairs*r_0)
    else:
        rupture = False
        delta_t = 0
    print ("---Info: RKMC Transition accepted? " + str(rupture) + "---")
    logging.info("--- RKMC Transition accepted? = " + str(rupture) + "---")

    return rupture, breakpair, atomtypes, delta_t


def check_if_error_occured(filepath):
    #Returns True if simulation ended or error occured
    error = False
    data_all, data_array = get_data_from_file(filepath)
    for i in range(len(data_array)):
        if 'error' in data_all[i]:
            error = True
            print ("Warning: Simulation seems to have stopped due to an error: " + str(data_all[i:i+4]))
            logging.warning("Simulation seems to have stopped due to an error:" + str(data_all[i:i+4]))
            
    return error


def del_backup_files_and_step_files ():
    #deletes backup_files created by gromacs, since more than 99 files can not be handeld at once

    import os

    files = os.listdir('.')
    for i in range(len(files)):
        if files[i][0] == '#':
            os.remove(files[i])
            print (files[i] + " deleted")
            logging.info(files[i] + " deleted")
        elif "".join(files[i][0:4]) == 'step':
            os.remove(files[i])
            print (files[i] + " deleted")
            logging.info(files[i] + " deleted")


def write_conditions_in_plumedfile(topfile, indexfile, indexgroup, parameterfile):
    #use once to create index and plumed file with all necessary entries / conditons    # to be adjusted system specific.
    #uses atoms from given indexgroup (and, hardcoded, crosslinks, if not commented out). Skips bonds that include hydrogens or oxygens since they are not break-relevant.
    
    data_all, data_array = get_data_from_file(topfile)

    #create dic_of_atoms with atomtype and number to sort out hydrogens afterwards
    dic_of_atoms = {}
    for i in range(len(data_all)):
        if "[ bonds ]" in data_all[i]:   #only go until atoms ended and continue here later
            bonds_pos = i        
            break
        
        if len(data_array[i]) == 0: #skip empty lines
            pass
        elif data_all[i][0] in (';', '#', '[', 'P'):  #skip comments, types, includes and Protein_Chains
            pass
        else:
            atom_nbr = data_array[i][0]
            atom_type = data_array[i][1]
            dic_of_atoms.update({atom_nbr : atom_type})
    dic_of_indx = {}
    list_of_non_h_bonds = []
    index_nbr = 1 #start at 1 since cond-stop has default group 0 alreay

    #collect backbone atoms from indexfile
    index_all, index_array = get_data_from_file(indexfile)
    relevant_atoms = []
    found_group = False
    for k in range(len(index_array)): 
        if indexgroup in index_all[k]:
            found_group = True
        if found_group:
            for element in index_array[k]:
                relevant_atoms.append(element)

        if found_group and '[ ' in index_all[k+1]:  #stop at next group
            break
        

    for j in range(bonds_pos+2, len(data_all)): #go through all bonds

        if "[ pairs ]" in data_all[j]:   #stop when done with all bonds
            break
        
        if len(data_array[j]) > 0:
            nbr1 = data_array[j][0]
            nbr2 = data_array[j][1]
            #skip irrelevant bonds (e.g. side chains which are not under force)
            if nbr1 not in relevant_atoms:
                continue
            if nbr2 not in relevant_atoms:
                continue
        #leave out stronger C-N bond (chemically not prone to rupture due to electron resonance) 
        if ('C', 'N') in [(dic_of_atoms[nbr1], dic_of_atoms[nbr2]), (dic_of_atoms[nbr2], dic_of_atoms[nbr1])]:  
            pass
        elif "H" in dic_of_atoms[nbr1] or "H" in dic_of_atoms[nbr2]: #leave out hydrogen bonds
            pass
        elif "O" in dic_of_atoms[nbr1] or "O" in dic_of_atoms[nbr2]: #leave out oxygen bonds
            pass        
        else:
            bond = (nbr1, nbr2)
            list_of_non_h_bonds.append(bond)
            if nbr1 not in dic_of_indx.keys():
                index_name = str(index_nbr) #+ '_' + dic_of_atoms[nbr1]
                dic_of_indx.update( {nbr1 : index_name} )
                index_nbr += 1
            if nbr2 not in dic_of_indx.keys():
                index_name = str(index_nbr) #+ '_' + dic_of_atoms[nbr2]
                dic_of_indx.update( {nbr2 : index_name} )
                index_nbr += 1

    cond_ngroups = len(dic_of_indx.keys())
    cond_nconds = len(list_of_non_h_bonds)

    #write plumed-file
    print_arg = ''
    file = open('plumed.dat', "wr")   #open in  append mode
    file.write ('#Define distances \n')
    for pair_nbr in range(cond_nconds):
        nbr1 = list_of_non_h_bonds[pair_nbr][0]
        nbr2 = list_of_non_h_bonds[pair_nbr][1]
        file.write('d' + str(pair_nbr) + ': DISTANCE ATOMS=' + str(nbr1) + ',' +str(nbr2) + ' \n')
        print_arg += 'd' + str(pair_nbr) + ','
    file.write (' \n#Print distances ARG to FILE every STRIDE steps \n')
    file.write ('PRINT ARG=' + str(print_arg) + ' STRIDE=100 FILE=distances.dat')    
    file.close()
    
    print ("finished writing plumed-file")



if __name__ == "__main__":
                   
    #write_conditions_in_plumedfile('topol.top', 'index.ndx', '[ Backbone ]', 'ffbonded.itp')
    #modify_hbond('topol.top', 'reaction_topol_test_v5.top', 'reaction2_topol_test_v3.top',['17','19'], ['15','16'], ['16','19'])
    #modify_atomtype('reaction2_topol_test_v3_atomtype.top', 'topol.top', '19', '16')
    #delete_bonded_interactions('reaction_topol_test_v5.top', 'reaction_topol_test_v5.top', ['43','16'])
    pass
    
