import numpy as np


def get_data_from_file(filepath):
    file = open(filepath, "r")
    data_all = []  # array with each entry corresponding to one (string) line
    data_array = []  # array of all lines with each subarray containing one value
    settings = []
    for line in file:
        data_all.append(line)
        line_array = np.asarray(line.split())
        data_array.append(line_array)
    file.close()

    return data_all, data_array


def write_conditions_in_plumedfile(topfile, indexfile, indexgroup, outplumed):
    # use once to create index and plumed file with all necessary entries / conditons.
    # to be adjusted system specific.
    # uses atoms from given indexgroup (and, hardcoded, crosslinks, if not commented out). Skips bonds that include hydrogens or oxygens since they are not break-relevant.

    data_all, data_array = get_data_from_file(topfile)

    # create dic_of_atoms with atomtype and number to sort out hydrogens afterwards
    dic_of_atoms = {}
    for i in range(len(data_all)):
        if (
            "[ bonds ]" in data_all[i]
        ):  # only go until atoms ended and continue here later
            bonds_pos = i
            break

        if len(data_array[i]) == 0:  # skip empty lines
            pass
        elif data_all[i][0] in (
            ";",
            "#",
            "[",
            "P",
        ):  # skip comments, types, includes and Protein_Chains
            pass
        else:
            atom_nbr = data_array[i][0]
            atom_type = data_array[i][1]
            dic_of_atoms.update({atom_nbr: atom_type})
    dic_of_indx = {}
    list_of_non_h_bonds = []
    index_nbr = 1  # start at 1 since cond-stop has default group 0 alreay

    # collect backbone atoms from indexfile
    index_all, index_array = get_data_from_file(indexfile)
    relevant_atoms = []
    found_group = False
    for k in range(len(index_array)):
        if indexgroup in index_all[k]:
            found_group = True
        if found_group:
            for element in index_array[k]:
                relevant_atoms.append(element)

        if found_group and "[ " in index_all[k + 1]:  # stop at next group
            break

    for j in range(bonds_pos + 2, len(data_all)):  # go through all bonds
        if "[ pairs ]" in data_all[j]:  # stop when done with all bonds
            break

        if len(data_array[j]) > 0:
            nbr1 = data_array[j][0]
            nbr2 = data_array[j][1]
            # skip irrelevant bonds (e.g. side chains which are not under force)
            if nbr1 not in relevant_atoms:
                continue
            if nbr2 not in relevant_atoms:
                continue
        # leave out stronger C-N bond (chemically not prone to rupture due to electron resonance)
        if ("C", "N") in [
            (dic_of_atoms[nbr1], dic_of_atoms[nbr2]),
            (dic_of_atoms[nbr2], dic_of_atoms[nbr1]),
        ]:
            pass
        elif (
            "H" in dic_of_atoms[nbr1] or "H" in dic_of_atoms[nbr2]
        ):  # leave out hydrogen bonds
            pass
        elif (
            "O" in dic_of_atoms[nbr1] or "O" in dic_of_atoms[nbr2]
        ):  # leave out oxygen bonds
            pass
        else:
            bond = (nbr1, nbr2)
            list_of_non_h_bonds.append(bond)
            if nbr1 not in dic_of_indx.keys():
                index_name = str(index_nbr)  # + '_' + dic_of_atoms[nbr1]
                dic_of_indx.update({nbr1: index_name})
                index_nbr += 1
            if nbr2 not in dic_of_indx.keys():
                index_name = str(index_nbr)  # + '_' + dic_of_atoms[nbr2]
                dic_of_indx.update({nbr2: index_name})
                index_nbr += 1

    cond_ngroups = len(dic_of_indx.keys())
    cond_nconds = len(list_of_non_h_bonds)

    # write plumed-file
    print_arg = ""
    file = open(outplumed, "a")  # open in  append mode
    file.write("#Define distances \n")
    for pair_nbr in range(cond_nconds):
        nbr1 = list_of_non_h_bonds[pair_nbr][0]
        nbr2 = list_of_non_h_bonds[pair_nbr][1]
        file.write(
            "d"
            + str(pair_nbr)
            + ": DISTANCE ATOMS="
            + str(nbr1)
            + ","
            + str(nbr2)
            + " \n"
        )
        print_arg += "d" + str(pair_nbr) + ","
    file.write(" \n#Print distances ARG to FILE every STRIDE steps \n")
    file.write("PRINT ARG=" + str(print_arg) + " STRIDE=100 FILE=distances.dat")
    file.close()

    print("finished writing plumed-file")


topfile = "example/example_hexala/hexala_out.top"
indexfile = "example/example_hexala/index.ndx"
indexgroup = "Backbone"
outplumed = "example/example_hexala/plumed.dat"
write_conditions_in_plumedfile(topfile, indexfile, indexgroup, outplumed)
