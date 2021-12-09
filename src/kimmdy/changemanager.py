import logging
from kimmdy.utils import store_linelist_to_file, get_data_from_file


class ChangeManager:
    def __init__(self, breakpair):
        self.breakpair = breakpair

    def modify_top(self, oldtop, newtop):
        data_all, data_array = get_data_from_file(oldtop)

        # print('Info: Start modification of Topology')
        # logging.info('Start modification of Topology')

        proper_set = False
        list_of_pairs = []  # pairs
        list_of_dihedrals = []
        # number_of_atm = 0
        # possible_lower_pairpartner = []
        # possible_higher_pairpartner = []
        deleted = 0
        reached_pairs = False
        data_new = data_all[:]  # get an independent copy with slice

        # remove bond, angles and dihedrals where breakpair was involved
        # save pairs, dihedrals and their respective positions in new data set for next step
        for i in range(len(data_all)):
            if (
                self.breakpair[0] in data_array[i][0:4]
                and self.breakpair[1] in data_array[i][0:4]
            ):  # [0:4] since value "func" should not be considered (leads e.g. to problem if '1' in breakpair)
                # print ("Deleted: " + data_new[i - deleted])  #shift -
                # logging.debug ('Deleted in .top: '+  data_new[i - deleted])
                data_new.remove(data_new[i - deleted])
                deleted = deleted + 1
            if "[ bonds ]" in data_all[i]:
                bonds = i
                number_of_atms = data_array[bonds - 2][0]
            if "[ pairs ]" in data_all[i]:
                pairs = i  # Note: This is the position in the NEW data since removel of all above interaction has already happend
                reached_pairs = True
            if reached_pairs and i > pairs + 1:
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
            if proper_set and i > dihedrals + 1:
                list_of_dihedrals.append(data_array[i])

        list_of_dihedrals = list_of_dihedrals[:-1]
        deleted_pairs = 0

        # go through all dihedrals and thus find pairs to be deleted if breakpair is in middle of dihedral
        # logging.debug('Deleted bonds, angles and dihedrals. Now starting deleting pairs: ')
        pairs_to_be_deleted = []
        for j in range(len(list_of_dihedrals)):
            if (self.breakpair[0] in list_of_dihedrals[j][0:4]) and (
                self.breakpair[1] in list_of_dihedrals[j][0:4]
            ):
                if float(list_of_dihedrals[j][0]) < float(
                    list_of_dihedrals[j][3]
                ):  # pairs are sorted
                    pair = [list_of_dihedrals[j][0], list_of_dihedrals[j][3]]
                else:
                    pair = [list_of_dihedrals[j][3], list_of_dihedrals[j][0]]
                pairs_to_be_deleted.append(pair)
        for k in range(len(list_of_pairs)):
            if list_of_pairs[k] in pairs_to_be_deleted:
                # print ("Deleted: " + data_new[k + pairs - deleted_pairs + 1])  #shift -
                # logging.debug('Deleted in .top:' + data_new[k + pairs - deleted_pairs + 1])
                data_new.remove(data_new[k + pairs - deleted_pairs + 1])
                deleted_pairs += 1

        store_linelist_to_file(data_new, newtop)
        return newtop

    def modify_plumed(self, oldplumeddat, newplumeddat, newplumeddist):  # essential
        broken_distnbr = []

        file = open(newplumeddat, "w")  # open in  append mode
        with open(oldplumeddat) as f:
            for line in f:
                if "ATOMS=" in line:
                    split1 = line.split(",")
                    atom2 = split1[-1][:-2]
                    split2 = split1[0].split("=")
                    atom1 = split2[-1]
                    split3 = split2[0].split(":")
                    dist_nbr = split3[0]
                    if [atom1, atom2] in self.breakpair:
                        line = "# --already broken-- " + line
                        broken_distnbr.append(dist_nbr)
                if "PRINT" in line:
                    line.find(dist_nbr)
                    for dist_nbr in broken_distnbr:
                        line = line.replace(dist_nbr + ",", "")
                    split = line.split()
                    file_old = split[-1]
                    line = line.replace(file_old, "FILE=" + str(newplumeddist))

                file.write(line)
        file.close()
        return newplumeddat, newplumeddist
        # print ('Adjusted Plumedfile: Commented out distances '+ str(broken_distnbr) + ' corresponding to breakpairs: ' + str(rupture_pairs) + '. \n Will now write distances to: ' + str(distancefile_new))
        # logging.info('Adjusted Plumedfile: Commented out distances '+ str(broken_distnbr) + ' corresponding to breakpairs: ' + str(rupture_pairs) + '. \n Will now write distances to: ' + str(distancefile_new))
