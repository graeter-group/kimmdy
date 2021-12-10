import logging
from kimmdy.reaction import ConversionRecipe
from kimmdy.parsing import read_topol, write_topol


class ChangeManager:
    def __init__(self, recipe : ConversionRecipe):
        self.breakpair = recipe.atom_idx

    def modify_top(self, oldtop, newtop):
        logging.info(f"Reading: {oldtop} and writing modified topology to {newtop}.")
        topology = read_topol(oldtop)

        # remove bond, angles and dihedrals where breakpair was involved
        topology["bonds"] = [
            bond for bond in topology["bonds"]
            if not bond[0] in self.breakpair and not bond[1] in self.breakpair
            ]
        topology["pairs"] = [
            pair for pair in topology["pairs"]
            if not pair[0] in self.breakpair and not pair[1] in self.breakpair
            ]
        topology["angles"] = [
            angle for angle in topology["angles"]
            if not angle[0] in self.breakpair and not angle[1] in self.breakpair
            ]
        topology["dihedrals"] = [
            dihedral for dihedral in topology["dihedrals"]
            if not dihedral[1] in self.breakpair and not dihedral[2] in self.breakpair
            ]

        write_topol(topology, newtop)


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

