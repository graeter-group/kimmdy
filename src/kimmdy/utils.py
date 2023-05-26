import subprocess as sp
import numpy as np
import logging
from pathlib import Path
import MDAnalysis as MDA
from scipy.spatial.transform import Rotation

from kimmdy.config import Config


def increment_logfile(f: Path) -> Path:
    backup_file_prefix = "#"
    backup_file_suffix = "#"
    logfile = f
    if logfile.exists():
        backup_count = 1
        backup_file = (
            f"{backup_file_prefix}{logfile}_{backup_count}{backup_file_suffix}"
        )
        while Path(backup_file).exists():
            backup_count += 1
            backup_file = (
                f"{backup_file_prefix}{logfile}_{backup_count}{backup_file_suffix}"
            )
        logfile.rename(backup_file)
    return logfile


def find_bond_param(atomtypes, filepath):
    # reads bond parameters vom gromacs forcefield file #based on ff.bonded.itp from amber99sb-ildn.ff
    _, data_array = get_data_from_file(filepath)
    comb1 = [atomtypes[0], atomtypes[1]]
    comb2 = [atomtypes[1], atomtypes[0]]

    found_first_entry = (
        False  # only use first entry correspondnig to [ bondtypes ] section
    )
    for i in range(len(data_array)):
        if (
            list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2
        ):  # or data_array[i][:2] == (atomtype2, atomtype1)) :
            r_0 = data_array[i][3]
            k_f = data_array[i][4]
            found_first_entry = True
            break  # stop when first entry found

    if not found_first_entry:
        # print("Warning: No bond parameters found")
        # logging.warning('No bond parameters found')
        # r_0 = 0.145   #default value in same order of magnitude
        # k_f = 280000   #default value in same order of magnitude
        print("Found no entry for " + str(atomtypes))
        logging.error("Found no entry for " + str(atomtypes))

    return r_0, k_f


def find_Edis(atomtypes, filepath):
    # reads dissociation energy from edissoc.dat (gromacs)

    _, data_array = get_data_from_file(filepath)
    comb1 = [atomtypes[0], atomtypes[1]]
    comb2 = [atomtypes[1], atomtypes[0]]
    Edis = 0
    for i in range(len(data_array)):
        if (
            list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2
        ):  # or data_array[i][:2] == (atomtype2, atomtype1)) :
            Edis = data_array[i][2]
    if not Edis:
        # print("Info: Used simplified atom types to determine dissociaton energy for morse potential")
        # logging.debug("Used simplified atom types to determine dissociaton energy for morse potential")
        comb1 = [atomtypes[0][0], atomtypes[1][0]]
        comb2 = [atomtypes[1][0], atomtypes[0][0]]
        for i in range(len(data_array)):
            if (
                list(data_array[i][:2]) == comb1 or list(data_array[i][:2]) == comb2
            ):  # or data_array[i][:2] == (atomtype2, atomtype1)) :
                Edis = data_array[i][2]
    # print ('E_dis = ' + str(Edis))
    # logging.debug('E_dis = ' + str(Edis))
    if not Edis:
        # print("Warning: No morse dissociation energy found. Used 350 as default value")
        # logging.warning("No morse dissociation energy found. Used 350 as default value")
        Edis = 350  # default value in same order of magnitude
    return Edis


def find_distances(plumedfile, datafile):
    data_all, data_array = get_data_from_file(plumedfile)

    # get all pairs from plumedfile
    list_of_pairs_and_distances = []
    nbr_of_pairs = 0
    for i in range(len(data_all)):
        if "DISTANCE" in data_all[i]:
            if "broken" in data_all[i]:  # leave out already broken distances
                continue
            split1 = data_all[i].split(",")
            atom2 = split1[-1][:-2]
            split2 = split1[0].split("=")
            atom1 = split2[-1]
            nbr_of_pairs += 1
            list_of_pairs_and_distances.append([str(atom1), str(atom2)])

    # get all distances from datafile
    # Note: Open file line by line i.o.t. to avoid memory error
    list_of_distances = []
    header = True  # used to skip header
    ctr = 0
    with open(datafile) as f:
        for line in f:
            line_split = []
            ctr += 1
            if header:
                header = False
                continue
            line_split = np.asarray(line.split())
            for k in range(
                1, nbr_of_pairs + 1
            ):  # shift by +1 i.o.t. leave out time (first entry)
                distance = float(line_split[k])
                list_of_pairs_and_distances[k - 1].append(distance)
    nbr_of_data_points = ctr - 1

    # print ('Collected distances from ' + str(datafile) + ' for ' + str(nbr_of_pairs) + ' pairs with ' + str(nbr_of_data_points) + ' distances per pair.')
    # logging.info('Collected distances from ' + str(datafile) + ' for ' + str(nbr_of_pairs) + ' pairs with ' + str(nbr_of_data_points) + ' distances per pair.')

    return list_of_pairs_and_distances


def morse_transition_rate(
    r_curr: list[float],
    r_0: float,
    E_dis: float,
    k_f: float,
    k_0: float = 0.288,
    kT: float = 2.479,
):
    """calculates energy barrier crossing rate [in ps]; barrier based on the model V = V_morse - F*X"""
    r_curr = np.asarray(r_curr)
    beta = np.sqrt(k_f / (2 * E_dis))

    Fs = (
        2
        * beta
        * E_dis
        * np.exp(-beta * (r_curr - r_0))
        * (1 - np.exp(-beta * (r_curr - r_0)))
    )

    # inflection calculation
    r_infl = (beta * r_0 + np.log(2)) / beta
    F_infl = (
        2
        * beta
        * E_dis
        * np.exp(-beta * (r_infl - r_0))
        * (1 - np.exp(-beta * (r_infl - r_0)))
    )
    Fs_mask = r_curr > r_infl
    Fs[Fs_mask] = F_infl

    # calculate extrema of shifted potential i.o.t. get barrier hight
    rmin = r_0 - 1 / beta * np.log(
        (beta * E_dis + np.sqrt(beta**2 * E_dis**2 - 2 * E_dis * beta * Fs))
        / (2 * beta * E_dis)
    )
    rmax = r_0 - 1 / beta * np.log(
        (beta * E_dis - np.sqrt(beta**2 * E_dis**2 - 2 * E_dis * beta * Fs))
        / (2 * beta * E_dis)
    )

    Vmax = E_dis * (1 - np.exp(-beta * (rmax - r_0))) ** 2 - Fs * (rmax - r_0)
    Vmin = E_dis * (1 - np.exp(-beta * (rmin - r_0))) ** 2 - Fs * (rmin - r_0)
    # Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero

    delta_V = Vmax - Vmin
    k = k_0 * np.exp(-delta_V / kT)  # [1/ps]
    k_mask = Fs < 0
    k[k_mask] = 0

    return k


def calc_transition_rate(r_curr, r_0, E_dis, k_f):
    # parameters
    kT = 2.479  # k_B T at 310K #in Gromacs units kj *mol^-1
    # tau =  0.16    #unfitted / theoretical pre-exponential factor #h/kT = 0.16 ps from transition state theory
    k_0 = 0.288  # pre-exponential factor #from fitting averaged C_a - N data to gromacs data, see paper  #or: 1/2pi sqrt(k/m)

    #

    beta = np.sqrt(k_f / (2 * E_dis))  # [beta] =1/nm since beta = sqrt(k/2D)

    # calculate inflection point corresponding to point with maximal force
    r_infl = (beta * r_0 + np.log(2)) / beta

    # calculate current force in bond F = -del V / del x
    if r_curr > r_infl:
        # logging.debug('Used maximum force for bond Evans model since position behind inflection point found')
        F = (
            2
            * beta
            * E_dis
            * np.exp(-beta * (r_infl - r_0))
            * (1 - np.exp(-beta * (r_infl - r_0)))
        )
    else:
        F = (
            2
            * beta
            * E_dis
            * np.exp(-beta * (r_curr - r_0))
            * (1 - np.exp(-beta * (r_curr - r_0)))
        )
    # logging.debug('Calculated force in bond F = ' + str(F))

    # calculate extrema of shifted potential i.o.t. get barrier hight
    rmin = r_0 - 1 / beta * np.log(
        (beta * E_dis + np.sqrt(beta**2 * E_dis**2 - 2 * E_dis * beta * F))
        / (2 * beta * E_dis)
    )
    rmax = r_0 - 1 / beta * np.log(
        (beta * E_dis - np.sqrt(beta**2 * E_dis**2 - 2 * E_dis * beta * F))
        / (2 * beta * E_dis)
    )

    Vmax = E_dis * (1 - np.exp(-beta * (rmax - r_0))) ** 2 - F * (
        rmax - r_0
    )  # Note: F*r should lead to same result as F*(r-r_0) since the shifts in Vmax-Vmin adds up to zero
    Vmin = E_dis * (1 - np.exp(-beta * (rmin - r_0))) ** 2 - F * (rmin - r_0)

    delta_V = Vmax - Vmin

    k = k_0 * np.exp(-delta_V / kT)  # [1/ps]

    if (
        float(r_curr) > rmax
    ):  # already jumped over the barrier? Even if not "open" in gromacs morse potential?
        pass
    if F <= 0.0:  # negative force: Vmax -> infinity impliying k -> 0
        k = 0.0
        # logging.info('Found negative force, most likely due to compression. Rate replaced with zero.')

    return k, F  # [0,1]


def calc_av_rate(distances, r_0, E_dis, k_f):
    # average distances first, if necessary
    if len(distances) > 1:
        r_av = sum(distances[2:]) / len(distances[2:])
    else:
        r_av = float(distances[0])
        print(r_av)
    k, F = calc_transition_rate(r_av, r_0, E_dis, k_f)
    print(r_av, k, F)

    return k


# utils
def run_shell_cmd(s, cwd=None) -> sp.CompletedProcess:
    return sp.run(s, shell=True, cwd=cwd)


def get_shell_stdout(s):
    process = sp.run(s, shell=True, capture_output=True, encoding="utf-8")
    return process.stdout


def check_gmx_version(config: Config):
    """Check for an existing gromacs installation.

    If PLUMED is meant to be used it additionally checks for the keyword
    'MODIFIED' in the version name.
    """
    try:
        version = [
            l
            for l in get_shell_stdout(
                f"{config.gromacs_alias} --quiet --version"
            ).split("\n")
            if "GROMACS version:" in l
        ][0]
    except Exception as e:
        m = "No system gromacs detected. With error: " + str(e)
        logging.error(m)
        raise SystemError(m)

    # i hate this
    if (
        any(
            "plumed" in y
            for y in [
                config.mds.attr(x).get_attributes() for x in config.mds.get_attributes()
            ]
        )
        and not "MODIFIED" in version
    ):
        m = "GROMACS version does not contain MODIFIED, aborting due to lack of PLUMED patch."
        logging.error(m)
        logging.error("Version was: " + version)
        if not config.dryrun:
            raise SystemError(m)
    return version


def float_or_str(elem):
    try:
        return float(elem)
    except:
        return elem


## helpers for changemanager
def str_to_int_or_0(elem):
    return int(elem) if elem.isdigit() else 0


def sort_bond(entry):
    return sorted((str_to_int_or_0(entry[0]), str_to_int_or_0(entry[1])))


def sort_angle(entry):
    return (
        str_to_int_or_0(entry[1]),
        str_to_int_or_0(entry[0]),
        str_to_int_or_0(entry[2]),
    )


def sort_dihedral(entry):
    return (
        str_to_int_or_0(entry[1]),
        str_to_int_or_0(entry[2]),
        str_to_int_or_0(entry[0]),
        str_to_int_or_0(entry[3]),
    )


def sort_improper(entry):
    return (
        str_to_int_or_0(entry[2]),
        str_to_int_or_0(entry[0]),
        str_to_int_or_0(entry[1]),
        str_to_int_or_0(entry[3]),
    )


def check_idx(object):
    try:
        return str_to_int_or_0(object.idx)
    except:
        raise ValueError("Non Atom object in AtomList")


## from kimmdy
def store_linelist_to_file(data, filepath):
    file = open(filepath, "w")
    for line in data:
        file.write(line)
    file.close()


def identify_atomtypes(filepath):
    dic_of_atoms_to_groups = {}
    file = open(filepath, "r")
    atoms = False

    for line in file:
        # start collecting data when [ atoms ] reached
        if "atoms" in line:
            atoms = True
            continue
        if atoms == False:
            continue

        # reached [ bonds ], stop
        if "bonds" in line:
            break

        # leave out comments, includes, and empty lines
        if len(line) <= 1:
            continue
        elif line[0] == ";":
            continue
        elif line[0] == "#":
            continue

        # hopefully, these are all and exclusively the lines now with atomnbrs, types etc.
        line_array = np.asarray(line.split())
        nbr = line_array[0]
        atomtype = line_array[1]
        dic_of_atoms_to_groups.update({nbr: atomtype})
    file.close()

    return dic_of_atoms_to_groups


# TODO: replace with our new parsers
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


# TODO: replace with our new parsers
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


## coordinate changer utils


def find_radical_pos(
    center: MDA.core.groups.Atom, bonded: MDA.core.groups.AtomGroup, tetrahedral=False
):
    """Calculates possible radical positions of a given radical atom

    Parameters
    ----------
    center : MDA.core.groups.Atom
        Radical atom
    bonded : MDA.core.groups.AtomGroup
        Atom group of bonded atoms. From its length the geometry is predicted.
    tetrahedral : bool
        Whether to assume a tetrahedral conformation around C and N

    Returns
    -------
    list
        List of radical positions, three dimensional arrays
    """
    scale_C = 1.10
    scale_N = 1.04
    scale_O = 0.97
    scale_S = 1.41

    if center.element == "C" and len(bonded) == 2:
        c_alphas = center.bonded_atoms
        c_o = c_alphas[np.nonzero(c_alphas.elements == "O")]
        assert len(c_o) in [0, 1], "Carboxyl radical?"
        if len(c_o) > 0 and len(c_o[0].bonded_atoms) > 1:
            tetrahedral = True

    if tetrahedral:
        # MAKE SP3 from C with 2 bonds
        # invert bonded positions
        inv_1 = (center.position - bonded[0].position) + center.position
        inv_2 = (center.position - bonded[1].position) + center.position
        # construct rotation axis
        midpoint = (inv_2 - inv_1) * 0.5 + inv_1
        rot_ax = midpoint - center.position
        # 90 degree rotation
        r = Rotation.from_rotvec((np.pi / 2) * (rot_ax / np.linalg.norm(rot_ax)))
        # rotated bonds relative to center
        rad_1 = r.apply(inv_1 - center.position)
        rad_2 = r.apply(inv_2 - center.position)

        # scale to correct bond length, make absolute
        if center.element == "N":
            scale = scale_N
        elif center.element == "C":
            scale = scale_C
        else:
            raise NotImplementedError("H Bondlength to central atom missing")

        rad_1 = (rad_1 / np.linalg.norm(rad_1)) * scale + center.position
        rad_2 = (rad_2 / np.linalg.norm(rad_2)) * scale + center.position
        return [rad_1, rad_2]

    if len(bonded) in [2, 3]:
        assert center.element in [
            "C",
            "N",
        ], "Element element does not match bond number"

        # prediction: inverse midpoint between bonded
        # scale to correct bond length, make absolute
        if center.element == "N":
            scale = scale_N
        elif center.element == "C":
            scale = scale_C

        b_normed = []
        for b in bonded:
            b_vec = b.position - center.position
            b_vec_norm = b_vec / np.linalg.norm(b_vec)
            b_normed.append(b_vec_norm)

        midpoint = sum(b_normed)

        v = midpoint / np.linalg.norm(midpoint)
        rad = center.position + (-1 * v * scale)
        return [rad]

    # Radicals w/ only one bond:
    elif len(bonded) == 1:
        assert center.element in [
            "O",
            "S",
        ], "Element element does not match bond number"
        if center.element == "O":
            scale = scale_O
        elif center.element == "S":
            scale = scale_S

        b = bonded[0]
        b_vec = b.position - center.position
        b_vec = b_vec / np.linalg.norm(b_vec)
        rnd_vec = [1, 1, 1]  # to find a vector perpendicular to b_vec

        rnd_rot_ax = np.cross(b_vec, rnd_vec)
        rnd_rot_ax = rnd_rot_ax / np.linalg.norm(rnd_rot_ax)

        r1 = Rotation.from_rotvec(1.911 * rnd_rot_ax)  # 109.5 degree
        r2 = Rotation.from_rotvec(0.785 * b_vec)  # 45 degree

        ends = [r1.apply(b_vec)]  # up to 109.5

        for i in range(8):
            ends.append(r2.apply(ends[-1]))  # turn in 45d steps

        # norm and vec --> position
        ends = [(e / np.linalg.norm(e)) * scale + center.position for e in ends]

        return ends

    else:
        raise ValueError(f"Weired count of bonds: {list(bonded)}\nCorrect radicals?")
