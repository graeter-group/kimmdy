import subprocess as sp
import numpy as np


def run_shell_cmd(s, cwd=None):
    sp.run(s, shell=True, cwd=cwd)


def get_shell_stdout(s):
    process = sp.run(s, shell=True, capture_output=True, encoding="utf-8")
    return process.stdout


def check_gmx_version():
    return get_shell_stdout("gmx --quiet --version")


## from kimmdy
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
