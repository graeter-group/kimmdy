import argparse
from pathlib import Path

from kimmdy.topology.topology import Topology
from kimmdy.parsing import read_top


def create_plumed_dat(
    top_path: Path, idx_path: Path, idx_group: str, entries_to_remove: list
):
    print(
        "Creating plumed file with following settings:\n\n"
        f"top_path: \t\t{top_path}\n"
        f"idx_path: \t\t{idx_path}\n"
        f"idx_group: \t\t{idx_group}\n"
        f"entries_to_remove: \t\t{entries_to_remove}\n"
    )

    # init variables
    top = Topology(read_top(top_path))
    idx = read_top(idx_path)
    plumed_stride = 100
    plumed_distances_out = "distances.dat"

    # prepare remove variables
    atoms_to_remove = []
    atoms_start_to_remove = []
    bonds_to_remove = []
    for entry in entries_to_remove:
        entrysplit = entry.split(sep="-")
        if len(entrysplit) == 1:
            if entrysplit[0].endswith("*"):
                atoms_start_to_remove.append(entrysplit[0][:-1])
            else:
                atoms_to_remove.append(entrysplit[0])
        elif len(entrysplit) == 2:
            bonds_to_remove.append(set(entrysplit))
        else:
            print(f"Don't know how to deal with entry {entry}, skipping it!")

    # create set of indices from the specified index group
    try:
        idxs_parsed = idx[idx_group]["content"]
    except KeyError as e:
        raise KeyError(
            f" Did not find index group '{idx_group}' in index file groups: {list(idx.keys())}"
        ) from e

    idxs_set = set([entry for l in idxs_parsed for entry in l])

    # find bonds with both indices in the index group
    plumed_distance_idxs = []
    for bond in top.bonds.keys():
        if bond[0] in idxs_set and bond[1] in idxs_set:
            ai_name = top.atoms[bond[0]].atom
            aj_name = top.atoms[bond[1]].atom
            if set([ai_name, aj_name]) in bonds_to_remove:
                continue
            if any([b in atoms_to_remove for b in [ai_name, aj_name]]):
                continue
            if any([b[0] in atoms_start_to_remove for b in [ai_name, aj_name]]):
                continue
            plumed_distance_idxs.append(bond)

    # write plumed configuration file
    print_arg = ""
    file = open("plumed.dat", "w")  # open in  append mode
    file.write("#Define distances \n")
    for i, dist in enumerate(plumed_distance_idxs):
        file.write(f"d{i}: DISTANCE ATOMS={dist[0]},{dist[1]}\n")
        print_arg += f"d{i},"
    file.write(" \n#Print distances ARG to FILE every STRIDE steps \n")
    file.write(
        "PRINT ARG="
        + print_arg
        + f" STRIDE={plumed_stride} FILE={plumed_distances_out}"
    )
    file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build examples for KIMMDY.")
    parser.add_argument("topfile", type=str, help="Gromacs topology file path.")
    parser.add_argument("indexfile", type=str, help="Gromacs index file path.")
    parser.add_argument(
        "indexgroup",
        type=str,
        help="Index group name out of which bonds will be written to the plumed configuration file.",
    )
    parser.add_argument(
        "--entries-to-remove",
        nargs="*",
        default=["C-N", "H*", "O*"],
        help="Either atomnames or bonds, i.e. atomnames separated by a hyphen that should not be written to the plumed configuration file.",
    )
    p = parser.parse_args()

    top_path = Path(p.topfile)
    idx_path = Path(p.indexfile)

    create_plumed_dat(
        top_path,
        idx_path,
        idx_group=p.indexgroup,
        entries_to_remove=p.entries_to_remove,
    )
