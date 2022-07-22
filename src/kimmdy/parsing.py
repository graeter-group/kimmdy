from pathlib import Path
from collections.abc import Iterable
from typing import Generator
import pandas as pd
from copy import deepcopy

# time to experiment
# what a topplogy is
Topology = dict[str, list[list[str]]]


def get_sections(
    seq: Iterable[str], section_marker: str
) -> Generator[list[str], None, None]:
    data = [""]
    for line in seq:
        if line.startswith(section_marker):
            if data:
                # first element will be empty
                # because newlines mark sections
                data.pop(0)
                # only yield section if non-empty
                if data:
                    yield data
                data = []
        data.append(line.strip("\n"))
    if data:
        yield data


def extract_section_name(ls: list[str]) -> tuple[str, list[str]]:
    """takes a list of lines and return a tuple
    with the name and the lines minus the
    line that contained the name.
    Returns the empty string of no name was found.
    """
    for i, l in enumerate(ls):
        if l and l[0] != ";" and "[" in l:
            name = l.strip("[] \n")
            ls.pop(i)
            return (name, ls)
    else:
        return ("", ls)


def read_topol(path: Path) -> Topology:
    # TODO look into following #includes
    # TODO look into [ intermolecule ] section
    with open(path, "r") as f:
        sections = get_sections(f, "\n")
        d = {}
        for i, s in enumerate(sections):
            # skip empty sections
            if s == [""]:
                continue
            name, content = extract_section_name(s)
            content = [c.split() for c in content if c]
            if not name:
                name = f"BLOCK {i}"
            # sections can be duplicated.
            # append values in such cases
            if name not in d:
                d[name] = content
            else:
                d[name] += content
        return d


def write_topol(d: Topology, outfile: Path) -> None:
    with open(outfile, "w") as f:
        for title, content in d.items():
            if title.startswith("BLOCK "):
                f.write(f"\n")
            else:
                f.write(f"[ {title} ]\n")
            s = (
                "\n".join(
                    [
                        " ".join([x.ljust(8) if l[0] != ";" else x for x in l])
                        for l in content
                    ]
                )
                + "\n\n"
            )
            f.write(s)

def topol_split_dihedrals(d: Topology):
    if 'dihedrals' in d.keys():
        d['propers'] = deepcopy(d['dihedrals'])
        d['impropers'] = []
        for i, dih in enumerate(d['propers'][::-1]):
            if dih[4] == '9':
                break
            else:
                d['impropers'].insert(0,(d['propers'].pop(-1)))
    return d

def topol_merge_propers_impropers(d: Topology):
    if all([x in d.keys() for x in ['propers','impropers']]):
        d['dihedrals'].clear()
        d['dihedrals'].extend(d.pop('propers'))
        d['dihedrals'].extend(d.pop('impropers'))
    return d



def read_plumed(path: Path) -> dict:
    """Read a plumed.dat configuration file."""
    with open(path, "r") as f:
        distances = []
        prints = []
        for l in f:
            if "DISTANCE" in l:
                d = l.split()
                distances.append(
                    {
                        "id": d[0].strip(":"),
                        "keyword": d[1],
                        "atoms": d[2].strip("ATOMS=").split(","),
                    }
                )
            elif "PRINT" in l[:5]:
                l = l.split()
                d = {}
                d["PRINT"] = l[0]
                for x in l[1:]:
                    key, value = x.split("=")
                    d[key] = value

                # TODO this could be better. the keys may not exist
                # and break the program here.
                d["ARG"] = d["ARG"].split(",")
                d["STRIDE"] = int(d["STRIDE"])
                d["FILE"] = Path(d["FILE"])

                prints.append(d)

        return {"distances": distances, "prints": prints}


def write_plumed(d, path: Path) -> None:
    """Write a plumed.dat configuration file."""
    with open(path, "w") as f:
        for l in d["distances"]:
            f.write(f"{l['id']}: {l['keyword']} ATOMS={','.join(l['atoms'])} \n")
        f.write("\n")
        for l in d["prints"]:
            f.write(
                f"{l['PRINT']} ARG={','.join(l['ARG'])} STRIDE={str(l['STRIDE'])} FILE={str(l['FILE'])}\n"
            )


def read_distances_dat(distances_dat: Path):
    """Read a distances.dat plumed output file."""
    with open(distances_dat, "r") as f:
        colnames = f.readline()[10:].split()

    return pd.read_csv(distances_dat, delim_whitespace=True, skiprows=1, names=colnames)


def read_plumed_distances(plumed_dat: Path, distances_dat: Path):
    plumed = read_plumed(plumed_dat)
    # plumed['distances']: [{'id': 'd0', 'keyword': 'DISTANCE', 'atoms': ['5', '7']}

    distances = read_distances_dat(distances_dat)
    # distances is a pd DataFrame with time and id's -> distance

    atoms = {
        x["id"]: x["atoms"] for x in plumed["distances"] if x["keyword"] == "DISTANCE"
    }

    return atoms
