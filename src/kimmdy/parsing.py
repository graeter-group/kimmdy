import os
from pathlib import Path
from collections.abc import Iterable
from typing import Generator, Optional, Union
from copy import deepcopy
import xml.etree.ElementTree as ET
from itertools import takewhile
from kimmdy.utils import pushd

TopologyDict = dict[str, Union[list[list[str]],dict[str, list[str]]]]


def is_not_comment(c: str) -> bool:
    return c != ";"


def get_sections(
    seq: Iterable[str], section_marker: str
) -> Generator[list[str], None, None]:
    data = [""]
    for line in seq:
        line = "".join(takewhile(is_not_comment, line))
        if line.strip(" ").startswith(section_marker):
            if data:
                # first element will be empty
                # because newlines mark sections
                data.pop(0)
                # only yield section if non-empty
                if len(data) > 0:
                    yield data
                data = [""]
        data.append(line.strip("\n"))
    if data:
        yield data


def extract_section_name(ls: list[str]) -> tuple[str, list[str]]:
    """takes a list of lines and return a tuple
    with the name and the lines minus the
    line that contained the name.
    Returns the empty string as the name if no name was found.
    """
    for i, l in enumerate(ls):
        if l and l[0] != ";" and "[" in l:
            name = l.strip("[] \n")
            ls.pop(i)
            return (name, ls)
    else:
        return ("", ls)


def create_subsections(ls: list[list[str]]):
    d = {}
    subsection_name = "other"
    for _, l in enumerate(ls):
        if l[0] == "[":
            subsection_name = l[1]
        else:
            if subsection_name not in d:
                d[subsection_name] = []
            d[subsection_name].append(l)

    return d


def read_rtp(path: Path) -> dict:
    # TODO: make this more elegant and performant
    with open(path, "r") as f:
        sections = get_sections(f, "\n")
        d = {}
        for i, s in enumerate(sections):
            # skip empty sections
            if s == [""]:
                continue
            name, content = extract_section_name(s)
            content = [c.split() for c in content if len(c.split()) > 0]
            if not name:
                name = f"BLOCK {i}"
            d[name] = create_subsections(content)
            # d[name] = content

        return d

def resolve_includes(path: Path) -> list[str]:
    """Resolve #include statements in a (top/itp) file."""
    dir = path.parent
    fname = path.name
    cwd = Path.cwd()
    os.chdir(dir)
    ls_prime = []
    with open(fname, "r") as f:
        for l in f:
            l = "".join(takewhile(is_not_comment, l)).strip()
            if not l: continue
            if l.startswith("#include"):
                path = Path(l.split('"')[1])
                try:
                    ls_prime.extend(resolve_includes(path))
                except Exception as _:
                    # keep the line if file path can't be resolved
                    # might be a path that only exists for a later step
                    # such as posres.itp created only later
                    ls_prime.append(l)
            else:
                ls_prime.append(l)

    os.chdir(cwd)
    return ls_prime

SECTIONS_WITH_SUBSECTIONS = ("moleculetype",)
NESTABLE_SECTIONS = ("atoms", "bonds", "pairs", "angles", "dihedrals", "impropers", "exclusions", "virtual_sites", "settles", "position_restraints")


def parse_topol(ls: list[str]) -> TopologyDict:
    """Parse a list of lines from a topology file.

    Assumptions:
    - `#include` statements have been resolved
    - comments have been removed
    - empty lines have been removed
    - all lines are stripped of leading and trailing whitespace
    - `#ifdef .. #endif` statements only surround a full section or subsection,
      not individual lines within a section
    - a section within `ifdef` may be a subsection of a section that was started
      outside of the `ifdef`
    """
    d = {}
    d['define'] = {}
    parent_section = None
    subsection_name = None
    section = None
    for l in ls:
        print(l)
        if l.startswith('*'):
            # comments at the start of forcefield.itp
            continue
        elif l.startswith("#define"):
            l = l.split()
            name = l[1]
            values = l[2:]
            d['define'][name] = values
        elif l.startswith("#if"):
            # TODO
            l = l.split()
            condition_type = l[0].removeprefix('#')
            condition = l[1]
        elif l.startswith("#endif"):
            # TODO
            continue
        elif l.startswith("["):
            # start a new section
            section = l.strip("[] \n").lower()
            if section in SECTIONS_WITH_SUBSECTIONS:
                parent_section = section
                section = None
                if d.get(parent_section) is None:
                    d[parent_section] = {
                        'content': [],
                        'extra': [],
                        'subsections': {}
                    }
            else:
                if parent_section is not None:
                    # in a parent_section that can have subsections
                    if section in NESTABLE_SECTIONS:
                        assert subsection_name is not None
                        if d[parent_section]['subsections'].get(subsection_name) is None:
                            d[parent_section]['subsections'][subsection_name] = {}
                        if d[parent_section]['subsections'][subsection_name].get(section) is None:
                            d[parent_section]['subsections'][subsection_name][section] = {
                                'content': [],
                                'extra': []
                            }
                    else:
                        # exit parent_section by setting it to None
                        parent_section = None
                        if d.get(section) is None:
                            d[section] = {
                                'content': [],
                                'extra': []
                            }
                else:
                    if d.get(section) is None:
                        # regular section that is not a subsection
                        d[section] = {
                            'content': [],
                            'extra': []
                        }
        else:
            if parent_section is not None:
                # line is in a section that can have subsections
                if section is None:
                    # but no yet in a subsection
                    l = l.split()
                    d[parent_section]['content'].append(l)
                    # the following subsections will be grouped
                    # under the name of the first line
                    # of the content of the current_section
                    subsection_name = l[0].lower()
                else:
                    print(parent_section)
                    print(subsection_name)
                    print(section)
                    # line is in in a subsection
                    d[parent_section]['subsections'][subsection_name][section]['content'].append(l.split())
            else:
                d[section]['content'].append(l.split())
    return d


def read_topol(path: Path) -> TopologyDict:
    # TODO look into following #includes
    # TODO look into [ intermolecule ] section
    ls = resolve_includes(path)
    return parse_topol(ls)


def write_topol(top: TopologyDict, outfile: Path):
    with open(outfile, "w") as f:
        for title, content in top.items():
            if content == []:
                continue
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


def split_dihedrals(top: TopologyDict):
    if "dihedrals" in top.keys():
        top["propers"] = deepcopy(top["dihedrals"])
        top["impropers"] = []
        for dih in top["propers"][::-1]:
            if dih[4] == "9":
                break
            else:
                top["impropers"].insert(0, (top["propers"].pop(-1)))


def merge_propers_impropers(top: TopologyDict):
    if set(["propers", "impropers"]).issubset(top.keys()):
        top["dihedrals"].clear()
        top["dihedrals"].extend(top.pop("propers"))
        top["dihedrals"].extend(top.pop("impropers"))


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
        d = {c: [] for c in colnames}
        for l in f:
            values = l.split()
            for k, v in zip(colnames, values):
                d[k].append(v)

    return d


def read_plumed_distances(plumed_dat: Path, distances_dat: Path):
    plumed = read_plumed(plumed_dat)
    distances = read_distances_dat(distances_dat)

    atoms = {
        x["id"]: x["atoms"] for x in plumed["distances"] if x["keyword"] == "DISTANCE"
    }

    return atoms


def read_xml_ff(path: Path) -> ET.Element:
    tree = ET.parse(path)
    root = tree.getroot()
    return root
