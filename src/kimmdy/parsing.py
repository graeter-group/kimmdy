import os
from pathlib import Path
from collections.abc import Iterable
from typing import Generator, Optional, Union
from copy import deepcopy
import xml.etree.ElementTree as ET
from itertools import takewhile
from kimmdy.utils import pushd
from typing import TYPE_CHECKING
import logging

TopologyDict = dict[str, dict[str, list[str]]]


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
                    # drop line if path can't be resolved
                    logging.warn('top include could not be resolved.')
                    continue
            else:
                ls_prime.append(l)

    os.chdir(cwd)
    return ls_prime

SECTIONS_WITH_SUBSECTIONS = ("moleculetype",)
NESTABLE_SECTIONS = ("atoms", "bonds", "pairs", "angles", "dihedrals", "impropers", "exclusions", "virtual_sites", "settles", "position_restraints")


def parse_topol(ls: Iterable[str]) -> TopologyDict:
    d = {}
    d['define'] = {}
    parent_section_index = 0
    parent_section = None
    subsection_name = None
    section = None
    condition = None

    def empty_section(condition):
        return {
            'content': [],
            'extra': [],
            'condition': condition
        }

    for i,l in enumerate(ls):
        if l.startswith("#define"):
            l = l.split()
            name = l[1]
            values = l[2:]
            d['define'][name] = values
        elif l.startswith("#if"):
            l = l.split()
            condition_type = l[0].removeprefix('#')
            condition_value = l[1]
            condition = {
                'type': condition_type,
                'value': condition_value
            }
        elif l.startswith("#endif"):
            condition = None
            continue
        elif l.startswith("["):
            # start a new section
            section = l.strip("[] \n").lower()
            if section in SECTIONS_WITH_SUBSECTIONS:
                parent_section = section
                section = None
                if d.get(parent_section) is None:
                    parent_section = f"{parent_section}_{parent_section_index}"
                    parent_section_index += 1
                    d[parent_section] = {
                        'content': [],
                        'extra': [],
                        'subsections': {},
                        'condition': condition
                    }
            else:
                if parent_section is not None:
                    # in a parent_section that can have subsections
                    if section in NESTABLE_SECTIONS:
                        assert subsection_name is not None
                        if d[parent_section]['subsections'].get(subsection_name) is None:
                            d[parent_section]['subsections'][subsection_name] = {}
                        if d[parent_section]['subsections'][subsection_name].get(section) is None:
                            d[parent_section]['subsections'][subsection_name][section] = empty_section(condition)
                    else:
                        # exit parent_section by setting it to None
                        parent_section = None
                        if d.get(section) is None:
                            d[section] = empty_section(condition)
                else:
                    if d.get(section) is None:
                        # regular section that is not a subsection
                        d[section] = empty_section(condition)
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
                    # line is in in a subsection
                    d[parent_section]['subsections'][subsection_name][section]['content'].append(l.split())
            else:
                d[section]['content'].append(l.split())
    return d


def read_topol(path: Path) -> TopologyDict:
    """Parse a list of lines from a topology file.

    Assumptions and limitation:
    - `#include` statements have been resolved
    - comments have been removed
    - empty lines have been removed
    - all lines are stripped of leading and trailing whitespace
    - `#if .. #endif` statements only surround a full section or subsection,
      not individual lines within a section
    - `#undef` is not supported
    - a section within `ifdef` may be a subsection of a section that was started
      outside of the `ifdef`
    - a section may either be contained within if ... else or it may not be,
      but it can not be duplicated with one part inside and one outside.
    - `if .. else` can't be nested
    - `#include`s that don't resolve to a valid file path are silently dropped 
    - sections that can have subsections can also exist multiple, separate times
      e.g. moleculetype will appear multiple times and they should not be merged
    """
    ls = resolve_includes(path)
    return parse_topol(filter(lambda l: not l.startswith('*'), ls))


def read_edissoc(path: Path) -> dict:
    """reads a edissoc file and turns it into a dict.
    the tuple of bond atoms make up the key,
    the dissociation energy E_dissoc [kJ mol-1] is the value
    """
    with open(path, "r") as f:
        edissocs = {}
        for l in f:
            at1, at2, edissoc, *_ = l.split()
            edissocs[frozenset((at1, at2))] = float(edissoc)
    return edissocs


def write_topol(top: TopologyDict, outfile: Path):
    with open(outfile, "w") as f:
        # extract sections that have to be written first
        define = top.pop('define')
        # extract sections that have to be written last
        system = top.pop('system')
        molecules = top.pop('molecules')
        if system is None or molecules is None:
            raise ValueError("Invalid topology, no [ system ] or no [ molecules ] section found")

        f.write('\n')
        for name, value in define.items():
            f.writelines('#define ' + name + ' '.join(value))
            f.write('\n')


        for name, section in top.items():
            subsections = section.get('subsections')
            if section['condition'] is None:
                f.write('\n')
                f.write(f"[ {name} ]\n")
                for l in section['content']:
                    f.writelines(' '.join(l))
                    f.write('\n')
            if subsections is not None:
                subsections = section.get('subsections')
                for subsection_name, subsection in subsections.items():
                    f.write('\n')
                    for name, section in subsection.items():
                        condition = section.get('condition')
                        f.write('\n')
                        if condition is not None:
                            condition_type = condition['type']
                            condition_value = condition['value']
                            f.write(f"#{condition_type} {condition_value}")
                            f.write("\n")
                        f.write(f"[ {name} ]\n")
                        for l in section['content']:
                            f.writelines(' '.join(l))
                            f.write('\n')
                        if condition is not None:
                            f.write(f"#endif\n")



        for name, section in [('system', system), ('molecules', molecules)]:
            f.write('\n')
            f.write(f"[ {name} ]\n")
            for l in section['content']:
                f.writelines(' '.join(l))
                f.write('\n')



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
                        "atoms": [int(x) for x in d[2].strip("ATOMS=").split(",")],
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
            f.write(
                f"{l['id']}: {l['keyword']} ATOMS={l['atoms'][0]},{l['atoms'][1]} \n"
            )
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
                d[k].append(float(v))

    return d


# def read_plumed_distances(plumed_dat: Path, distances_dat: Path):
#     plumed = read_plumed(plumed_dat)
#     distances = read_distances_dat(distances_dat)

#     atoms = {
#         x["id"]: x["atoms"] for x in plumed["distances"] if x["keyword"] == "DISTANCE"
#     }

#     return atoms


def read_xml_ff(path: Path) -> ET.Element:
    tree = ET.parse(path)
    root = tree.getroot()
    return root
