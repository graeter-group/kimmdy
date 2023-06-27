import os
import re
from pathlib import Path
from collections.abc import Iterable
from typing import Generator, Optional
import xml.etree.ElementTree as ET
from itertools import takewhile
import logging

from kimmdy.utils import get_gmx_dir

TopologyDict = dict
GMX_BUILTIN_FF_DIR = get_gmx_dir() / "top"


def is_not_comment(c: str) -> bool:
    return c != ";"


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
    # TODO combine with top parser?
    # would need a way to tell the parser
    # that here, all sections have subsections
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


def resolve_includes(path: Path) -> tuple[list[str], Optional[Path]]:
    """Resolve #include statements in a (top/itp) file.
    Returns a tuple with the list of lines and the Path to the ff directory.
    """
    path = path.resolve()
    dir = path.parent
    fname = path.name
    cwd = Path.cwd()
    if not dir.exists() or not path.exists():
        return ([], None)
    os.chdir(dir)
    ls = []
    ffdir = None
    with open(fname, "r") as f:
        for l in f:
            l = "".join(takewhile(is_not_comment, l)).strip()
            if not l:
                continue
            if l.startswith("#include"):
                include_path = Path(l.split('"')[1])
                # if the include path contains `ff` it is a file in a force field directory
                # get the path to the force field directory
                if '.ff' in l:
                    # e.g. #include "amber99.ff/forcefield.itp"
                    ffdir = include_path.parent
                    # test if the path is in the cwd, otherwise search in gmx ff dir
                    if not ffdir.exists():
                        ffdir = GMX_BUILTIN_FF_DIR / ffdir
                    ffdir = ffdir.resolve()
                ls_prime, _ = resolve_includes(include_path)
                if not ls_prime:
                    ls_prime, _ = resolve_includes(GMX_BUILTIN_FF_DIR / include_path)
                if not ls_prime:
                    logging.warning(
                        f"top include {include_path} could not be resolved. Line was dropped."
                    )
                ls.extend(ls_prime)
            else:
                ls.append(l)

    os.chdir(cwd)
    return (ls, ffdir)


def read_top(path: Path) -> TopologyDict:
    """Parse a list of lines from a topology file.

    Assumptions and limitation
    -----
    - ``#include`` statements will be resolved
    - comments will be removed
    - all lines are stripped of leading and trailing whitespace
    - ``#undef`` is not supported
    - a section within ``ifdef`` may be a subsection of a section that was started
      outside of the ``ifdef``
    - ``#if..#endif`` statements only surround a full section or subsection,
      not individual lines within a section and
      a section may either be contained within if ... else or it may not be,
      but it can not be duplicated with one part inside and one outside.
    - `if .. else` can't be nested
    - ``#include`` s that don't resolve to a valid file path are silently dropped
    - sections that can have subsections can also exist multiple, separate times
      e.g. moleculetype will appear multiple times and they should not be merged
    """
    SECTIONS_WITH_SUBSECTIONS = ("moleculetype",)
    NESTABLE_SECTIONS = (
        "atoms",
        "bonds",
        "pairs",
        "angles",
        "dihedrals",
        "impropers",
        "exclusions",
        "virtual_sites",
        "settles",
        "position_restraints",
        "dihedral_restraints",
    )

    ls, ffdir = resolve_includes(path)
    ls = filter(lambda l: not l.startswith("*"), ls)
    d = {}
    d['ffdir'] = ffdir
    d["define"] = {}
    parent_section_index = 0
    parent_section = None
    section = None
    condition = None
    is_first_line_after_section_header = False
    content_key = "content"

    def empty_section(condition):
        return {"content": [], "else_content": [], "extra": [], "condition": condition}

    for l in ls:
        # where to put lines dependign on current context
        if l.startswith("#define"):
            l = l.split()
            name = l[1]
            values = l[2:]
            d["define"][name] = values
            continue
        elif l.startswith("#if"):
            if is_first_line_after_section_header:
                raise NotImplementedError(
                    f"""Errro parsing {path}: #if ... #endif can only be used to surround a section, not within."""
                )
            l = l.split()
            condition_type = l[0].removeprefix("#")
            condition_value = l[1]
            condition = {"type": condition_type, "value": condition_value}
            continue
        elif l.startswith("#else"):
            content_key = "else_content"
            continue
        elif l.startswith("#endif"):
            condition = None
            content_key = "content"
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
                    d[parent_section] = empty_section(condition)
                    d[parent_section]["subsections"] = {}
            else:
                if parent_section is not None:
                    # in a parent_section that can have subsections
                    if section in NESTABLE_SECTIONS:
                        if d[parent_section]["subsections"].get(section) is None:
                            d[parent_section]["subsections"][section] = empty_section(
                                condition
                            )
                    else:
                        # exit parent_section by setting it to None
                        parent_section = None
                        if d.get(section) is None:
                            d[section] = empty_section(condition)
                else:
                    if d.get(section) is None:
                        # regular section that is not a subsection
                        d[section] = empty_section(condition)
            is_first_line_after_section_header = True
        else:
            if parent_section is not None:
                # line is in a section that can have subsections
                if section is None:
                    # but no yet in a subsection
                    l = l.split()
                    d[parent_section][content_key].append(l)
                else:
                    # part of a subsection
                    d[parent_section]["subsections"][section][content_key].append(
                        l.split()
                    )
            else:
                if section is None:
                    raise ValueError(
                        f"topology file {path} contains lines outside of a section"
                    )
                d[section][content_key].append(l.split())
            is_first_line_after_section_header = False

    if len(d) <= 1:
        raise ValueError(f"topology file {path} does not contain any sections")

    return d


def write_top(top: TopologyDict, outfile: Path):
    def write_section(f, name, section):
        printname = re.sub(r"_\d+", "", name)
        condition = section.get("condition")
        f.write("\n")
        if condition is None:
            f.write(f"[ {printname} ]\n")
            for l in section["content"]:
                f.writelines(" ".join(l))
                f.write("\n")
        else:
            condition_type = condition["type"]
            condition_value = condition["value"]
            f.write(f"#{condition_type} {condition_value}")
            f.write("\n")
            f.write(f"[ {printname} ]\n")
            for l in section["content"]:
                f.writelines(" ".join(l))
                f.write("\n")
            else_content = section.get("else_content")
            if else_content:
                f.write("#else\n")
                f.write(f"[ {printname} ]\n")
                for l in else_content:
                    f.writelines(" ".join(l))
                    f.write("\n")
            f.write("#endif\n")

    with open(outfile, "w") as f:
        define = top.get("define")
        for name, value in define.items():
            f.writelines("#define " + name + " ".join(value))
            f.write("\n")

        for name, section in top.items():
            if name == "define":
                continue
            f.write("\n")
            subsections = section.get("subsections")
            write_section(f, name, section)
            if subsections is not None:
                for name, section in subsections.items():
                    write_section(f, name, section)


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


def read_xml_ff(path: Path) -> ET.Element:
    tree = ET.parse(path)
    root = tree.getroot()
    return root
