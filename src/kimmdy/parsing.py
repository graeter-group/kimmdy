"""
All read_<...> and write_<...> functions.
"""
import os
import re
from pathlib import Path
from collections.abc import Iterable
from typing import Generator, Optional, Union
import xml.etree.ElementTree as ET
from itertools import takewhile
import logging
import json
import numpy as np
from typing import TypedDict

from kimmdy.utils import get_gmx_dir

logger = logging.getLogger(__name__)
TopologyDict = dict
"""A raw representation of a topology file.

Every section, apart from `ffdir` and `define`,
comes with a condition that can be checked against the
`define`s by the helper functions to determine if the content
(a list of lists) should come from `content` or `else_content`.
Some sections such as `moleculetype` also come with `subsections`.

Examples
-------
```python
raw_top = 
{'ffdir': PosixPath('/usr/share/gromacs/top/amber99.ff'),
 'define': {'_FF_AMBER': [], '_FF_AMBER99': []},
 'defaults': {'content': [['1', '2', 'yes', '0.5', '0.8333']],
  'else_content': [],
  'extra': [],
  'condition': None},
 'atomtypes': {'content': [
   ['C', '6', '12.01', '0.0000', 'A', '3.39967e-01', '3.59824e-01'],
   ['MNH3', '0', '0.0000', '0.0000', 'A', '0.00000e+00', '0.00000e+00']],
  'else_content': [],
  'extra': [],
  'condition': None},
 'moleculetype_0': {'content': [['Urea', '3']],
  'else_content': [],
  'extra': [],
  'condition': None,
  'subsections': {'atoms': {'content': [['1',
      'C',
      '1',
      'URE',
      'C',
      '1',
      '0.880229',
      '12.01000'],
     ['2', 'O', '1', 'URE', 'O', '2', '-0.613359', '16.00000'],
     ['3', 'N', '1', 'URE', 'N1', '3', '-0.923545', '14.01000'],
     ['4', 'H', '1', 'URE', 'H11', '4', '0.395055', '1.00800'],
     ['5', 'H', '1', 'URE', 'H12', '5', '0.395055', '1.00800'],
     ['6', 'N', '1', 'URE', 'N2', '6', '-0.923545', '14.01000'],
     ['7', 'H', '1', 'URE', 'H21', '7', '0.395055', '1.00800'],
     ['8', 'H', '1', 'URE', 'H22', '8', '0.395055', '1.00800']],
    'else_content': [],
    'extra': [],
    'condition': None},
   'bonds': {'content': [['1', '2'],
     ['1', '3'],
     ['1', '6'],
     ['3', '4'],
     ['3', '5'],
     ['6', '7'],
     ['6', '8']],
    'else_content': [],
    'extra': [],
    'condition': None},
   'dihedrals': {'content': [['2', '1', '3', '4', '9'],
     ['2', '1', '3', '5', '9'],
     ['2', '1', '6', '7', '9'],
     ['2', '1', '6', '8', '9'],
     ['3', '1', '6', '7', '9'],
     ['3', '1', '6', '8', '9'],
     ['6', '1', '3', '4', '9'],
     ['6', '1', '3', '5', '9'],
     ['3', '6', '1', '2', '4'],
     ['1', '4', '3', '5', '4'],
     ['1', '7', '6', '8', '4']],
    'else_content': [],
    'extra': [],
    'condition': None},
   'position_restraints': {'content': [['1', '1', '1000', '1000', '1000'],
     ['2', '1', '1000', '0', '1000'],
     ['3', '1', '1000', '0', '0']],
    'else_content': [],
    'extra': [],
    'condition': None},
   'dihedral_restraints': {'content': [['3',
      '6',
      '1',
      '2',
      '1',
      '180',
      '0',
      '10'],
     ['1', '4', '3', '5', '1', '180', '0', '10']],
    'else_content': [],
    'extra': [],
    'condition': None}}},
}
```
"""


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

        return d


def resolve_includes(
    path: Path, gmx_builtin_ffs: Optional[Path] = None
) -> tuple[list[str], Optional[Path]]:
    """Resolve #include statements in a (top/itp) file.

    Arguments
    ---------
    path :
        Filepath to read.

    Returns
    -------
    lines:
        List of lines.
    ffdir:
        Path to the ff directory if one of the includes used a file from it.
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
                if ".ff" in l:
                    # e.g. #include "amber99.ff/forcefield.itp"
                    ffdir = include_path.parent
                    # test if the path is in the cwd, otherwise search in gmx ff dir
                    if not ffdir.exists() and gmx_builtin_ffs is not None:
                        ffdir = gmx_builtin_ffs / ffdir
                    ffdir = ffdir.resolve()
                ls_prime, _ = resolve_includes(include_path, gmx_builtin_ffs)
                if not ls_prime and gmx_builtin_ffs is not None:
                    ls_prime, _ = resolve_includes(
                        gmx_builtin_ffs / include_path, gmx_builtin_ffs
                    )
                if not ls_prime:
                    logger.warning(
                        f"top include {include_path} could not be resolved. Line was dropped."
                    )
                ls.extend(ls_prime)
            else:
                ls.append(l)

    os.chdir(cwd)
    return (ls, ffdir)


def read_top(
    path: Path, ffdir: Optional[Path] = None, use_gmx_dir: bool = True
) -> TopologyDict:
    """Read a topology file into a raw TopologyDict represenation.

    Parameters
    ----------
    path :
        Path to the topology file.

    Assumptions and limitation
    -----
    - `#include` statements will be resolved
    - comments will be removed
    - all lines are stripped of leading and trailing whitespace
    - `#undef` is not supported
    - a section within `ifdef` may be a subsection of a section that was started
      outside of the `ifdef`
    - `#if..#endif` statements only surround a full section or subsection,
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
        "pairs_nb",
        "exclusions",
        "constraints" "virtual_sites",
        "virtual_sites1",
        "virtual_sites2",
        "virtual_sites3",
        "virtual_sites4",
        "virtual_sitesn",
        "settles",
        "position_restraints",
        "dihedral_restraints",
        "angle_restraints",
        "angle_restraints_z",
        "orientation_restraints",
    )

    gmx_builtin_ffs = None
    if use_gmx_dir:
        gmxdir = get_gmx_dir()
        if gmxdir is not None:
            gmx_builtin_ffs = gmxdir / "top"
    ls, parsed_ffdir = resolve_includes(path, gmx_builtin_ffs)
    if ffdir is None and parsed_ffdir is not None:
        ffdir = parsed_ffdir
    if ffdir is None:
        logger.warning(f"No #include for a forcefield directory found in {path}.")

    ls = filter(lambda l: not l.startswith("*"), ls)
    d = {}
    d["ffdir"] = ffdir
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
                    f"""Error parsing {path}: #if ... #endif can only be used to surround a section, not within."""
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

    if len(d) <= 2:
        raise ValueError(f"topology file {path} does not contain any sections")

    return d


def write_top(top: TopologyDict, outfile: Path):
    """Write a TopologyDict to a topology file.

    Parameters
    ----------
    top:
        Raw dictionary represenation of the topology.
    outfile :
        Path to the topology file to write to.
    """

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
        if define is not None:
            for name, value in define.items():
                f.writelines("#define " + name + " ".join(value))
                f.write("\n")

        for name, section in top.items():
            if name in ["define", "ffdir"]:
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

    Parameters
    ----------
    path :
        Path to the file. E.g. Path("edissoc.dat")
    """
    with open(path, "r") as f:
        edissocs = {}
        for l in f:
            at1, at2, edissoc, *_ = l.split()
            edissocs[frozenset((at1, at2))] = float(edissoc)
    return edissocs


class Plumed_dict(TypedDict):
    """Dict representation of a plumed.dat file."""

    distances: list[dict]
    prints: list[dict]


def read_plumed(path: Path) -> Plumed_dict:
    """Read a plumed.dat configuration file.

    Parameters
    ----------
    path :
        Path to the file. E.g. "plumed.dat"

    Returns
    -------
    dict :
        dict with keys: 'distances' and 'prints'
        Each is a list of dicts containing plumed keywords
    """
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
                        "atoms": [str(x) for x in d[2].strip("ATOMS=").split(",")],
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
                # TODO: output should be consolidated dicts instad of lists of dicts
                d["ARG"] = d["ARG"].split(",")
                d["STRIDE"] = int(d["STRIDE"])
                d["FILE"] = Path(d["FILE"])

                prints.append(d)

        return Plumed_dict(distances=distances, prints=prints)


def write_plumed(d: Plumed_dict, path: Path) -> None:
    """Write a plumed.dat configuration file.

    Parameters
    ----------
    d :
        Dictionary containing 'distances' and 'prints'
    path :
        Path to the file. E.g. "plumed.dat"
    """
    with open(path, "w") as f:
        for l in d["distances"]:
            f.write(
                f"{l['id']}: {l['keyword']} ATOMS={l['atoms'][0]},{l['atoms'][1]} \n"
            )
        f.write("\n")
        for l in d["prints"]:
            f.write(
                f"{l['PRINT']} ARG={','.join(l['ARG'])} STRIDE={str(l['STRIDE'])} FILE={str(l['FILE'].name)}\n"
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


class JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(JSONEncoder, self).default(obj)


def write_json(d: dict, path: Path) -> None:
    logger.debug(f"writing dictionary to json: {d}")
    with open(path, "w") as f:
        json.dump(d, f, cls=JSONEncoder)


def read_json(path: Union[str, Path]) -> dict:
    with open(path, "r") as f:
        data = json.load(f)
    return data
