"""
All read_<...> and write_<...> functions.
"""

import csv
import json
import logging
import os
from datetime import datetime
from itertools import takewhile
from pathlib import Path
from typing import Optional, TypedDict, Union

import numpy as np

from kimmdy.utils import get_gmx_dir

logger = logging.getLogger(__name__)
TopologyDict = dict
"""A raw representation of a topology file returned by [](`~kimmdy.parsing.read_top`).
"""


## convenience functions
def is_not_comment(c: str) -> bool:
    """Returns whether a string is not a comment.

    Used for topology like files that use ';' for comments.
    """
    return c != ";"


def resolve_includes(
    path: Path, gmx_builtin_ffs: Optional[Path] = None
) -> tuple[list[str], Optional[Path]]:
    """Resolve #include statements in a (top/itp) file.

    Parameters
    ----------
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
    cwd = Path.cwd()
    if not dir.exists() or not path.exists():
        return ([], None)
    os.chdir(dir)
    ls = []
    ffdir = None

    with open(path, "r") as f:
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
                    logger.info(
                        f"top include {include_path} could not be resolved. Line was dropped."
                    )
                ls.extend(ls_prime)
            else:
                ls.append(l)

    os.chdir(cwd)
    return (ls, ffdir)


## parsing functions ##
## topology
def empty_section(condition: Optional[dict] = None) -> dict:
    return {"content": [], "else_content": [], "extra": [], "condition": condition}


def read_top(
    path: Path,
    ffdir: Optional[Path] = None,
    use_gmx_dir: bool = True,
) -> TopologyDict:
    """Read a topology file (*.top,*.itp,*.rtp) into a raw TopologyDict represenation.

    Assumptions and limitation:

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

    Parameters
    ----------
    path :
        Path to the topology file.

    Returns
    ------
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
    'moleculetype_Urea': {'content': [['Urea', '3']],
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

    CHILD_SECTIONS = (
        "atoms",
        "bonds",
        "pairs",
        "angles",
        "dihedrals",
        "impropers",
        "pairs_nb",
        "exclusions",
        "constraints",
        "virtual_sites",
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
        logger.debug(f"No #include for a forcefield directory found in {path}.")

    ls = [l for l in ls if not l.startswith("*")]
    d = {}
    d["ffdir"] = ffdir
    d["define"] = {}
    parent_section = None
    section = None
    condition = None
    is_first_line_after_section_header = False
    content_key = "content"

    for i, l in enumerate(ls):
        # decide where to put lines depending on current context
        # deal with statements
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
        # deal with a new section
        elif l.startswith("["):
            section = l.strip("[] \n")
            # add nested structure for non-nestable sections
            if section not in CHILD_SECTIONS:
                parent_section = section
                section = None
                if d.get(parent_section) is None:
                    if parent_section == "moleculetype":
                        # append the name of the moleculetype to the section name
                        # but the name comes on the next line
                        moleculetype_name = ls[i + 1].split()[0]
                        parent_section = f"{parent_section}_{moleculetype_name}"
                    d[parent_section] = empty_section(condition)
                    d[parent_section]["subsections"] = {}
            else:
                if parent_section is not None:
                    # add empty section as subsection
                    if d[parent_section]["subsections"].get(section) is None:
                        d[parent_section]["subsections"][section] = empty_section(
                            condition
                        )
                else:
                    # add regular section
                    if d.get(section) is None:
                        d[section] = empty_section(condition)
            is_first_line_after_section_header = True
        # deal with content lines
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
                # adding content to dict here; this is where the loop goes 98% of the time
                d[section][content_key].append(l.split())
            is_first_line_after_section_header = False

    if len(d) <= 2:
        m = f"topology file {path} does not contain any sections"
        logger.error(m)
        raise ValueError(m)

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

    def write_if_nonempty(f, printname, ls):
        if ls:
            f.write(f"[ {printname} ]\n")
            for l in ls:
                f.writelines(" ".join(l))
                f.write("\n")

    def write_section(f, name: str, section):
        if name.startswith("moleculetype_"):
            printname = "moleculetype"
        else:
            printname = name
        condition = section.get("condition")
        f.write("\n")
        if condition is None:
            write_if_nonempty(f, printname, section["content"])
        else:
            condition_type = condition["type"]
            condition_value = condition["value"]
            f.write(f"#{condition_type} {condition_value}")
            f.write("\n")
            write_if_nonempty(f, printname, section["content"])
            else_content = section.get("else_content")
            if else_content:
                f.write("#else\n")
                write_if_nonempty(f, printname, else_content)
            f.write("#endif\n")

    with open(outfile, "w") as f:
        # write defines first
        define = top.get("define")
        if define is not None:
            for name, value in define.items():
                f.writelines("#define " + name + " ".join(value))
                f.write("\n")

        for name, section in top.items():
            if name in ["define", "ffdir", "system", "molecules"]:
                continue
            f.write("\n")
            subsections = section.get("subsections")
            write_section(f, name, section)
            if subsections is not None:
                for name, section in subsections.items():
                    write_section(f, name, section)

        # write system and molecules last
        for name in ["system", "molecules"]:
            section = top.get(name)
            if section is not None:
                f.write("\n")
                write_section(f, name, section)


## Plumed
class Plumed_dict(TypedDict):
    """Dict representation of a plumed.dat file."""

    other: list
    labeled_action: dict
    prints: list[dict]


def read_plumed(path: Path) -> Plumed_dict:
    """Read a plumed.dat configuration file.

    Follows the plumed naming scheme of label, keyword, action.

    Parameters
    ----------
    path :
        Path to the file. E.g. "plumed.dat"

    Returns
    -------
    dict :
        dict with keys: 'distances' and 'prints'
        Each is a dict/list of dicts containing plumed keywords
    """
    with open(path, "r") as f:
        other = []
        labeled_action = {}
        prints = []
        for l in f:
            # deal with PRINT action
            if l.startswith("PRINT"):
                l = l.split()
                d = {}
                d["PRINT"] = l[0]
                for x in l[1:]:
                    key, value = x.split("=")
                    d[key] = value

                # STRIDE is a mandatory keyword, others optional
                d["STRIDE"] = int(d["STRIDE"])
                if d.get("ARG"):
                    d["ARG"] = d["ARG"].split(",")
                if d.get("FILE"):
                    d["FILE"] = Path(d["FILE"])

                prints.append(d)
            # deal with labeled action
            elif ": " in l:
                ll = l.split(sep=":")

                action = ll[1].split()
                if action[0] == "DISTANCE":
                    labeled_action[ll[0]] = {
                        "keyword": action[0],
                        "atoms": [str(x) for x in action[1].strip("ATOMS=").split(",")],
                    }
                else:
                    labeled_action[ll[0]] = {"action": action}
            # deal with other lines in config
            else:
                if not any(l.strip().startswith(x) for x in ["#", "\n"]) and l.strip():
                    other.append(l)

        return Plumed_dict(other=other, labeled_action=labeled_action, prints=prints)


def write_plumed(d: Plumed_dict, path: Path) -> None:
    """Write a plumed.dat configuration file.

    Parameters
    ----------
    d :
        Dictionary containing 'labeled_action', 'other' and 'prints'
    path :
        Path to the file. E.g. "plumed.dat"
    """
    with open(path, "w") as f:
        for l in d["other"]:
            f.write(f"{l}\n")
        for k, v in d["labeled_action"].items():
            if v.get("keyword") and v.get("atoms"):
                f.write(f"{k}: {v['keyword']} ATOMS={','.join(v['atoms'])} \n")
            elif v.get("action"):
                f.write(f"{k}:{' '.join(v['action'])}\n")
            else:
                logger.debug(f"Skipping plumed config file line '{k}:{v}'\n")
        f.write("\n")
        for l in d["prints"]:
            f.write(
                f"{l['PRINT']} ARG={','.join(l['ARG'])} STRIDE={str(l['STRIDE'])} FILE={str(l['FILE'].name)}\n"
            )


def read_distances_dat(distances_dat: Path) -> dict:
    """Read a distances.dat plumed output file."""
    with open(distances_dat, "r") as f:
        colnames = f.readline()[10:].split()
        d = {c: [] for c in colnames}
        for l in f:
            values = l.split()
            for k, v in zip(colnames, values):
                d[k].append(float(v))

    return d


## JSON
class JSONEncoder(json.JSONEncoder):
    """Encoder that enables writing JSONs with numpy types."""

    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        elif isinstance(o, np.floating):
            return float(o)
        elif isinstance(o, np.ndarray):
            return o.tolist()
        else:
            return super(JSONEncoder, self).default(o)


def write_json(
    d: dict,
    path: Path,
) -> None:
    """Write dict to file according to JSON format."""
    logger.debug(f"writing dictionary to json: {path}")
    with open(path, "w") as f:
        json.dump(d, f, cls=JSONEncoder, indent=4)


def read_json(path: Union[str, Path]) -> dict:
    """Return JSON file content as dict."""
    with open(path, "r") as f:
        data = json.load(f)
    return data


def read_csv_to_list(csv_file: Path) -> list:
    data = []
    with open(csv_file, "r") as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            data.extend(row)
    return data


## Miscellaneous files
def read_edissoc(path: Path) -> dict:
    """Reads a edissoc file and turns it into a dict.

    The dissociation energy is assigned per pair of atom names. Atom names are unique to a residue, and the dict is nested by residues.
    The set of bond atoms make up the key, the dissociation energy E_dissoc [kJ mol-1] is the value.


    Parameters
    ----------
    path :
        Path to the file. E.g. Path("edissoc.dat")
    """
    with open(path, "r") as f:
        key = "_"
        edissocs = {key: {}}
        for l in f:
            if l.startswith(";"):
                continue
            elif l.strip().startswith("[") and l.strip().endswith("]"):
                key = l.strip().strip("[]")
                if key != key.strip():
                    m = f"Unexpected whitespace in edissoc file {path} for key: {l}. Please remove it and try again."
                    logger.error(m)
                    raise ValueError(m)
                edissocs[key] = {}
            elif len(l.split(sep=";")[0].split()) == 3:
                at1, at2, edissoc, *_ = l.split()
                edissocs[key][frozenset([at1, at2])] = float(edissoc)
            else:
                logger.debug(f"Unexpected line in edissoc file: {l}")
    return edissocs


def write_time_marker(p: Path, marker: str):
    with open(p, "a") as f:
        f.write(f"{marker},{datetime.now().isoformat()}\n")


def read_time_marker(p: Path) -> tuple[list[str], list[datetime]]:
    with open(p) as f:
        lines = f.readlines()
    events = []
    times = []
    for line in lines:
        e, t = line.strip().split(",")
        try:
            events.append(e)
            times.append(datetime.fromisoformat(t))
        except ValueError as err:
            logger.error(f"Error trying to read time for event {e}: {err}\nfile: {p}")
    return events, times
