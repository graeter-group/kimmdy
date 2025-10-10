"""
All read_<...> and write_<...> functions.
"""

import csv
import json
import logging
from datetime import datetime
from itertools import takewhile
from pathlib import Path
from typing import TypeAlias, TypedDict, Union
import numpy as np

from gmx_top4py.parsing import (
    TopologyDict,
    is_not_comment,
    resolve_includes,
    empty_section,
    read_top,
    write_top,
)

logger = logging.getLogger(__name__)


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


def read_distances_dat(path: Path, dt: float = 0) -> dict[float, dict[str, float]]:
    """Read a distances.dat plumed output file.

    A typical file looks like this:

    ```
    #! FIELDS time d0 d1 d2 d3 d4 d5 d6  ...
    0.000000 0.153211 0.157662 0.139923 ...
    ```
    """
    with open(path, "r") as f:
        colnames = f.readline()[10:].strip().split()
        d = {}
        for i, l in enumerate(f):
            if "#" in l:
                i = l.find("#")
                logger.warning(
                    f"Found second header in plumed file {path.name} in {path.parent.name}. Ignoring the rest of the line."
                )
                continue

            l = l.strip().split()
            time = l[0]
            # time is in ps
            # but needs to be truncated to
            # 3 decimal places (1fs) to avoid
            # floating point errors
            time = round(float(time), 3)

            # if we went back in time that's because PLUMED started from earlier during a restart
            # we want to keep only the latest data
            # this is achieved by using a dictionary indexed by time

            # if the time is not a multiple of dt, skip it
            if dt != 0 and round(i % dt, 3) != 0:
                continue

            # iterate over the rest of the columns
            d[time] = {}
            for k, v in zip(colnames[1:], l[1:]):
                d[time][k] = float(v)

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


EdissocDict: TypeAlias = dict[str, dict[tuple[str, str], float]]


## Miscellaneous files
def read_edissoc(path: Path) -> EdissocDict:
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
                interaction_key = tuple(sorted([at1, at2]))
                edissocs[key][interaction_key] = float(edissoc)
            elif len(l.strip()) > 0:
                logger.warning(f"Unexpected line in edissoc file: {l}")
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


def read_mdp(p: Path) -> dict[str, str]:
    """Reads a mdp file and returns a dict.

    MDP files are key-value pairs separated by '='.
    ; denotes a comment, the rest of the line is ignored.
    whitespace is ignored.
    """
    d = {}
    with open(p) as f:
        for l in f:
            l = "".join(takewhile(is_not_comment, l)).strip()
            if not l:
                continue
            k, v = l.split("=")
            v = v.strip()
            if v.lower() in ["yes", "true"]:
                v = True
            elif v.lower() in ["no", "false"]:
                v = False
            elif v.replace(".", "").isnumeric():
                v = float(v)
            elif v.isnumeric():
                v = int(v)
            d[k.strip().replace("_", "-")] = v

    return d
