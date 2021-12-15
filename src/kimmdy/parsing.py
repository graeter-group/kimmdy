from pathlib import Path
from collections.abc import Iterable
from typing import Generator

from kimmdy.reaction import Topology


def is_comment(l: str):
    return len(l) == 0 or l[0] in ["#", "\n", ";"]


def get_sections(
    seq: Iterable[str], section_marker: str
) -> Generator[list[str], None, None]:
    data = []
    for line in seq:
        if is_comment(line):
            continue
        if line.startswith(section_marker):
            if data:
                yield data
                data = []
        data.append(line.strip("\n"))
    if data:
        yield data


def read_topol(path: Path) -> Topology:
    with open(path, "r") as f:
        sections = get_sections(f, "[")
        return {
            title.strip("[] \n"): [c.split() for c in content]
            for title, *content in sections
        }


def write_topol(d: Topology, outfile: Path):
    with open(outfile, "w") as f:
        for title, content in d.items():
            s = f"[ {title} ]\n"
            f.write(s)
            s = "\n".join([" ".join(c) for c in content]) + "\n\n"
            f.write(s)


def read_plumed(path: Path):
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
            if "PRINT" in l[:5]:
                l = l.split()
                d = {
                    "PRINT": l[0]
                }
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


def write_plumed(d, path: Path):
    with open(path, "w") as f:
        for l in d["distances"]:
            f.write(f"{l['id']}: {l['keyword']} ATOMS={','.join(l['atoms'])}\n")
        for l in d["prints"]:
            f.write(f"{l['PRINT']} ARG={','.join(l['ARG'])} STRIDE={str(l['STRIDE'])} FILE={str(l['FILE'])}\n")

