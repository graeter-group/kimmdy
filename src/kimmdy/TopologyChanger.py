#%%
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from collections.abc import Iterable
from typing import Generator
from copy import deepcopy
import itertools
from operator import itemgetter

#%%
# from parsing
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

class ConversionType(Enum):
    BREAK = auto()
    MOVE = auto()

@dataclass
class ConversionRecipe:
    """A ConversionRecipe.
    encompasses a single transformation, e.g. moving one
    atom or braking one bond.
    Parameters
    ----------
    type : list[ConversionType.BREAK or .MOVE]
    atom_idx : list[(from, to)]
    """

    type: list[ConversionType] = field(default_factory=list)
    atom_idx: list[tuple[int, int]] = field(default_factory=list)








    
    # atom_terms = fromGraph.parameterize_around_radical()

    # print(f"\n---- Changing the Topology ----")
    # print(topoldict['dihedrals'])
    # #topoldict = topol_add_terms(topoldict,atom_terms)
    # topoldict = topol_remove_terms(topoldict,atom_terms)
    # print('\n',topoldict.keys())
    # write_topol(topoldict,outpath)
    # print('\n',topoldict['dihedrals'])

    # print('\n\n--- Wrapping Up ---')
    # print('Atomlist: ',fromGraph.AtomList)
    # print(fromGraph.get_AtomList_property('idx'))
    # print(fromGraph.get_AtomList_property('bound_to'))
    # print(fromGraph.AtomList[0].idx)
    # print('BondList: ',fromGraph.BondList)
    # print('PairList: ',fromGraph.PairList)
    # print('AngleList: ',fromGraph.AngleList)
    # print('DihedralList: ',fromGraph.DihedralList)
    # print('\n')
    # print(fromGraph.get_terms_with_atom(str(22)))
    # print(fromGraph.find_missing_terms())
    # print(atom_terms)



# %%
