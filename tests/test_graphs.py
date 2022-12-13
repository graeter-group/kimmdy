#%%
from copy import deepcopy
from itertools import takewhile
from pathlib import Path
import re
import os
from typing import Optional
from xml.etree.ElementTree import Element
import pytest
from kimmdy.parsing import read_topol, read_xml_ff
from kimmdy.changemanager import LocalGraph
from kimmdy.topology import Topology
import difflib as dl

#%%
def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_topology"
    except NameError:
        test_dir = Path("./tests/test_files/test_topology")
    os.chdir(test_dir)


set_dir()

#%%
col_top_dict = read_topol(
    Path(
        "/hits/fast/mbm/buhrjk/phd/col-hydrolysis/col-fibril-crosslinks/run1/topol.top"
    )
)
col_top = Topology(col_top_dict, ffdir)

col_top.bonds[:10]


#%%
hexala_graph = LocalGraph(hexala_top, ("1", "2"), ffdir)


#%%
xml_ff = read_xml_ff(Path("amber99sb_trunc.xml"))
xml_ff[0]

atomtypes = [type.attrib for type in xml_ff.findall("AtomTypes/")]

#%%
for x in xml_ff.iter():
    print(x)


#%%
patch = read_xml_ff(Path("amber99sb_patches.xml"))

#%%
for x in patch.iter():
    print(x)

#%%
p = patch.findall('HarmonicBondForce/Bond[@class1="*_R"]')
s = set(p)
p2 = patch.findall('HarmonicBondForce/Bond[@class1="*_R"]')
s2 = set(p2)
list(s.union(s2))

#%%
# need to find the most specific radical patch,
# a bit like in CSS
# e.g.
patch_names = [
    e.get("class1") for e in patch.findall("HarmonicBondForce/Bond[@class1]")
]
# ['*_R', '*_R', 'O+_R', 'O+', 'S+_R', 'S+', 'N', 'C']
# should match '*_R'


def match_attr(patches: list[Element], attr: str, m: str) -> Optional[Element]:
    matches = []
    for p in patches:
        if value := p.get(attr):
            if value == m:
                return p
            pattern = value.replace("*", r".*").replace("+", r"\+")
            if re.match(pattern, m):
                matches.append(p)
    if matches:
        matches.sort(key=lambda x: x.get(attr))
        return matches[0]
    else:
        return None


ms = match_attr(patch.findall("*/Bond"), "class1", "C")

ms.items()

#%%
match_attr(patch.findall("*/Bond"), "class1", "S+_R")

#%%
match_attr(patch.findall("*/Bond"), "class1", "N")
