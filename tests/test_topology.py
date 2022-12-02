#%%
from copy import deepcopy
from itertools import takewhile, permutations
from pathlib import Path
import os
from xml.etree.ElementTree import Element
from kimmdy.parsing import is_not_comment, read_topol, read_xml_ff
from kimmdy.changemanager import LocalGraph
from kimmdy.topology import FF, Topology, Atom, Bond, Dihedral, Pair, get_by_permutations, get_element_id

#%%
def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_topology"
    except NameError:
        test_dir = Path("./tests/test_files/test_topology")
    os.chdir(test_dir)


set_dir()

#%%
hexala_top = read_topol(Path('hexala.top'))
ffdir = Path("../assets/amber99sb-star-ildnp.ff")
ffpatch = Path('amber99sb_patches.xml')
top = Topology(hexala_top, ffdir, ffpatch)
oldtop = deepcopy(top)

# 9 = CA
# 10 = HA
# in residue 2 Ala
# 9        CT       2        ALA      CA       9        0.0337   12.01   
# 10       H1       2        ALA      HA       10       0.0823   1.008   

#%%
top.break_bond(('9', '10'))

#%%
for p in top.ffpatches.atompatches:
    print(get_element_id(p))

#%%
filter(lambda x: x.funct == '9', top.dihedrals.values()).__next__()


#%%

