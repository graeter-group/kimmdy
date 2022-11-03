#%%
from pathlib import Path
import os
from kimmdy.parsing import read_topol, read_xml_ff
from kimmdy.changemanager import LocalGraph
from kimmdy.topology import Topology

#%%
def set_dir():
    try:
        test_dir = Path(__file__).parent / "test_files/test_graphs"
    except NameError:
        test_dir = Path("./tests/test_files/test_graphs")
    os.chdir(test_dir)


set_dir()

#%%
hexala_top = read_topol(Path('hexala.top'))
ffdir = Path("../assets/amber99sb-star-ildnp.ff")

#%%
top = Topology(hexala_top, ffdir)
top.improper_dihedrals

#%%
col_top_dict = read_topol(Path('/hits/fast/mbm/buhrjk/phd/col-hydrolysis/col-fibril-crosslinks/run1/topol.top'))
col_top = Topology(col_top_dict, ffdir)





#%%

hexala_graph = LocalGraph(hexala_top, ('1', '2'), ffdir)


#%%
patch = read_xml_ff(Path('amber99sb_trunc.xml'))
patch[0]

atomtypes = [type.attrib for  type in patch.findall('AtomTypes/')]




