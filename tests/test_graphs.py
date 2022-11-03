#%%
from pathlib import Path
import os
from kimmdy.parsing import read_topol, read_xml_ff
from kimmdy.changemanager import LocalGraph

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
hexala_graph = LocalGraph(hexala_top, ('1', '2'), ffdir)
# hexala_graph

#%%
patch = read_xml_ff(Path('amber99sb_trunc.xml'))
patch[0]

atomtypes = [type.attrib for  type in patch.findall('AtomTypes/')]




