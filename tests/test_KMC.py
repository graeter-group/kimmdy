# %%
import MDAnalysis as mda
from pathlib import Path
import os
import matplotlib.pyplot as plt
from kimmdy.tasks import TaskFiles
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.config import Config
from kimmdy.runmanager import RunManager
from kimmdy.parsing import read_plumed, read_distances_dat, read_plumed_distances, read_topol, read_rtp, read_edissoc
from kimmdy.utils import calc_transition_rate
# %%

def test_KMC():
    tpr = Path("../example/example_triala/test_out_003/1_equilibration/equil.tpr")
    trr = Path("../example/example_triala/test_out_003/1_equilibration/equil.trr")

    u = mda.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

    print(u)

    for ts in u.trajectory:
        print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, u.trajectory.time))
    return

## this should be somewhere else
def test_homolysis():
    # %%
    from numpy.random import default_rng
    rng = default_rng(1)

    cwd = Path("tests/test_files/test_KMC")
    os.chdir(cwd)
    # %%
    rmgr = RunManager(Config(Path("kimmdy.yml")))
    files = TaskFiles(rmgr)
    files.input['top'] = Path('topol.top')
    for f_path in ['plumed.dat','distances.dat','edissoc.dat','ffbonded.itp']:
        files.input[f_path] = Path(f_path)
                      
    r = Homolysis(name='homolysis',runmng=rmgr)
    # %%
    RR = [r.get_reaction_result(files)]
    print(RR)

    # %%
    return

def test_parsers():
    cwd = Path("tests/test_files/test_KMC")
    os.chdir(cwd)
    top = read_topol(Path('topol.top'))
    plumed = read_plumed(Path('plumed.dat'))
    distances = read_distances_dat(Path('distances.dat'))
    edissoc_lookup = read_edissoc(Path('edissoc.dat'))
    ffbonded_dict = read_rtp(Path('ffbonded.itp'))
    ffbonded_lookup = {tuple(l[:2]):{'b0':float(l[3]),'kb':float(l[4])} for l in ffbonded_dict['bondtypes']['other']}
    return


def plot_time_evolution():
    cwd = Path("tests/test_files/test_KMC")
    os.chdir(cwd)
    distances = read_distances_dat(Path('distances.dat'))
    for i,key in enumerate(distances.keys()):
        #k, F = calc_transition_rate(r_av, r_0, E_dis, k_f)
        if key == 'time':
            continue
        plt.plot(distances['time'],distances[key],alpha=0.7)
        if i%20 == 19:
            plt.show()
            plt.figure()
        
# %%
test_homolysis()
#test_KMC()
#test_parsers()
#plot_time_evolution()

# %%
