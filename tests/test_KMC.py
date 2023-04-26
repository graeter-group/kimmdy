# %%
import MDAnalysis as mda
from pathlib import Path
from ast import literal_eval
import os
import dill
import dill.detect
import numpy as np
import matplotlib.pyplot as plt
from kimmdy.tasks import TaskFiles
from kimmdy.reaction import ReactionResult
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.config import Config
from kimmdy.runmanager import RunManager
from kimmdy.parsing import read_plumed, read_distances_dat, read_plumed_distances, read_topol, read_rtp, read_edissoc
from kimmdy.utils import calc_transition_rate

#%%
print(os.getcwd())
cwd = Path("tests/test_files/test_KMC")
os.chdir(cwd)
print(os.getcwd())
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
    from numpy.random import default_rng
    rng = default_rng(1)

    # cwd = Path("tests/test_files/test_KMC")
    # os.chdir(cwd)
    rmgr = RunManager(Config(Path("kimmdy.yml")))
    files = TaskFiles(rmgr)
    files.input['top'] = Path('topol.top')
    for f_path in ['plumed.dat','distances.dat','edissoc.dat','ffbonded.itp']:
        files.input[f_path] = Path(f_path)
                      
    r = Homolysis(name='homolysis',runmng=rmgr)
    RR = r.get_reaction_result(files)
    print(RR)
    for RO in RR:
        print(RO)

    #RR.to_csv(Path("RR_test.csv"))
    RR.to_dill(Path("RR_test.dill"))
    # RR_parsed = read_ReactionResult(Path("RR_test.csv"))
    # print(RR_parsed)
    # if RR == RR_parsed:
    #     print("IDENTITY!")
    # else:
    #     print("FAILURE")

    # for RO in RR:
    #     r_ts = RO.r_ts
    #     ts = RO.ts
    #     plt.plot(ts,r_ts)
    #     plt.plot([min(ts),max(ts)],np.full(2,RO.rate),color='k',linewidth=2,label='dist avg')
    #     plt.scatter(ts,r_ts)
    #     plt.yscale('log')
    #     plt.legend()
    #     plt.show()

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

def test_RR_parsing():
    RR_blank = ReactionResult()
    # RR_ref = read_ReactionResult(RR_blank, Path("RR_ref.csv"))

    # write_ReactionResult(RR_ref, Path("RR_test.csv"))
    # RR_blank = ReactionResult()
    # RR_test = read_ReactionResult(RR_blank, Path("RR_test.csv"))
    
    # print(RR_test)
    # if RR_ref == RR_test:
    #     print("IDENTITY!")
    # else:
    #     print("FAILURE")

def test_pandas_lists():
    import pandas as pd
    # df = pd.DataFrame(
    #     {'trial_num': [1, 2, 3, 1, 2, 3],
    #     'subject': [1, 1, 1, 2, 2, 2],
    #     'samples': [list(np.random.randn(3).round(2)) for i in range(6)]
    #     }
    # )   
    # df.to_csv(Path('pandas_test.csv'))
    df_new = pd.read_csv(Path('RR_test.csv'))
    print(df_new)
    return

def test_RR_dill():
    #resource.setrlimit(resource.RLIMIT_STACK, [0x100*500, resource.RLIM_INFINITY])
    print('testing...')

    with open(Path('RR_test.dill'),"rb") as f:
        RR = dill.load(f)
    print(RR)
    RR.to_csv(Path('RR_test2.csv'))
    return

# %%
test_homolysis()
#test_KMC()
#test_parsers()
#plot_time_evolution()
#test_RR_parsing()
#test_pandas_lists()
test_RR_dill()
print('x')
# %%
