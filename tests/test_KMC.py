# %%
import MDAnalysis as mda
from pathlib import Path
from ast import literal_eval
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
from matplotlib.pyplot import cm
from kimmdy.tasks import TaskFiles
from kimmdy.reaction import ReactionResult, ReactionOutcome, Conversion, ConversionRecipe, ConversionType
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.config import Config
from kimmdy.runmanager import RunManager
from kimmdy.parsing import read_plumed, read_distances_dat, read_plumed_distances, read_topol, read_rtp, read_edissoc
from kimmdy.utils import calc_transition_rate
from kimmdy.KMC import rfKMC, running_mean_uniform_filter1d
import seaborn as sns
from itertools import cycle

#%%
print(os.getcwd())
#cwd = Path("tests/test_files/test_KMC")
#cwd = Path("example/example_short/")
#cwd = Path("/hits/fast/mbm/hartmaec/workdir/KIMMDY_trjs/KMC_test/run6_50ns")
cwd = Path("/hits/fast/mbm/hartmaec/workdir/KIMMDY_trjs/KMC_test/run7_100ns")

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
    # for RO in RR:
    #     print(RO)

    RR.to_dill(Path("RR_ref.dill"))
    RR_new = ReactionResult.from_dill(Path("RR_ref.dill"))
    RR_new.to_csv(Path("RR_ref.csv"))
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
    print(os.getcwd())
    distances = read_distances_dat(Path('distances_5000fs.dat'))
    for i,key in enumerate(distances.keys()):
        #k, F = calc_transition_rate(r_av, r_0, E_dis, k_f)
        if key == 'time':
            continue
        plt.plot(distances['time'],distances[key],alpha=0.7,linewidth=1)
        plt.scatter(distances['time'],distances[key],s=2,alpha=0.7)
        if i==1:
            plt.show()
            plt.figure()
            break

def plot_running_avg_distances():
    print('starting...')
    distances = read_distances_dat(Path('distances.dat'))
    #distances = read_distances_dat(Path('distances_10000.dat'))
    print(len(distances['time']))

    plt.scatter(distances['time'], distances['d100'],s=2,alpha=0.7,label='d100',c='cornflowerblue')
    plt.scatter(distances['time'], distances['d200'],s=2,alpha=0.7,label='d200',c='forestgreen')
    plt.scatter(distances['time'], distances['d300'],s=2,alpha=0.7,label='d300',c='tomato')

    plt.plot([min(distances['time']),max(distances['time'])], [np.average(distances['d100']),np.average(distances['d100'])], c= 'darkblue',label='d100 average')
    plt.plot([min(distances['time']),max(distances['time'])], [np.average(distances['d200']),np.average(distances['d200'])], c= 'darkgreen',label='d200 average')
    plt.plot([min(distances['time']),max(distances['time'])], [np.average(distances['d300']),np.average(distances['d300'])], c= 'darkred',label='d300 average')

    #plt.legend()
    plt.xlabel('t [ns]')
    plt.ylabel('x [nm]',)
    plt.title("Distance of breakpairs in a tripelhelix")
    plt.tight_layout()
    plt.savefig(Path("/hits/fast/mbm/hartmaec/Labbook/figures_230515") / "running_avg_distances.png",dpi=300)
    # for key in distances.keys():
    #     if key == 'time':
    #         continue
    #     plt.plot([min(distances['time']),max(distances['time'])], [min(distances[key]),max(distances[key])], c= 'k')
    #     plt.scatter(distances['time'], distances[key],alpha=0.7)
        # for i in np.power(2,range(1,25)):
        #     print(i)
        #     if i < len(distances[key]):
        #         block_dists = running_mean_uniform_filter1d(np.asarray(distances[key]),i)
        #         print(block_dists)
        #         plt.plot(distances['time'], block_dists,label=f"block avg N: {i}",alpha=0.7)
                    
        #     else:
        #         plt.legend()
        #         plt.show()
        #         break

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
    print('testing...')

    RR = ReactionResult([ReactionOutcome(recipe=[Conversion(ConversionType.BREAK, ('0','1'))], rate=0.5, r_ts=[0.5], ts=[1] )])
    RR.to_dill(Path('RR_test.dill'))
    RR_new = ReactionResult.from_dill(Path('RR_test.dill'))
    print(RR_new)
    RR_new.to_csv(Path('RR_test2.csv'))
    return

def test_decision_strategy():
    result = ReactionResult.from_dill(Path('RR_ref.dill'))
    k_dist_avg = [outcome.rate for outcome in result.outcomes]

    _, k_integrate = rfKMC(result)

    # plt.scatter(range(len(k_dist_avg)),k_dist_avg,label='dist_avg',alpha=0.5)
    # plt.scatter(range(len(k_integrate)),k_integrate,label='k_integrate',alpha=0.5)
    # plt.yscale('log')
    # plt.show()
    # plt.figure()

    k_dist_avg = np.cumsum(k_dist_avg)/np.sum(k_dist_avg)
    k_integrate = np.cumsum(k_integrate)/np.sum(k_integrate)

    print(np.sum(k_integrate),np.sum(k_dist_avg))
    print(np.max(k_dist_avg),np.max(k_integrate))

    print(len(result.outcomes),len(k_dist_avg),len(k_integrate))
    for i in range(len(k_dist_avg)-1,0,-1): 
        plt.bar('k_dist',k_dist_avg[i])
        plt.bar('k_integrate',k_integrate[i])
    #
    plt.legend()
    plt.show()

def calc_rates():
    for i in [1,"state_avg_dist_ca"]:#[1]:#[100,1000,10000,100000]:
        rmgr = RunManager(Config(Path("kimmdy.yml")))
        files = TaskFiles(rmgr)
        files.input['top'] = Path('topol.top')
        files.input['distances.dat'] = Path(f"distances_10000fs_{i}.dat")
        if i == 1:
            files.input['distances.dat'] = Path(f"distances_10000fs.dat")
        for f_path in ['plumed.dat','edissoc.dat','ffbonded.itp']:
            files.input[f_path] = Path(f_path)
                        
        r = Homolysis(name='homolysis',runmng=rmgr)
        RR = r.get_reaction_result(files)       
        RR.to_dill(Path(f"RR_10000fs_{i}.dill"))

def write_running_avg_dist():
    def write_dists(block_dists_dict,N):
        with open(f"distances_50fs_{N}.dat", "w") as f:
            f.write('#! FIELDS ' + ' '.join(block_dists_dict.keys()))
            for v in zip(*block_dists_dict.values()):
                f.write("\n" + " ".join("{:9.7f}".format(x) for x in v))
            f.write('\n')
        return

    distances = read_distances_dat(Path('distances_50fs.dat'))
    print(len(distances['time']))
    for i in [100,1000,10000,100000]:
        block_dists_dict = {}
        block_dists_dict['time'] = np.asarray(distances['time'])
        for key in distances.keys():
            if key == 'time':
                continue
            block_dists = running_mean_uniform_filter1d(np.asarray(distances[key]),i)
            block_dists_dict[key] = block_dists
        write_dists(block_dists_dict,i)

def get_state_idxs():
    with open("run7_states_PCA_dist_ca.txt",'r') as f:
        foo = f.readlines()
    states = np.asarray([int(x) for x in foo[0].strip().split(sep=',')])
    stateidxs = []
    weights = []
    for i in range(-1,max(states)+1):
        idxs = np.argwhere(states == i)
        stateidxs.append(idxs)
        weights.append(len(states[stateidxs[-1]]))
    print(len(states),weights)
    return stateidxs, weights

def write_state_avg_dist():
    def write_dists(dists_dict,):
        with open(f"distances_10000fs_state_avg_dist_ca.dat", "w") as f:
            f.write('#! FIELDS ' + ' '.join(dists_dict.keys()))
            for v in zip(*dists_dict.values()):
                f.write("\n" + " ".join("{:9.7f}".format(x) for x in v))
            f.write('\n')
        return
    distances = read_distances_dat(Path('distances_10000fs.dat'))
    distances = {k:np.asarray(v) for k,v in distances.items()}
    stateidxs, weights = get_state_idxs()
    for i,key in enumerate(distances.keys()):
        if key == 'time':
            continue
        for i,idxs in enumerate(stateidxs):
            if i == 0:
                distances[key][idxs] = 0
            else:
                distances[key][idxs] = np.average(distances[key][idxs])
    write_dists(distances)
    return
    
def plot_running_avg_rates():
    for i in [1]:
        result = ReactionResult.from_dill(Path(f"RR_10000fs_{i}.dill"))
        print(result.outcomes[0].r_ts[0],result.outcomes[0].r_ts[-1])
        _, k_integrate = rfKMC(result)
        k_integrate = np.cumsum(k_integrate)/np.sum(k_integrate)
        colours = cycle(sns.color_palette("pastel", 20))
        for j in range(len(k_integrate)-1,0,-1): 
            color = next(colours)
            
            plt.bar(f"no_avg",k_integrate[j],color=color)
        
    result = ReactionResult.from_dill(Path(f"RR_10000fs_state_avg_dist_ca.dill"))
    k_dist_avg = [outcome.rate for outcome in result.outcomes]
    k_dist_avg = np.cumsum(k_dist_avg)/np.sum(k_dist_avg)
    print(result.outcomes[0].r_ts)
    _, k_integrate = rfKMC(result)
    k_integrate = np.cumsum(k_integrate)/np.sum(k_integrate)
    colours = cycle(sns.color_palette("pastel", 20))
    for j in range(len(k_integrate)-1,0,-1):  
        color = next(colours)
        plt.bar('state_avg_resdist',k_integrate[j],color=color)

    # result = ReactionResult.from_dill(Path(f"RR_10000fs_state_avg_rmsd.dill"))
    # k_dist_avg = [outcome.rate for outcome in result.outcomes]
    # k_dist_avg = np.cumsum(k_dist_avg)/np.sum(k_dist_avg)
    # print(result.outcomes[0].r_ts)
    # _, k_integrate = rfKMC(result)
    # k_integrate = np.cumsum(k_integrate)/np.sum(k_integrate)
    # colours = cycle(sns.color_palette("pastel", 20))
    # for j in range(len(k_integrate)-1,0,-1):  
    #     color = next(colours)
    #     plt.bar('state_avg_rmsd',k_integrate[j],color=color)

    colours = cycle(sns.color_palette("pastel", 20))
    for j in range(len(k_integrate)-1,0,-1):  
        color = next(colours)
        plt.bar('overall',k_dist_avg[j],color=color) 

    plt.ylabel("probability")
    plt.xlabel("rate averaging scheme")
    plt.title("homolysis rates")
    plt.xticks(rotation = 45)
    plt.tight_layout()
    plt.savefig(Path("/hits/fast/mbm/hartmaec/Labbook/figures_230515") / "barplot_rates_run7_10000fs_states.png",dpi=300)
    plt.show()



def check_result():
    result = ReactionResult.from_dill(Path(f"RR_ref.dill"))
    print(result.outcomes[0].r_ts)
    result.to_csv(Path(f"RR_ref_1.csv"))


# def cluster_DBscan():
#     return

# %%
#test_homolysis()
#test_KMC()
#test_parsers()
#plot_time_evolution()
#plot_running_avg_distances()
#test_RR_parsing()
#test_pandas_lists()
#test_RR_dill()
#test_decision_strategy()
#write_running_avg_dist()


#write_state_avg_dist()
calc_rates()
plot_running_avg_rates()
#get_state_idxs()
#check_result()

# %%
