import MDAnalysis as mda
from pathlib import Path
import os
import matplotlib.pyplot as plt
from kimmdy.tasks import TaskFiles
from kimmdy.reactions.homolysis import Homolysis
from kimmdy.config import Config
from kimmdy.runmanager import RunManager
from kimmdy.parsing import read_plumed, read_distances_dat, read_plumed_distances, read_rtp


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
    cwd = Path("tests/test_files/test_KMC")
    os.chdir(cwd)
    yml = "kimmdy.yml"
    print(yml)


    r = Homolysis(name='homolysis',runmng=RunManager(Config(yml)))
    print(r)

    return

def test_parsers():
    cwd = Path("tests/test_files/test_KMC")
    os.chdir(cwd)
    plumed = read_plumed(Path('plumed.dat'))
    distances = read_distances_dat(Path('distances.dat'))
    distances_v2 = read_plumed_distances(Path('plumed.dat'),Path('distances.dat'))
    edissoc = None
    ffbonded = read_rtp(Path('ffbonded.itp'))
    return


def plot_time_evolution():
    cwd = Path("tests/test_files/test_KMC")
    os.chdir(cwd)
    distances = read_distances_dat(Path('distances.dat'))
    for key in distances.keys():
        if key == 'time':
            continue
        plt.plot(distances['time'],distances[key])
        plt.show()

#test_KMC()
test_parsers()
#plot_time_evolution()



