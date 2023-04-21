import MDAnalysis as mda
from pathlib import Path

tpr = Path("/hits/fast/mbm/hartmaec/kimmdys/kimmdy/example/example_triala/test_out_003/1_equilibration/equil.tpr")
trr = Path("/hits/fast/mbm/hartmaec/kimmdys/kimmdy/example/example_triala/test_out_003/1_equilibration/equil.trr")

u = mda.Universe(str(tpr), str(trr), topology_format="tpr", format="trr")

print(u)

for ts in u.trajectory:
    print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, u.trajectory.time))

