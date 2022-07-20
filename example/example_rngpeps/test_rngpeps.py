import shutil
import os
from pathlib import Path
import subprocess as sp
import numpy as np
rng = np.random.default_rng()


basedir = Path("/hits/fast/mbm/hartmaec/kimmdys/kimmdy_topology/example/example_rngpeps")
del_old = True

for i in range(10):
    currdir = basedir / str(i)
    if del_old:
        if currdir.exists():
            print(f"Deleted dir {currdir}")
            shutil.rmtree(currdir)
    shutil.copytree(basedir/"files_md",basedir/str(i),dirs_exist_ok=False,symlinks=True)
    AAstring = rng.choice(['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','P','Y'],5) #extra P for Hydroxyproline, Y for DOPA
    print(AAstring,''.join(AAstring))
    sp.run(f"pymol -c -d 'fab B{''.join(AAstring)}Z, pep{i}, ss=2;set pdb_use_ter_records, 0;save pep{i}.pdb, pep{i}'",shell=True,cwd=currdir)

    P_pos = [x for x in range(len(AAstring)) if AAstring[x] == 'P']
    Y_pos = [x for x in range(len(AAstring)) if AAstring[x] == 'Y']
    print(P_pos,Y_pos)
    for pos in P_pos:
        if rng.choice([0,1],1):
            print(f"Modify {pos} P")
            sp.run(f"sed 's/pep_i/pep{i}/g' ProtoHyp_blank.pml | sed 's/pos/{pos+2}/g' > {i}/ProtoHyp_{i}_{pos+2}.pml",shell=True,cwd=basedir)
            sp.run(f"pymol -c ProtoHyp_{i}_{pos+2}.pml",shell=True,cwd=currdir)

    for pos in Y_pos:
        if rng.choice([0,1],1):
            print(f"Modify {pos} Y")
            sp.run(f"sed 's/pep_i/pep{i}/g' TyrtoDop_blank.pml | sed 's/pos/{pos+2}/g' > {i}/TyrtoDop_{i}_{pos+2}.pml",shell=True,cwd=basedir)
            sp.run(f"pymol -c TyrtoDop_{i}_{pos+2}.pml",shell=True,cwd=currdir)

    shutil.copy2(currdir/f"pep{i}.pdb",currdir/"pep.pdb")
    sp.run("bash quick_gmx.sh",shell=True,cwd=currdir)

    