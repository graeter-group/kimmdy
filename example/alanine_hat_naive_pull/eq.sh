#!/bin/env bash

# %%
cp topol_init.top topol.top
gmx editconf -f ./init.gro -o ./init_box.gro -c -d 4.0 -bt triclinic
gmx solvate -cp ./init_box.gro  -o ./init_solv.gro -p ./topol.top
gmx grompp -f ./ions.mdp -c ./init_solv.gro -p ./topol.top -o ./ions.tpr

echo '13' | gmx genion -s ./ions.tpr -o ./init_solv_ions.gro -p ./topol.top -pname NA -nname CL -neutral -conc 0.15

gmx grompp -f ./minim.mdp -c ./init_solv_ions.gro -p ./topol.top -o ./em.tpr
gmx mdrun -v -deffnm em

gmx grompp -f ./nvt.mdp -c ./em.gro -p ./topol.top -o ./nvt.tpr
gmx mdrun -v -deffnm nvt

gmx grompp -f ./npt.mdp -c ./nvt.gro -p ./topol.top -o ./npt.tpr
gmx mdrun -v -deffnm npt

echo -e 'q\n' | gmx make_ndx -f ./npt.gro

cat <<EOF >> index.ndx
[ ace ]
1 2 3 4 5 6
[ nme ]
16 17 18 19 20 21
EOF


rm \#*\#

