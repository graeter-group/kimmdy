#!/usr/bin/env bash
printf "1\n1\n " | gmx pdb2gmx -f pep.pdb -o pep.gro -p pep_out.top -ignh
gmx editconf -f pep.gro -o pep_box.gro -c -d 3.5 -bt dodecahedron
gmx solvate -cp pep_box.gro -p pep_out.top -o pep_solv.gro

gmx grompp -f ions.mdp -c pep_solv.gro -p pep_out.top -o pep_out_genion.tpr
echo "SOL" | gmx genion -s pep_out_genion.tpr -p pep_out.top -o pep_out_ion.gro -conc 0.15 -neutral

gmx grompp -f minim.mdp -c pep_out_ion.gro -p pep_out.top -o pep_out_min.tpr
gmx mdrun -deffnm pep_out_min -v
wait
gmx grompp -f nvt.mdp -c pep_out_min.gro -p pep_out.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
wait
gmx grompp -f npt.mdp -c nvt.gro -p pep_out.top -o npt.tpr
gmx mdrun -v -deffnm npt
wait
echo "q" | gmx make_ndx -f npt.gro
