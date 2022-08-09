#!/usr/bin/env bash
echo "1 0" | gmx trjconv -f 1_equilibration/equil.trr -s 1_equilibration/equil.tpr -o 1_equil.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 5_production/prod.trr -s 5_production/prod.tpr -o 5_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 9_production/prod.trr -s 9_production/prod.tpr -o 9_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 13_production/prod.trr -s 13_production/prod.tpr -o 13_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 17_production/prod.trr -s 17_production/prod.tpr -o 17_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 21_production/prod.trr -s 21_production/prod.tpr -o 21_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 25_production/prod.trr -s 25_production/prod.tpr -o 25_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 1_equilibration/equil.trr -s 1_equilibration/equil.tpr -o 1_equil.trr -center -pbc mol
gmx trjcat -f 1_equil.trr 5_prod.trr 9_prod.trr 13_prod.trr 17_prod.trr 21_prod.trr 25_prod.trr -cat


