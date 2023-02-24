#!/usr/bin/env bash
echo "1 0" | gmx trjconv -f 2_production/prod.trr -s 2_production/prod.tpr -o 2_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 7_production/prod.trr -s 7_production/prod.tpr -o 7_prod.trr -center -pbc mol
echo "1 0" | gmx trjconv -f 12_production/prod.trr -s 12_production/prod.tpr -o 12_prod.trr -center -pbc mol
gmx trjcat -f 2_prod.trr 7_prod.trr 12_prod.trr -cat


