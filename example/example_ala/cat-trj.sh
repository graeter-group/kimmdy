#!/usr/bin/env bash

for path in\
  hat_tf_000/1_equilibration/equil\
  hat_tf_000/2_equilibration/equil\
  hat_tf_000/5_relaxation/relax\
  hat_tf_000/6_equilibration/equil\
  hat_tf_000/9_relaxation/relax; do
echo "1 0" | gmx trjconv -f ${path}.trr -s ${path}.tpr -o ${path}-center.trr -center -pbc mol
done


gmx trjcat -cat -o cat-center.xtc -f \
  hat_tf_000/1_equilibration/equil-center.trr\
  hat_tf_000/2_equilibration/equil-center.trr\
  hat_tf_000/5_relaxation/relax-center.trr\
  hat_tf_000/6_equilibration/equil-center.trr\
  hat_tf_000/9_relaxation/relax-center.trr

echo "1 0" | gmx trjconv -dump 0 -f cat.trr -s ./hat_tf_000/1_equilibration/equil.tpr -center -pbc mol -o cat-center.gro

