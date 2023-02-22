#!/usr/bin/env bash
#
gmx trjcat -cat -o cat.trr  -f \
  hat_tf_000/1_equilibration/equil.trr\
  hat_tf_000/2_equilibration/equil.trr\
  hat_tf_000/5_relaxation/relax.trr\
  hat_tf_000/6_equilibration/equil.trr\
  hat_tf_000/9_relaxation/relax.trr\
  hat_tf_000/10_equilibration/equil.trr\
  hat_tf_000/13_relaxation/relax.trr\
  hat_tf_000/14_equilibration/equil.trr\
  hat_tf_000/17_relaxation/relax.trr\
  hat_tf_000/18_equilibration/equil.trr\
  hat_tf_000/21_relaxation/relax.trr\
  hat_tf_000/22_equilibration/equil.trr\
  hat_tf_000/25_relaxation/relax.trr\
  hat_tf_000/26_equilibration/equil.trr\
  hat_tf_000/29_relaxation/relax.trr\
  hat_tf_000/30_equilibration/equil.trr\
  hat_tf_000/33_relaxation/relax.trr\
  hat_tf_000/34_equilibration/equil.trr\
  hat_tf_000/37_relaxation/relax.trr\
  hat_tf_000/38_equilibration/equil.trr\
  hat_tf_000/41_relaxation/relax.trr


echo "1 0" | gmx trjconv -dump 0 -f cat.trr -s ./hat_tf_000/1_equilibration/equil.tpr -center -pbc mol -o cat-center.gro
echo "1 0" | gmx trjconv -f cat.trr -s ./hat_tf_000/1_equilibration/equil.tpr -center -pbc mol -o cat-center.xtc



popd

