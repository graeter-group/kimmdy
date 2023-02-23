#!/usr/bin/env bash

dir='002'

gmx trjcat -cat -o cat.trr  -f \
  hat_tf_${dir}/1_equilibration/equil.trr\
  hat_tf_${dir}/2_equilibration/equil.trr\
  hat_tf_${dir}/5_relaxation/relax.trr\
  hat_tf_${dir}/6_equilibration/equil.trr\
  hat_tf_${dir}/9_relaxation/relax.trr\
  hat_tf_${dir}/10_equilibration/equil.trr\
  hat_tf_${dir}/13_relaxation/relax.trr\
  hat_tf_${dir}/14_equilibration/equil.trr\
  hat_tf_${dir}/17_relaxation/relax.trr\
  hat_tf_${dir}/18_equilibration/equil.trr\
  hat_tf_${dir}/21_relaxation/relax.trr\
  hat_tf_${dir}/22_equilibration/equil.trr\
  hat_tf_${dir}/25_relaxation/relax.trr\
  hat_tf_${dir}/26_equilibration/equil.trr\
  hat_tf_${dir}/29_relaxation/relax.trr\
  hat_tf_${dir}/30_equilibration/equil.trr\
  hat_tf_${dir}/33_relaxation/relax.trr\
  hat_tf_${dir}/34_equilibration/equil.trr\
  hat_tf_${dir}/37_relaxation/relax.trr\
  hat_tf_${dir}/38_equilibration/equil.trr\
  hat_tf_${dir}/41_relaxation/relax.trr


echo "1 1" | gmx trjconv -dump 0 -f cat.trr -s ./hat_tf_000/1_equilibration/equil.tpr -center -pbc mol -o cat-center.gro
echo "1 1" | gmx trjconv -f cat.trr -s ./hat_tf_000/1_equilibration/equil.tpr -center -pbc mol -o cat-center.xtc


