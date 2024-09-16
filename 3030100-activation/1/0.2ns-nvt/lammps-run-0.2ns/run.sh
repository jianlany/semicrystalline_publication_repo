#!/bin/bash

f=step_109991
ts="300 310 320 330 340 350 360"
p=1
for t in ${ts}; do
  cd ${f}-${t}K
  cp ../in.lammps .
  sed -i "s/PRESSURE_VALUE/${p}/g" in.lammps
  sed -i "s/TEMPERATURE_VALUE/${t}/g" in.lammps
  sed -i "s/FNAME_VALUE/${f}/g" in.lammps

  if [ -n "`which squeue`" ]; then
    echo "running full test through slurm"
    sbatch ../slurm.submit
  else 
    echo "running short test"
    mpiexec -n 4 lmp_sparky -in in.lammps
  fi
  cd ..
done
