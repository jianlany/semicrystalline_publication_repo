#!/bin/bash

## slurm input settings - ignored if not submitted through slurm.
#SBATCH --job-name semicrystalline
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --time 500:00:00
#SBATCH --error=job.%J.out --output=job.%J.out


if [ $(hostname) == "sparky" ]; then
  echo "running full test through torque"
  workdir=$(pwd)
  echo " 
  #!/bin/bash
  #PBS -l walltime=200:00:00
  #PBS -l nodes=1:ppn=16
  #PBS -N semicrystalline
  #PBS -j oe
  #PBS -o o_semicrystalline
  cd $workdir
  ./make_system.py build input-full.sc
  " | qsub
elif [ -n "`which squeue`" ]; then
  echo "running full test through slurm"
  ./make_system.py build input-full.sc
  #./make_system.py post  input-full.sc
else 
  echo "running short test"
  #./make_system.py build input-full.sc
  #./make_system.py resume input-full.sc
  ./make_system.py post  input-full.sc
  #./make_system.py slip  input-full.sc
fi


