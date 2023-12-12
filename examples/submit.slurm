#!/bin/bash -l
#SBATCH -J main
#SBATCH -t 02-00:00:00
#SBATCH -o slurm.%x.%j.out
#SBATCH -e slurm.%x.%j.err
#SBATCH --nodes 8
#SBATCH --exclusive

module load anaconda mpi
source activate mpipyscf

cd /home/gridsan/jburns/kinetics/hamdb/examples

mpirun -np 8 python main.py
