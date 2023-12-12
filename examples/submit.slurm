#!/bin/bash -l
#SBATCH -J main_gpu
#SBATCH -t 02-00:00:00
#SBATCH -o slurm.%x.%j.out
#SBATCH -e slurm.%x.%j.err
#SBATCH --nodes 1
#SBATCH --partition xeon-g6-volta
#SBATCH --exclusive

module load anaconda
source activate gpupyscf

cd /home/gridsan/jburns/kinetics/hamdb/examples

python main.py
