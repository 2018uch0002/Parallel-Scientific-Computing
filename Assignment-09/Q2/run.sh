#!/bin/bash -x
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH -p el8
#SBATCH --time=00:10:00

spack load gcc@12.4.0
module load spectrum-mpi

# abeille input file
echo "Job started at: $(date)"
FOLDER=/gpfs/u/home/PCM2/PCM2nghp/scratch/Assignment_09/

mpirun -np $SLURM_NPROCS /gpfs/u/home/PCM2/PCM2nghp/scratch/Assignment_09/helloword