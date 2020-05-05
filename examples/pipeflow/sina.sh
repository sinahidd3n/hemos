#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --job-name=sina_job
#SBATCH --mem=1000
#SBATCH --output=/scratch/m/maftoon/sanvarin/hemo1/examples/pipeflow/hemo_5.out
#SBATCH --error=/scratch/m/maftoon/sanvarin/hemo1/examples/pipeflow/hemo_5.err
#SBATCH --time=00:15:00



module load gcc
module load hdf5
module load cmake
module load openmpi
srun mpirun -np 40 ./pipeflow config.xml
