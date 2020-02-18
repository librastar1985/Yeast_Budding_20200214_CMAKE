#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --output=YB_cellwall3_l0d5b0d05a0d1_detgrowth_strain0d2_maxt400_esfreq200_longsim400_# This affects the print out of the "std::cout" in the script, make sure this is changed for different jobs.
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="0500501"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load extra
module load GCC
module load cuda/9.1


srun -p gpu --gres=gpu:1 ./virus-model -dt=0.001 Data_structure.xml

