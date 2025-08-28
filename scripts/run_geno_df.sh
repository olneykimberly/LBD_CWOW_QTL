#!/bin/bash
#SBATCH --job-name=transpose                         
##SBATCH --tasks=1     
#SBATCH --cpus-per-task=4                                   
#SBATCH --time=48:00:00 # 8 hours   
#SBATCH -n 8 # threaded 
#SBATCH --mem=500G 
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# source your bach profile to get your specific settings  
source $HOME/.bash_profile

module load python3
conda activate QTL 

python get_geno_df.py

