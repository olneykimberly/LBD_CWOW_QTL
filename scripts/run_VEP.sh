#!/bin/bash
#SBATCH --job-name=VEP                                                              
#SBATCH --time=7-00:00:00                               
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH -o slurm.VEP.job.%j.out
#SBATCH -e slurm.VEP.job.%j.err

# activate conda environment
source $HOME/.bash_profile
module load python3
conda activate QTL

# change directory to where Snakefile is located
cd /tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/snp_array/reference
/tgen_labs/jfryer/kolney/tools/ensembl-vep/./vep -i snp_rs.txt --cache

