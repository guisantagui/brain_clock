#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=parse_perts
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 8
#SBATCH --time=01-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/parse_perts/output_parse_perts.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/parse_perts/error_parse_perts.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

########################################################################################################################
# Brain Clock - parse_perts module: Perform integration of the perturbation datasets.                                  #
########################################################################################################################

# 1. Parsing the compilation of perturbations
conda activate r-4.3.1
Rscript parse_pert_compilation.R
conda deactivate

# 2. Download, extraction and parsing LINCS L1000 data
bash generate_lincs_matrix.sh

# 3. Merge LINCS L1000 and compilation of perturbation data
conda activate r-4.3.1
Rscript merge_perts.R

# 4. Normalization of the merged perturbation dataset
pert_file="../../results/parsed/merged_perts/merged_perts_exprsn.csv"
clin_file="../../results/preproc/second_round/merged_counts_mod_alpha1_coefs_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
train_test_file="../../results/parsed/merged/train_test.csv"

Rscript norm_perts_data.R $pert_file --clin_dat $clin_file --train_test $train_test_file

conda deactivate