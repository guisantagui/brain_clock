#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=bigPCA
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 4
#SBATCH --time=00-01:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/output_bigPCA_exe.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/error_bigPCA_exe.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

input="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/combined_postMergeTPMs_log_varClean.csv"
nPCs=20
outDir="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/"

# Convert ROSMAP FPKM to TPM
Rscript bigPCA_exe.R $input --nPCs $nPCs --stand --outDir $outDir

conda deactivate