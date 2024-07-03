#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=FPKM2TPM
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10GB
#SBATCH -c 1
#SBATCH --time=00-06:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/output_ROSMAP_FPKM2TPM.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/error_ROSMAP_FPKM2TPM.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

ROSMAP_FPKM="/home/users/gsantamaria/projects/brain_clock/data/ROSMAP_RNAseq_FPKM_gene.tsv"

# Convert ROSMAP FPKM to TPM
Rscript ROSMAP_FPKM2TPM.R $ROSMAP_FPKM

conda deactivate