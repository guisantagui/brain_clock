#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=scGeneDirChngs
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-12:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_testing/output_scGeneDirChngs.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_testing/error_scGeneDirChngs.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

# Run the simulations
########################################################################################################################

Rscript get_sc_gene_dir_changes_prop.R

conda deactivate