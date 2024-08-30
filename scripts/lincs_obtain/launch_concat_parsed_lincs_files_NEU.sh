#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=lincsConcatNEU
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 4
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/output_concat_parsed_lincs_files_NEU.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/error_concat_parsed_lincs_files_NEU.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
parsLincsDir="/mnt/lscratch/users/gsantamaria/test_large_files/NEU/parsed_mats/"
cellType="NEU"
outDir="/mnt/lscratch/users/gsantamaria/test_large_files/NEU/parsed_mats/"

# Run
########################################################################################################################

Rscript concat_parsed_lincs_files.R $parsLincsDir --cellType $cellType --whichDF "expMat" --outDir $outDir &
Rscript concat_parsed_lincs_files.R $parsLincsDir --cellType $cellType --whichDF "metDat" --outDir $outDir &
wait

conda deactivate