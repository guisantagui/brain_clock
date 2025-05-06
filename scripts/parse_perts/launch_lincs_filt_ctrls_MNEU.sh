#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=lincsFiltMNEU
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH -c 8
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/output_lincs_filtMNEU.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/lincs_obtain/error_lincs_filtMNEU.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
gctxFile="/mnt/lscratch/users/gsantamaria/test_large_files/ctl_predicted_RNAseq_profiles.gctx"
cellType="MNEU"
batchSize=8000
outDir="/mnt/lscratch/users/gsantamaria/test_large_files/MNEU/gctx/"

# Run
########################################################################################################################
Rscript lincs_filt.R $gctxFile --cellType $cellType --batchSize $batchSize --outDir $outDir

conda deactivate