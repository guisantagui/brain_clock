#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=svaModFit
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 1
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/sva_model/output_sva_model_fit_allButLINCS.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/sva_model/error_sva_model_fit_allButLINCS.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

inFile="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast_excldLINCS/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_noCerebell.csv"
svobj="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast_excldLINCS/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_noCerebell_onlyAge_svobj.rds"
outDir="/home/users/gsantamaria/projects/brain_clock/results/svaMod/allButLINCS/"

Rscript sva_model_fit.R $inFile --svobj $svobj --outDir $outDir

conda deactivate