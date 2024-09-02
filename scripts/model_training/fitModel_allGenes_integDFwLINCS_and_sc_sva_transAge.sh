#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=modTrTrns
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50GB
#SBATCH -c 1
#SBATCH --time=00-01:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/output_modTrain_lincsAll_sva_oAge_trans.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/error_modTrain_lincsAll_sva_oAge_trans.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
data="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_noCerebell_onlyAge_svaAdj.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS.csv"
respVar="age_trans"
ageTransPars="/home/users/gsantamaria/projects/brain_clock/data/for_model_files/GompertzMakehamParameters.rds"
alpha="0.5"
mem="24G"
braakThrshld="4"
outDir="/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWLincs_and_sc_sva_trans_age/"

# Analysis
########################################################################################################################
Rscript mod_train_and_test.R $data --metDat $metDat --respVar $respVar --ageTransPars $ageTransPars --alpha $alpha --mem $mem --braakThrshld $braakThrshld --outDir $outDir

modFile="${outDir}modFuncsAlpha${alpha}.rds"

Rscript model_getCoefs.R $modFile --outDir $outDir

conda deactivate