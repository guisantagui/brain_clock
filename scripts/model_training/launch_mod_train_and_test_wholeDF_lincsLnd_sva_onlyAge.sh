#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=modTrain
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-00:30:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/output_model_training_lincsLnd_sva_oAge.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/error_model_training_lincsLnd_sva_oAge.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1


# Variables for the pipeline
########################################################################################################################
data="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_lincsLndmrk_sva_fast/combined_counts_wTBI_wPert111_log2_quantNorm_preproc_wLINCS_lincsLndmrkGenes_noCerebell_onlyAge_svaAdj.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wLINCS.csv"
ageTransPars="/home/users/gsantamaria/projects/brain_clock/data/for_model_files/GompertzMakehamParameters.rds"
alpha="0.5"
mem="50G"
braakThrshld="4"
outDir="/home/users/gsantamaria/projects/brain_clock/results/models/modLincsLnd_ingegWLincs_svaOnlyAge/"

# Train the model
########################################################################################################################
Rscript mod_train_and_test.R $data --metDat $metDat --ageTransPars $ageTransPars --alpha $alpha --mem $mem --braakThrshld $braakThrshld --outDir $outDir

conda deactivate