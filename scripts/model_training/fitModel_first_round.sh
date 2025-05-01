#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=fitModel_first_round
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-00:30:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/output_model_training_first_round.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/error_model_training_first_round.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
data="../../results/preproc/test_no_lincs/merged_counts_log2_quantNorm_noCerebell_combat_onlyAge_svaAdj.csv"
metDat="../../results/parsed/merged/merged_metdat.csv"
respVar="age_chron"
ageTransPars="/home/users/gsantamaria/projects/brain_clock/data/for_model_files/GompertzMakehamParameters.rds"
alpha="1"
mem="24G"
preFiltGenes="none"
braakThrshld="4"
outDir="../../results/models/first_round/"

# Train the model
########################################################################################################################
Rscript mod_train_and_test.R $data
    --metDat $metDat \
    --respVar $respVar \
    --ageTransPars $ageTransPars \
    --alpha $alpha \
    --mem $mem \
    --preFiltGenes $preFiltGenes \
    --lambda \
    --braakThrshld $braakThrshld \
    --refit \
    --overwrite_mod \
    --outDir $outDir

modFile="${outDir}modFuncsAlpha${alpha}.rds"

Rscript model_getCoefs.R $modFile --outDir $outDir

conda deactivate