#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=fitModel_secnd_sva_strat
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-00:30:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/output_model_training_secnd_sva_strat.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_training/error_model_training_secnd_sva_strat.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
data="../../results/preproc/second_round/merged_counts_mod_alpha1_coefs_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
metDat="../../results/parsed/merged/merged_metdat.csv"
train_test="../../results/parsed/merged/train_test.csv"
alpha="1"
mem="24G"
preFiltGenes="none"
ext_substudy="brainSeq_pI"
braakThrshld="4"
CVfolds=10
#outDir="../../results/models/first_round_sva_strat_wBSPIext_qnormAllTogether/"
outDir="../../results/models/secnd_round_sva_strat_wBSPIext/"

# Train the model
########################################################################################################################
Rscript mod_train_and_test.R $data \
    --metDat $metDat \
    --trainSetFile $train_test \
    --alpha $alpha \
    --mem $mem \
    --preFiltGenes $preFiltGenes \
    --lambda \
    --ext_substudy $ext_substudy \
    --braakThrshld $braakThrshld \
    --CVfolds $CVfolds \
    --stratifiedCV \
    --overwrite_mod \
    --outDir $outDir #\
    #--refit

#modFile="${outDir}modFuncsAlpha${alpha}.rds"

#Rscript model_getCoefs.R $modFile --outDir $outDir

conda deactivate