#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=fitModel_first
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-00:30:00
#Define sdout path
#SBATCH --output=output_model_training_first.txt
#Define sderr path
#SBATCH --error=error_model_training_first.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
data="../../results/preproc/first_round/merged_counts_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
metDat="../../results/parsed/merged/merged_metdat.csv"
train_test="../../results/parsed/merged/train_test.csv"
alpha="1"
mem="24G"
preFiltGenes="none"
ext_substudy="brainSeq_pI"
braakThrshld="4"
CVfolds=10
outDir="../../results/models/first_round/"

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
    --outDir $outDir \
    --refit

conda deactivate