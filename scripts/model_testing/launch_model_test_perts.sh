#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=modTstPrts
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-03:00:00
#Define sdout path
#SBATCH --output=output_model_test_perts.txt
#Define sderr path
#SBATCH --error=error_model_test_perts.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# This version takes the integrated dataset/model with all brain cells in LINCS (NEU, NPCs and MIC).

# Variables for the pipeline
########################################################################################################################
datFile="../../results/parsed/merged_perts/merged_perts_exprsn_qnorm.csv"
modDir="../../results/models/secnd_round/mod_alpha1/"
metDat="../../results/parsed/merged_perts/merged_perts_metdat.csv"
alpha=0.05
what_test="t.test"
lincs_drugs="/home/users/gsantamaria/projects/brain_clock/data/perturbation/lincs/LINCS_small_molecules.tsv"
batchSize=4000
whatSampsTest="perturbation"
mem="50G"
outDir="../../results/model_test_perts/"

# Run the simulations
########################################################################################################################

Rscript model_test.R $datFile \
    --modDir $modDir \
    --metDat $metDat \
    --sizeBatch $batchSize \
    --whatSampsTest $whatSampsTest \
    --mem $mem \
    --outDir $outDir

pert_stats_in="${outDir}pred_ages.csv"

Rscript model_perts_stats.R $pert_stats_in --metDat $metDat

pert_stats="${outDir}pred_ages_stats.csv"

Rscript examine_pert_stats.R $pert_stats \
    --alpha $alpha \
    --what_test $what_test \
    --lincs_drugs $lincs_drugs

conda deactivate