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
#SBATCH --output=/home/users/gsantamaria/projects/brain_clock/scripts/model_testing/output_model_test_perts.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/brain_clock/scripts/model_testing/error_model_test_perts.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# This version takes the integrated dataset/model with all brain cells in LINCS (NEU, NPCs and MIC).

# Variables for the pipeline
########################################################################################################################
datFile="../../results/parsed/merged_perts/merged_perts_counts_quantNorm.csv"
modDir="../../results/models/first_round_sva_strat_wBSPIext_qnormAllTogether/mod_alpha1/"
metDat="../../results/parsed/merged_perts/merged_perts_metdat.csv"
respVar="age_chron"
#ageTransPars="../../data/for_model_files/GompertzMakehamParameters.rds"
batchSize=4000
whatSampsTest="perturbation"
mem="50G"
outDir="../../results/model_test_perts/mod_first_round_sva_strat_wBSPIext_qnormAllTogether/"

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

conda deactivate