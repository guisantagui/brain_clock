#!/bin/bash


data="/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_noCerebell_onlyAge_svaAdj.csv"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS.csv"
respVar="age_chron"
ageTransPars="/home/users/gsantamaria/projects/brain_clock/data/for_model_files/GompertzMakehamParameters.rds"
alpha="1"
mem="24G"
braakThrshld="4"
outDir="/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWLincs_and_sc_sva_chron_age/"


#Rscript mod_train_and_test.R $data --metDat $metDat --respVar $respVar --ageTransPars $ageTransPars --alpha $alpha --mem $mem --braakThrshld $braakThrshld --outDir $outDir

modFile="${outDir}modFuncsAlpha${alpha}.rds"

Rscript model_getCoefs.R $modFile --outDir $outDir