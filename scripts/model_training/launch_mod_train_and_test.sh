#!/bin/bash


data="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/integ_LINCSSamps_all/combined_counts_wTBI_wPert111_log2_quantNorm_preproc_wLINCS_noCerebell_combat.csv"
metDat="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/integ_LINCSSamps_all/combined_metDat_wTBI_wPert111_wLINCS.csv"
ageTransPars="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/GompertzMakehamParameters.rds"
alpha="0.5"
mem="24G"
braakThrshld="4"
outDir="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/modAllGenes_ingegWLincs/"


Rscript mod_train_and_test.R $data --metDat $metDat --ageTransPars $ageTransPars --alpha $alpha --mem $mem --braakThrshld $braakThrshld --outDir $outDir