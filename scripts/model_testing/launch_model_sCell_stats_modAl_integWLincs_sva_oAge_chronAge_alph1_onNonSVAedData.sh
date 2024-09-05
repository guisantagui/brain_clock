#!/bin/bash

preds="/home/users/gsantamaria/projects/brain_clock/results/model_test_sCell/modAllGenes_integWLincs_and_sc_oAge_chronAge_alph1_onNotSVAedData/pred_ages.csv"
respVar="age_chron"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS.csv"

Rscript model_sCell_stats.R $preds --respVar $respVar --metDat $metDat