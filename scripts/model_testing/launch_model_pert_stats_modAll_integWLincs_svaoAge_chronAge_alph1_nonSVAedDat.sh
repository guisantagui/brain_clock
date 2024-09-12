#!/bin/bash

perts="/home/users/gsantamaria/projects/brain_clock/results/model_test_perts/modAllGenes_integWLincs_and_sc_oAge_chronAge_alph1_onNonSVAedDat/pred_ages.csv"
respVar="age_chron"
metDat="/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS_NPC_NEU_MIC.csv"

Rscript model_perts_stats.R $perts --respVar $respVar --metDat $metDat