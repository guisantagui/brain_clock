#!/bin/bash

perts="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/model_pert_tests/modAllGenes_integWLincs_sva_oAge_chronAge_alph1/pred_ages.csv"
respVar="age_chron"
metDat="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/integ_LINCSSamps_all/combined_metDat_wTBI_wPert111_wLINCS.csv"

Rscript model_perts_stats.R $perts --respVar $respVar --metDat $metDat