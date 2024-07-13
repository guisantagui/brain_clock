#!/bin/bash

perts="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/model_pert_tests/modAllGenes_ingegWLincs_combat/pred_ages.csv"
metDat="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/integ_LINCSSamps_all/combined_metDat_wTBI_wPert111_wLINCS.csv"

Rscript model_perts_stats.R $perts --metDat $metDat