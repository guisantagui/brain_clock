#!/bin/bash

perts="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/model_pert_tests/modLndGenes_integWLincs_sva_oAge/pred_ages.csv"
metDat="/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/integ_LINCSSamps_all/combined_metDat_wTBI_wPert111_wLINCS.csv"

Rscript model_perts_stats.R $perts --metDat $metDat