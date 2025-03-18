# Brain clock: Code for the analyses performed in the study "A machine-learning approach identifies rejuvenating interventions in the human brain"

This repository contains the scripts used for preprocessing the transcriptomics data we used to train and test the transcriptomic brain clock that was used in our study "A machine-learning approach identifies rejuvenating interventions in the human brain", as well as to perform the rest of the analyses described in the study. This trained model underlies our platform for detection of brain-rejuvenating interventions, [`brainAgeShiftR`](https://gitlab.lcsb.uni.lu/CBG/brainAgeShiftR). The scripts within the scripts directory are divided into the following modules:
1. Dataset parsing (`scripts/dataset_parsing`)
2. Preprocessing (`scripts/preprocessing`).
3. Model training (`scripts/model_training`).
4. Model testing (`scripts/model_testing`).
5. Mice analyses (`scripts/mice_analyses`).
6. Functional enrichment (`scripts/func_enrich`).

## Preprocessing the data and training the model
To train the model and reproduce our analyses, the following steps need to be performed
## 1. Dataset parsing
In this module, bulk RNA-seq datasets from several sources and their metadata were manually curated in order to homogenize the gene nomenclature and metadata variables, and were merged into the `data.csv` and `metadata.csv` files the next steps use as input.  
