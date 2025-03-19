# Brain clock: Code for the analyses performed in the study "A machine-learning approach identifies rejuvenating interventions in the human brain"

This repository contains the scripts used for preprocessing the transcriptomics data we used to train and test the transcriptomic brain clock that was used in our study "A machine-learning approach identifies rejuvenating interventions in the human brain", as well as to perform the rest of the analyses described in the study. This trained model underlies our platform for detection of brain-rejuvenating interventions, [`brainAgeShiftR`](https://gitlab.lcsb.uni.lu/CBG/brainAgeShiftR). The scripts within the `scripts` directory are structured into the following modules:
- Dataset parsing ([`scripts/dataset_parsing`](scripts/dataset_parsing))
- Preprocessing ([`scripts/preprocessing`](scripts/preprocessing)).
- Model training ([`scripts/model_training`](scripts/model_training)).
- Model testing ([`scripts/model_testing`](scripts/model_testing)).
- Mice analyses ([`scripts/mice_analyses`](scripts/mice_analyses)).
- Functional enrichment ([`scripts/func_enrich`](scripts/func_enrich)).

## Preprocessing the data and training the model
To train the model and reproduce our analyses, the following steps need to be performed
## 1. Dataset parsing
In this module, bulk RNA-seq datasets from several sources and their metadata were manually curated in order to homogenize the gene nomenclature and metadata variables, and were merged into the `data.csv` and `metadata.csv` files the next steps use as input.  

## 2. First training round
As mentioned in the **Methods** section of our paper, the training of the model was performed in two steps:
1. A first round performed on a preprocessed matrix with all the genes intersected with the integrated datasets.
2. A second round performed on a matrix that was filtered to keep the genes that were selected in the first fitting round before preprocessing was performed.
Therefore, a fitting round was divided in preprocessing and fittinf. 
### 2.1 Preprocessing
Here, the whole `data.csv` matrix was preprocessed. More information on the preprocessing process can be found in [`scripts/preprocessing`](scripts/preprocessing).
For running the first preprocessing round the following commands need to be executed:
```bash
cd scripts/preprocessing
sbatch launch_preproc_pipe.sh
cd ../..
```
Results will be saved in `results/preprocessing/first_round/`. This pipeline generates several PCA plots of the steps outlined in [`scripts/preprocessing`](scripts/preprocessing), some intermediate `rds` files used for batch removal, and the `data_noCerebell_onlyAge_svaAdj.csv` file that will be used for the first round of training. 
### 2.2 Model training
Train the model on the `data_noCerebell_onlyAge_svaAdj.csv` stored into `results/preprocessing/first_round/`. Tests the model on samples suffering from neurodegeneration, and computes significance tests on the differences of predicted ages between healthy individuals and neurodegeneration-affected individuals. For doing so, run the following commands:
```bash
cd scripts/model_training
sbatch fitModel_first_round.sh
cd ../..
```
Results will be saved in `results/models/first_round_model/`. They consist on plots of the performance metrics of the model, plots and `TXT` files of the comparisons between samples originating from healthy individuals and from individuals affected by neurodegeneration, the `h2o` model object, stored in `results/models/first_round_model/modFuncsAlpha1/`, and the coefficient file `modFuncsAlpha1_coefs.csv`. The latter will be used to filter genes during the 2nd round of preprocessing.
## 3. Second training round
### 3.1 Preprocessing
Run:
```bash
cd scripts/preprocessing
sbatch launch_preproc_pipe_filtsigngenes.sh
cd ../..
```
Analogously to 2.1 section, preprocessed data for the second round of training will be saved into `results/preprocessing/second_round/`.
### 3.2 Model training
Train the model on the `data_noCerebell_onlyAge_svaAdj.csv` stored into `results/preprocessing/second_round/`.
```bash
cd scripts/model_training
sbatch fitModel_second_round.sh
cd ../..
```
Results of model training will be saved in `results/models/second_round_model/`.
## 4. Model testing
Here we will use the final model to compute the transcriptional ages of:
- The perturbation data.
- The single cell pseudobulk data.
These samples are included in `results/preprocessing/second_round/data_noCerebell_onlyAge_svaAdj.csv`. After computing the transcriptional ages, this pipeline computes significance tests for differences in transcriptional ages either between control and perturbed samples, and between young and old samples.
To run the model testing run:
```bash
cd scripts/model_testing
# Perturbation analysis
sbatch launch_model_test_perts.sh
# Single cell pseudo-bulk analysis
sbatch launch_model_test_sCell.sh
cd ../..
```
Results will be saved in `results/model_testing/`.