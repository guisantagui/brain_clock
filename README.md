# Brain clock: Code for the analyses performed in the study "A machine-learning approach identifies rejuvenating interventions in the human brain"

This repository contains the scripts used for preprocessing the transcriptomics data we used to train and test the transcriptomic brain clock that was used in our study "A machine-learning approach identifies rejuvenating interventions in the human brain", as well as to perform the rest of the analyses described in the study. This trained model underlies our platform for detection of brain-rejuvenating interventions, [`brainAgeShiftR`](https://gitlab.lcsb.uni.lu/CBG/brainAgeShiftR). The scripts within the `scripts` directory are structured into the following modules:
- Clinical data curation and parsing ([`scripts/parse_clinc`](scripts/parse_clinc))
- Preprocessing ([`scripts/preprocessing`](scripts/preprocessing)).
- Model training ([`scripts/model_training`](scripts/model_training)).
- Perturbation data curation and parsing [`scripts/parse_perts`](scripts/parse_perts).
- Model testing ([`scripts/model_testing`](scripts/model_testing)).
- Mice analyses ([`scripts/mice_analyses`](scripts/mice_analyses)).
- Functional enrichment ([`scripts/func_enrich`](scripts/func_enrich)).

## Preprocessing the data and training the model
To train the model and reproduce our analyses, the following steps need to be performed in this exact order:
### 1. Clinical data curation and parsing
In this module, brain bulk RNA-seq datasets from several sources and their metadata are manually curated in order to homogenize the gene nomenclature and metadata variables, and are merged into the `merged_counts.csv` and `merged_metdat.csv` files the next step uses as input. To perform this step run the following commands:
```bash
cd scripts/parse_clinc
sbatch dataset_parsing.sh
cd ../..
```
More information on the curation and parsing process can be found in [`scripts/parse_clinc`](scripts/parse_clinc).
### 2. First training round
As mentioned in the **Methods** section of our paper, the training of the model was performed in two steps:
1. A first round performed on a preprocessed matrix with all the genes intersected with the integrated datasets.
2. A second round performed on a matrix that was filtered to keep the genes that were selected in the first fitting round before preprocessing was performed.
Therefore, a fitting round was divided in preprocessing and fittinf. 
#### 2.1 Preprocessing
Here, the whole `merged_counts.csv` matrix is preprocessed. More information on the preprocessing process can be found in [`scripts/preprocessing`](scripts/preprocessing).
For running the first preprocessing round the following commands need to be executed:
```bash
cd scripts/preprocessing
sbatch launch_preproc_pipe.sh
cd ../..
```
Results of this step are saved in `results/preproc/first_round/`. This pipeline generates several PCA plots of the steps outlined in [`scripts/preprocessing`](scripts/preprocessing), some intermediate `rds` files used for batch removal, and the `merged_counts_log2_qnorm_noCerebell_onlyAge_svaAdj.csv` file that will be used for the first round of training. 
#### 2.2 Model training
Train the model on the `merged_counts_log2_qnorm_noCerebell_onlyAge_svaAdj.csv` stored into `results/preproc/first_round/`. Tests the model on samples suffering from neurodegeneration, and computes significance tests on the differences of predicted ages between healthy individuals and neurodegeneration-affected individuals. For doing so, run the following commands:
```bash
cd scripts/model_training
sbatch fitModel_first_round.sh
cd ../..
```
Results are saved in `results/models/first_round/`. They consist on plots of the performance metrics of the model, plots and `TXT` files of the comparisons between samples originating from healthy individuals and from individuals affected by neurodegeneration, the `h2o` model object, stored in `results/models/first_round/mod_alpha1/`, and the coefficient file `mod_alpha1_coefs.csv`. The latter is used to filter genes during the 2nd round of preprocessing.
### 3. Second training round
#### 3.1 Preprocessing
Run:
```bash
cd scripts/preprocessing
sbatch launch_preproc_pipe_filtsigngenes.sh
cd ../..
```
Analogously to 2.1 section, preprocessed data for the second round of training will be saved into `results/preprocessing/second_round/`.
#### 3.2 Model training
Train the model on the `data_noCerebell_onlyAge_svaAdj.csv` stored into `results/preprocessing/second_round/`. Run:
```bash
cd scripts/model_training
sbatch fitModel_second_round.sh
cd ../..
```
Results of model training will be saved in `results/models/secnd_round/`.

### 4. Perturbation data curation and parsing
The perturbation data used in this study consisted in a compendium of in vitro perturbations individual cell types, whose accession numbers are indicated in the Table S5 of our publication, and LINCS L1000 Level 3 perturbation data. In this step, we filtered the perturbation compendium to include only cell types found in the central nervous system, downloaded LINCS perturbation and controls files, extracted brain cell types, merged the matrices and normalized them. Run:
```bash
cd scripts/parse_perts
Rscript launch_parse_perts.sh
cd ../..
```
Results of this step will be saved in `results/parsed/merged_perts/`. More information on this step can be found in [`scripts/parse_perts`](scripts/parse_perts).
### 5. Model testing
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