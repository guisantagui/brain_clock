# Perturbation data curation and parsing
In this module, compendium of perturbation expression datasets was curated and formatted, the expresion data corresponding to central nervous system cells was extracted from LINCS L1000 Level 3 `GCTX` files, and they were combined into single expression and metadata matrices. Finally, the resulting expression matrix was quantile-normalized using as target distribution the one extracted from the SVA-processed clinical matrix.
For running this module, the following command needs to be executed:
```bash
sbatch launch_parse_perts.sh
```

The folloging are the detailed steps of this module.

## 1. Parsing the compilation of perturbations
Here, the datasets we compiled from several sources, which are indicated in the Supplementary Table 5 of our manuscript. This step is performed by `parse_pert_compilation.R`, and requires the following files, all of them already included in `./data/perturbation/`:
- `111_sascha_complete_counts_df.csv`
- `javier_111_counts.csv`
- `Sascha_metadata.csv`
- `javier_metadata.csv`

The output files are `counts_111.csv` and `metDat_111.csv`, both saved at `./results/parsed/`.

## 2. Download, extraction and parsing LINCS L1000 data
Here, GTCX files of LINCS L1000 Level 3 are downloaded from [here](https://lincs-dcic.s3.amazonaws.com/LINCS-data-2020/RNA-seq/cp_predicted_RNAseq_profiles.gctx) and [here](https://lincs-dcic.s3.amazonaws.com/LINCS-data-2020/RNA-seq/ctl_predicted_RNAseq_profiles.gctx), and expression data and metadata corresponding to brain cell types are extracted from them in batches. Finally, they are concatenated into a single expression and metadata files per cell type. This step is performed `generate_lincs_matrix.sh`.

## 3. Merge LINCS L1000 and compilation of perturbation data
Here, the perturbation compilation was combined with the cell-type matrices extracted from the GCTX files during the previous step. The metadata files are homogenized to follow the same structure as the clinical metadata. This step is performed by `merge_perts.R`, and the resulting files are `merged_perts_exprsn.csv` and `merged_perts_metdat.csv`, both saved at `./parsed/merged_perts/`. The files required by this script are the following:
- LINCS concatenated expression matrices and metadata files of NPC (Neural Progenitor Cells), NEU (neurons) and MICROGLIA-PSEN1, which should have been placed by `generate_lincs_matrix.sh` in the corresponding subdirectory within `./data/perturbation/lincs/`.
- `counts_111.csv` and `metDat_111.csv`, which should have been placed in `./results/parsed/` by `parse_pert_compilation.R`.
- `merged_metdat.csv`, which should have been placed in `./results/parsed/merged/` during the execution of the [dataset_parsing module](./scripts/dataset_parsing).

## 4. Normalization of the merged perturbation dataset
Here, the combined dataset obtained is quantile-normalized, using the distribution of the SVA-processed clinical data (generated during [dataset_parsing](./scripts/dataset_parsing) module) that was used as input for second round of training as a target for normalization (stored in `./results/preproc/second_round/`). For doing so, we only take into account the samples assigned to the training set. This step is performed by `norm_perts_data.R`. The resulting file is `merged_perts_exprsn_qnorm.csv`, and is saved at `./parsed/merged_perts/`.