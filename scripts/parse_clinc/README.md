# Clinical data curation and parsing
This directory contains the code used for parsing the different clinical datasets into unified metadata and counts matrices. This step in the pipeline was rather manual: we inspected the metadata provided by each substudy, categorized the samples as controls based on the available variables, and homogenized the units and categories included across the different datasets. The datasets integrated were the following:
- [**The RNAseq Harmonization Study (rnaSeqReprocessing)**](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn9702085)
- [**The Living Brain Project (LBP)**](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn26337520)
- [**The Aging, Dementia and Traumatic Brain Injury Study(TBI)**](https://aging.brain-map.org)
- [**The Genotype-Tissue Expression (GTEx) project (GTEx)**](https://www.gtexportal.org/home/downloads/adultgtex/bulk_tissue_expression)
- [**BrainSeq Phase 1**](http://eqtl.brainseq.org/phase1/)
- [**BrainSeq Phase 2**](https://eqtl.brainseq.org/phase2/)
- Pseudobulk data of brain samples from [**ageAnno**](https://relab.xidian.edu.cn/AgeAnno/#/)

To run this module the following command needs to be run:
```bash
sbatch dataset_parsing.sh
```

## Steps in the module
This manual parsing can be divided in the following steps:

### 1. Parsing of AMP-AD datasets
#### 1.1. Parsing of AMP-AD individual metadata
Here we explored diagnosis variables to determine which individuals are healthy, which suffered from neurodegeneration and which needed to be removed due to having ambiguous or non-existent age information or diagnosis out of the scope of our study. In this step we homogenized information such as ethnicity and sex accross the different metadata datasets. This step is performed by `amp_ad_metadata_parsing.R`. To run, this script needs to have the following files, whose access need to be requested in synapse.org:
- `LBP_individual_metadata.csv` (syn58589260), obtained from [The Living Brain Project (LBP) metadata directory](https://www.synapse.org/Synapse:syn27127033). This file has to be placed in `./data/metadata/LBP`.
- `RNAseq_Harmonization_Mayo_combined_metadata.csv` (syn27000373), `RNAseq_Harmonization_MSBB_combined_metadata.csv` (syn27000243) and `RNAseq_Harmonization_ROSMAP_combined_metadata.csv` (syn27034471), obtained from [The RNAseq Harmonization Study (RNAseq Harmonization) metadata directory](https://www.synapse.org/Synapse:syn27000096). These files need to be placed in `./data/metadata/RNAseq_Harmonization`.
- `MayoRNAseq_individual_metadata.csv` (syn23277389), obtained from [The MayoRNAseq metadata directory](https://www.synapse.org/Synapse:syn23634010). This file needs to be placed in `./data/metadata/rnaSeqReprocessing`. The reason for including this file is that there are several samples that in the RNAseq Harminozation version have the Thal and Braak scores missing but in this version are not, so we included this metadata file to complete them.

The output of this script are a CSV metadata file for each dataset, and a metadata file comining all the selected samples from the RNAseq Harmonization study (`RNAseq_Harmonization_ind_all.csv`), all written to `./results/metadata_parsed/`.

#### 1.2. Homogenizing RNAseq harmonization metadata and counts matrices
Here the counts files from the RNAseq Harmonization studies are merged together, and the clinical information of ROSMAP study is added to the RNAseq Harmonization comined metadata file. This step is performed by `rnaHarm_parseCounts.R`. To run, it requires the files generated in the previous step, and other files obtained from synapse.org. The whole set of files needed by this script and their required location  is the following:
- The metadata files generated in the previous step.
- `ROSMAP_clinical.csv` (syn3191087), obtained from the [ROSMAP study metadata subdirectory](https://www.synapse.org/Synapse:syn3157322).
- `Mayo_gene_all_counts_matrix_clean.txt` (syn21544635), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Mayo subdirectory](https://www.synapse.org/Synapse:syn20825471). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `MSBB_gene_all_counts_matrix_clean.txt` (syn21544666), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) MSSM subdirectory](https://www.synapse.org/Synapse:syn20957610). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch1_gene_all_counts_matrix_clean.txt` (syn22283382), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 1 subdirectory](https://www.synapse.org/Synapse:syn22279877). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch2_gene_all_counts_matrix_clean.txt` (syn22301601), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 2 subdirectory](https://www.synapse.org/Synapse:syn22296752). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch3_gene_all_counts_matrix_clean.txt` (syn22314230), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 3 subdirectory](https://www.synapse.org/Synapse:syn22300974). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch4_gene_all_counts_matrix_clean.txt` (syn25817661), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 4 subdirectory](https://www.synapse.org/Synapse:syn25810549). This file needs to be placed in `./data/expression/RNAseqHarm/`.

The output of this script are a counts matrix of all the RNAseq_Harmonization study datasets (`RNAseqHarm_allCounts.csv`) and a metadata file including clinical variables (), which are saved in `./results/parsed/`.

#### 1.3. Homogenizing LBP metadata.
Here, the individual, biospecimen and assay metadata of the LBP study are combined into a single file, following the same structure as the metadata of the RNAseq harmonization study. This step is performed by `parse_LBP_metdat.R`. The files required by this script are:
- `LBP_ind_all.csv`, which should have been generated by `amp_ad_metadata_parsing.R` in step 1.1, and placed into `./results/metadata_parsed/`.
- `LBP_individual_metadata.csv` (syn58589260), `LBP_assay_RNAseq_metadata.csv` (syn64555267) and `LBP_biospecimen_metadata.csv` (syn64504156), obtained from [LBP's Metadata subdirectory](https://www.synapse.org/Synapse:syn27127033). These files need to be placed into `./data/metadata/LBP/`.

The output of this script is a metadata file of the LBP samples following the same structure of `RNAseq_Harmonization_ind_all.csv` (`LBP_metadata_unified.csv`), which is saved in `./results/parsed/`.

### 2. Parsing of GTEx datasets.
Here, the counts files originated from brain tissues belonging to GTEx were merged into a single file, the metadata files were explored to select the samples that fit the requirements of our study (having an unambiguous age and not having suffered psychiatric or drug abuse disorders), and to divide the samples in healthy controls or affected by neurodegeneration, and homogenized the file structure to resemble that of `RNAseq_Harmonization_ind_all.csv`. This step is performed by `GTEx_parsing.R`. The required files for this step are:
- GTEx brain bulk counts files. These files are open access, and can be obtained from [The GTEx portal bulk tissue expression section](https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression). These files need to be placed in `./data/gtex/counts/`. In particular, the files we used were:
    - `gene_reads_2017-06-05_v8_brain_amygdala.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_anterior_cingulate_cortex_ba24.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_caudate_basal_ganglia.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_cerebellar_hemisphere.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_cerebellum.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_cortex.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_hippocampus.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_hypothalamus.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_nucleus_accumbens_basal_ganglia.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_putamen_basal_ganglia.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_spinal_cord_cervical_c-1.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_substantia_nigra.gct.gz`
- The GTEx metadata files `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` and `GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt`, available upon request at dbGaP (protected data access information can be found [here](https://www.gtexportal.org/home/protectedDataAccess)). These files need to be placed in `./data/gtex/`.

The output of this script are a counts matrix integrating all the selected GTEx datasets (`GTEx_allCounts.csv`), and the corresponding metadata file (`GTEx_metadata_unified.csv`). Both files will be written at `./results/parsed/`.

### 3. Merge GTEx, LBP and RNAseq Harmonization datasets
Here, the parsed counts and metadata matrices coming from GTEx, LBP and RNAseq Harmonization that were generated in the previous steps were merged. This step is performed by `merge_GTEx_LBP_rnaSeqHarm.R`. The files required for this step are the following:
- `LBP_FlagshipPaper_featureCounts.csv` (syn64567831), which was obtained from the [LBP Gene Expression (RNAseq) counts subdirecotry](https://www.synapse.org/Synapse:syn52132890) and needs to be placed in `./data/expression/LBP/`.
- `LBP_metadata_unified.csv`, which should have been previously generated by `rnaHarm_parseCounts.R` and placed in `./results/parsed/` in step 1.2.
- RNAseq Harmonization counts matrix file `RNAseqHarm_allCounts.csv`, which should have been previously generated by `rnaHarm_parseCounts.R`, and placed in `./results/parsed/` in step 1.2.
- RNAseq Harmonization metadata file `RNAseq_Harmonization_ind_all_ROSMAPBtch.csv`, which should have been created by `rnaHarm_parseCounts.R` and placed in `./results/metadata_parsed/` in step 1.2.
- GTEx counts matrix file `GTEx_allCounts.csv` and its corresponding metadata file `GTEx_metadata_unified.csv`, which should have been previously generated by `GTEx_parsing.R` and placed in `./results/parsed/` in step 2.

The output of this script are combined counts and metadata matrices: `combined_counts.csv` and `combined_metDat.csv`, respecively. Both files are saved in `./results/parsed/`.

### 4. Obtention of TBI data and merge with the combined dataset
Here, the data from [the Aging, Dementia and TBI study](https://aging.brain-map.org) was downloaded, merged into a single counts matrix (keeping only samples that didn't suffered a TBI, and have non-ambiguous ages), and with the `combined_counts.csv` matrix previously generated in step 3. Metadata is processed the same way. This step is performed by `tbi_get.R`. The files required by this step are the following:
- `tbi_data_files.csv`, which is freely available [here](https://aging.brain-map.org/data/tbi_data_files.csv). This file needs to be placed in `./data/dataset_parsing/`.
- `DonorInformation.csv`, which is freely available [here](https://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]). This file needs to be placed in `./data/dataset_parsing/`.
- `combined_counts.csv` and `combined_metDat.csv`, which should have been previously generated by `merge_GTEx_LBP_rnaSeqHarm.R` in `./results/parsed/`.

The output of this script are `TBI_metadata_unified.csv`, `combined_metDat_wTBI.csv` and `combined_counts_wTBI.csv`, which are saved into `./results/parsed/`.

### 5. Obtention of BrainSeq Phase 2 data
Here, data from [BrainSeq Phase 2 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7000204/#S10) was parsed into counts and metadata matrices, following the same format of the rest of the files and including only control individuals without ambiguous ages and of 20 years old or older. This step is performed by `brainSeq_pII_parsing.R`. The only file required by this script is `rse_gene_unfiltered.Rdata`, which can be freely accessed [here](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata) and needs to be placed in `./data/brainSeq_pII/`. The outputs of this script are `brainseq_pII_counts.csv` and `brainseq_pII_metadata.csv`, the two of them written into `./results/parsed/`.

### 6. Obtention of BrainSeq Phase 1 data
Analogously to previous step, we parsed data from [BrainSeq Phase 1 study](https://doi.org/10.1038/s41593-018-0197-y) into counts and metadata matrices to fit the same format as the rest. This dataset has much less samples than the Phase 2. Because of this, we used entirely this dataset as independent validation dataset. In [BrainSeq Phase 1](http://eqtl.brainseq.org/phase1/) there are some samples that come from individuals also present in BrainSeq Phase 2. We removed those samples, and retained only the controls of 20 years old or older, finally retaining 22 samples. This step is performed by `brainSeq_pI_parsing.R`, and only requires one file: `rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda`. This file can be freely accessed [here](https://s3.us-east-2.amazonaws.com/jaffe-nat-neuro-2018/rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda), and needs to be placed in `./data/brainSeq_pI/`. The outputs of this script are `brainseq_pI_counts.csv` and `brainseq_pI_metadata.csv`, the two of them written into `./results/parsed/`.

### 7. Pre-processing, integration and annotation of AgeAnno brain scRNA-seq data
Here, the brain scRNA-seq Seurat object from [AgeAnno](https://relab.xidian.edu.cn/AgeAnno/#/) was pre-processed to remove low quality samples and cells, the different samples were integrated and the cell types annotated, and a matrix of pseudo-bulk counts per cell type within each sample was generated. This step is performed by `ageAnno_preproc.R`. The files required by this script are the following:
- `brain.rds`, which is the Seurat Object of AgeAnno's brain samples and can be downloaded [here](https://drive.google.com/drive/folders/160i9KmFJ0tEYP2QBT5IjbJqhqvuu5rBW). This file needs to be placed in `./data/ageAnno/`
- `age_info.csv`, which is the age of death of each one of the donors from AgeAnno samples. This file is already included in data/`ageAnno`.
- `CellCycleGenes_Human.csv`, which are the cell cycle genes obtained from [HBC Training](https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html). This file is already included in `./data/utility_files/`.
- `scRNAmarker.txt`, which are the cell markers used to annotate each one of the cell types in the brain in [AgeAnno's original publication](https://doi.org/10.1093/nar/gkac847) and can be obtained from [AgeAnno's GitHub page](https://github.com/vikkihuangkexin/AgeAnno?tab=readme-ov-file). This file is already included in `./data/utility_files/`.

The outputs of this script are be written in `./results/parsed/ageAnno/`, and consist in several plots of the pre-processing steps, seurat objects for each sample and for the integrated objects, plots of the annotation process and the two files that will be used in downstream analyses, which are:
- `ageanno_brain_pb_counts.csv`: the pseudo-bulk counts of each cell type per sample.
- `ageanno_brain_pb_metdat.csv`: the parsed metadata matching pseudo-bulk rows.

### 8. Integrate all clinical samples into a single file and pre-process
#### 8.1 Add brainSeq Phase II, AgeAnno and perturbation compilation to merged dataset
In this step, the counts and metadata datasets generated in steps 4, 5, 6 and 7 are merged with `combined_metDat_wTBI.csv` and `combined_counts_wTBI.csv`, which were generated in step 4. This step is performed by `merge_clinical.R`. If previous step have been ran correctly, the required files should be already placed in the corresponding directories. The output files are `merged_counts.csv` and `merged_metdat.csv`, both of them saved at `./results/parsed/merged/`.
#### 8.2 Filter clinical dataset
Here the merged counts and metadata datasets are filtered by:
- Keeping only genes present in LINCS L1000 level 3
- Removing genes with zero or near zero expression
- Removing samples with RIN < 6
This step is performed by `filt_dataset.R`. Unless specified, `merged_counts.csv` and `merged_metdat.csv` will be overwritten.