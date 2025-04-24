# Dataset parsing
This directory contains the code used for parsing the different datasets into unified metadata and quantile-normalized counts matrices. This step in the pipeline was rather manual: we inspected the metadata provided by each substudy, categorized the samples as controls based on the available variables, and homogenized the units and categories included across the different datasets.
This manual parsing can be divided in the following steps:
## 1. Parsing of AMP-AD datasets
### 1.1. Parsing of AMP-AD individual metadata
Here we explored diagnosis variables to determine which individuals are healthy, which suffered from neurodegeneration and which needed to be removed due to having ambiguous or non-existent age information or diagnosis out of the scope of our study. In this step we homogenized information such as ethnicity and sex accross the different metadata datasets. This step is performed by `amp_ad_metadata_parsing.R`. To run, this script needs to have the following files, whose access need to be requested in synapse.org:
- `LBP_individual_metadata.csv` (syn58589260), obtained from [The Living Brain Project (LBP) metadata directory](https://www.synapse.org/Synapse:syn27127033). This file has to be placed in `./data/metadata/LBP`.
- `RNAseq_Harmonization_Mayo_combined_metadata.csv` (syn27000373), `RNAseq_Harmonization_MSBB_combined_metadata.csv` (syn27000243) and `RNAseq_Harmonization_ROSMAP_combined_metadata.csv` (syn27034471), obtained from [The RNAseq Harmonization Study (RNAseq Harmonization) metadata directory](https://www.synapse.org/Synapse:syn27000096). These files need to be placed in `./data/metadata/RNAseq_Harmonization`.
- `MayoRNAseq_individual_metadata.csv` (syn23277389), obtained from [The MayoRNAseq metadata directory](https://www.synapse.org/Synapse:syn23634010). This file needs to be placed in `./data/metadata/rnaSeqReprocessing`. The reason for including this file is that there are several samples that in the RNAseq Harminozation version have the Thal and Braak scores missing but in this version are not, so we included this metadata file to complete them.
The output of this script are a CSV metadata file for each dataset, and a metadata file comining all the selected samples from the RNAseq Harmonization study (`RNAseq_Harmonization_ind_all.csv`), all written to `./results/metadata_parsed/`.
### 1.2. Homogenizing AMP-AD metadata and counts matrices
Here the counts files from the RNAseq Harmonization studies are merged together, and the individual, biospecimen and assay metadata of the LBP study are combined into a single file, following the same structure as the metadata of the RNAseq harmonization study. This step is performed by `rnaHarm_parseCounts.R`. To run, it requires the files generated in the previous step, and other files obtained from synapse.org. The whole set of files needed by this script and their required location  is the following:
- The metadata files generated in the previous step.
- `Mayo_gene_all_counts_matrix_clean.txt` (syn21544635), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Mayo subdirectory](https://www.synapse.org/Synapse:syn20825471). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `MSBB_gene_all_counts_matrix_clean.txt` (syn21544666), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) MSSM subdirectory](https://www.synapse.org/Synapse:syn20957610). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch1_gene_all_counts_matrix_clean.txt` (syn22283382), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 1 subdirectory](https://www.synapse.org/Synapse:syn22279877). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch2_gene_all_counts_matrix_clean.txt` (syn22301601), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 2 subdirectory](https://www.synapse.org/Synapse:syn22296752). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch3_gene_all_counts_matrix_clean.txt` (syn22314230), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 3 subdirectory](https://www.synapse.org/Synapse:syn22300974). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `ROSMAP_batch4_gene_all_counts_matrix_clean.txt` (syn25817661), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 4 subdirectory](https://www.synapse.org/Synapse:syn25810549). This file needs to be placed in `./data/expression/RNAseqHarm/`.
- `LBP_FlagshipPaper_featureCounts.csv` (syn64567831), obtained from [LBP Gene Expression (RNAseq) counts subdirecotry](https://www.synapse.org/Synapse:syn52132890). This file needs to be placed in `./data/expression/LBP/`.
- `LBP_individual_metadata.csv` (syn58589260), `LBP_assay_RNAseq_metadata.csv` (syn64555267) and `LBP_biospecimen_metadata.csv` (syn64504156), obtained from [LBP's Metadata subdirectory](https://www.synapse.org/Synapse:syn27127033).

The output of this script are a counts matrix of all the RNAseq_Harmonization study datasets (`RNAseqHarm_allCounts.csv`) and a metadata file of the LBP samples following the same structure of `RNAseq_Harmonization_ind_all.csv` (`LBP_metadata_unified.csv`), which are saved in `./results/parsed/`.
## 2. Parsing of GTEx datasets.
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
## 3. Merge GTEx, LBP and RNAseq Harmonization datasets
Here, the parsed counts and metadata matrices coming from GTEx, LBP and RNAseq Harmonization that were generated in the previous steps were merged. This step is performed by `merge_GTEx_LBP_rnaSeqHarm.R`. The files required for this step are the following:
- `LBP_FlagshipPaper_featureCounts.csv` (syn64567831), which was obtained from the [LBP Gene Expression (RNAseq) counts subdirecotry](https://www.synapse.org/Synapse:syn52132890) and needs to be placed in `./data/expression/LBP/`.
- `LBP_metadata_unified.csv`, which should have been previously generated by `rnaHarm_parseCounts.R` and placed in `./results/parsed/` in step 1.2.
- RNAseq Harmonization counts matrix file `RNAseqHarm_allCounts.csv`, which should have been previously generated by `rnaHarm_parseCounts.R`, and placed in `./results/parsed/` in step 1.2.
- RNAseq Harmonization metadata file `RNAseq_Harmonization_ind_all_ROSMAPBtch.csv`, which should have been created by `rnaHarm_parseCounts.R` and placed in `./results/metadata_parsed/` in step 1.2.
- GTEx counts matrix file `GTEx_allCounts.csv` and its corresponding metadata file `GTEx_metadata_unified.csv`, which should have been previously generated by `GTEx_parsing.R` and placed in `./results/parsed/` in step 2.

The output of this script 


