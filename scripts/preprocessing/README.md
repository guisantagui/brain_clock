# 2. Preprocessing
This directory contains the scripts used for preprocessing the RNA-seq data that was used for the training and testing of the brain clock. The input files consist in `data.csv`, a quantile-normalized log transformed counts (`quant_norm(log2(counts + 1))`) CSV matrix (genes in columns, samples in rows), and `metadata.csv`, a metadata CSV file. `data.csv` consists of samples from the rnaseqreprocessing, LBP, TBI and GTEx studies, pseudobulk data of brain samples from ageAnno, and brain cell types controls and perturbation data from LINCS1000, as well as several perturbation datasets collected from several sources (indicated in the Table S5). This matrix was generated in the module 1 (Dataset parsing), which already hade some preprocessing needed for merging the datasets. The steps carried out in this module are:
1. Plot PCA to evaluate batch effects and if there are sub-tissues of the brain with significant differences to others.
2. Removal of samples originating from the cerebellum, as we noticed that they have very different transcriptional profiles.
3. Plot PCA without the cerebellum samples.
4. Remove batch effects using surrogate variable analysis (SVA), via regressing out the surrogate variables.
5. Plot PCA after the SVA to see if batch structure is removed from the data.
6. Compute surrogate variables (with SVA again) on the SVA-adjusted dataset and plot the surrogate variables to ensure batch structure has been removed.
`launch_preproc_pipe.sh` coordinates the outlined preprocessing steps. The pipeline can be launched in a Slurm-based HPC system by running:

```bash
sbatch launch_preproc_pipe.csv
```
As indicated in the methods, our model was fit in two different steps. In the first step, we identified 664 genes with non-zero coefficients. `launch_preproc_pipe_filtsigngenes.csv` performs the same preprocessing steps than `launch_preproc_pipe.csv`, but with a prior filtering step to keep only the 664 genes selected by the first fitting. The dataset resulting from `launch_preproc_pipe.csv` is the one that was used used to perform the second training round, which selected 431 genes among the 664 initially provided. By this two-step process two things are accomplished:
- We obtain a model 35% smaller, with a really slight decline in predictive performance.
- We obtain a set of surrogate variables from the 664 genes obtained from the first fitting, which allowed us to train a frozen SVA model to correct for inferred surrogate variables based on the expression of these 664 genes, which is implemented in our package `brainAgeShiftR` and can be applied with the function `do_frozenSVA()`.