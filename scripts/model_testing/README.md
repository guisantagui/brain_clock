# Model testing
This directory contains the scripts for testing the model on:
- The chemical perturbation data, which comes from Level 3 of [**LINCS L1000**](https://maayanlab.cloud/sigcom-lincs/#/Download).
- The pseudo-bulk single cell data, derived from the brain samples included in [**ageAnno**](https://relab.xidian.edu.cn/AgeAnno/#/).
The analyses performed in this module consist in:
1. Computing the transcriptional age.
2. Assessing the significance of transcriptional age shifts.
For the commands needed for performing the analyses relative to this step are the following:
```bash
# For perturbation data
sbatch launch_model_test_perts.sh
# For single cell data
sbatch launch_model_test_sCell.sh
```
## Outputs
Results will be saved in `./results/model_testing/perts/` and `./results/model_testing/sCell/` for perturbations and single cell, respectively, and will consist on:
- `pred_ages.csv`: a CSV file with the transcriptional ages of each one of the samples.
- `pred_ages_stats.csv`: a CSV file with the p-values and log2 fold changes of the comparisons, which consist in:
    - Treated vs un-treated for the perturbation data
    - Old (age >= 70) vs young (age <= 30) for the pseudo-bulk of the single cell data.
- Plots representing the results.