sc_mod_coef_file <- "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_SCmod_chron_age/modFuncsAlpha1_coefs.csv"

full_mod_coef_file <- "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_chron_age/modFuncsAlpha1_coefs.csv"

sc_mod_coefs <- read.csv(sc_mod_coef_file)
full_mod_coefs <- read.csv(full_mod_coef_file)

merged_df <- data.frame(ensembl_gene_id = unique(c(full_mod_coefs$ensembl_gene_id,
                                                   sc_mod_coefs$ensembl_gene_id)))

write.csv(merged_df, file = "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_SCmod_chron_age/fullMod_andSCMod_genes.csv")