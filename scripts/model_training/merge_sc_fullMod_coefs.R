sc_mod_coef_file <- "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_SCmod_chron_age/modFuncsAlpha1_coefs.csv"

full_mod_coef_file <- "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_chron_age/modFuncsAlpha1_coefs.csv"

sc_mod_coefs <- read.csv(sc_mod_coef_file)
full_mod_coefs <- read.csv(full_mod_coef_file)

# From SC coefs, take the best 15% in terms of coefficient strength
sc_mod_coefs <- sc_mod_coefs[order(abs(sc_mod_coefs$coefficients),
                                   decreasing = T), ]

sc_mod_coefs <- sc_mod_coefs[1:round(nrow(sc_mod_coefs) * 0.15), ]

merged_df <- data.frame(ensembl_gene_id = unique(c(full_mod_coefs$ensembl_gene_id,
                                                   sc_mod_coefs$ensembl_gene_id)))

write.csv(merged_df, file = "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_SCmod_chron_age/fullMod_andSCMod_genes.csv")