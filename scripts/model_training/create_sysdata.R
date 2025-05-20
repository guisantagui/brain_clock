######################################################################
# Brain clock: create sysdata.rda object that is used by             #
# brainAgeShiftR.                                                    #
######################################################################


if (!require("BiocManager", quietly = T)){
        install.packages("BiocManager")
}
if (!require("biomaRt", quietly = T)){
        BiocManager::install("biomaRt", update = F)
}
if (!require(plotUtils)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
library(biomaRt)
library(plotUtils)

# Files
######################################################################

# Model coefficient file
mod_coefs_f <- "../../results/models/secnd_round/mod_alpha1_coefs_wIntercpt.csv"

# First round coefficient file, to get the background genes for
# normalization
coefs_first_f <- "../../results/models/first_round/mod_alpha1_coefs.csv"

# Training data reference distribution

train_target_f <- "../../results/parsed/merged/train_target.rds"

outDir <- dirname(mod_coefs_f)

# Load data
######################################################################

# Load model coefficients
mod_coef <- read.csv(mod_coefs_f, row.names = 1)
colnames(mod_coef) <- gsub("ensembl_gene_id", "names", colnames(mod_coef))
row.names(mod_coef) <- 1:nrow(mod_coef)

# Load clinical counts
coefs_first <- read_table_fast(coefs_first_f, row.names = 1)

# Create gene info dataframe
######################################################################

human <- useMart("ensembl", dataset="hsapiens_gene_ensembl",
                 host = "https://www.ensembl.org")

mod_gene_info <- getBM(attributes = c("ensembl_gene_id",
                                      "hgnc_symbol",
                                      "entrezgene_id",
                                      "description"),
                   filters = "ensembl_gene_id",
                   values = mod_coef$name[mod_coef$name != "Intercept"],
                   mart = human)

# This has 2 entrez ids, so remove the second one.
mod_gene_info <- gene_info[gene_info$entrezgene_id != 107080644, ]

mod_gene_info <- mod_gene_info[match(mod_coef$names[2:nrow(mod_coef)],
                                     mod_gene_info$ensembl_gene_id), ]
rownames(mod_gene_info) <- 1:nrow(mod_gene_info)

# Create bg_genes object
######################################################################
bg_genes <- coefs_first$ensembl_gene_id

# Load train_quant_means object
######################################################################

train_quant_means <- readRDS(train_target_f)


# Save the sysdata.rda file
######################################################################

save(mod_coef,
     mod_gene_info,
     bg_genes,
     train_quant_means,
     file = sprintf("%s/sysdata.rda",
                    outDir))