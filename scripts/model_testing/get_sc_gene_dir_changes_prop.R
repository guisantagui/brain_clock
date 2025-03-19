if(!require("ggplot2", quietly = T)){
        install.packages("ggplot2", repos='http://cran.us.r-project.org')
}
library(ggplot2)
# Functions
################################################################################

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Read csv faster
readCsvFast <- function(f){
        df <- data.frame(data.table::fread(f))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

# Directory stuff
################################################################################
datFile <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_NPC_NEU_MIC.csv"
datFile <- "/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast_allLINCSBrain_filtSignChron/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_NPC_NEU_MIC_modFuncsAlpha1_coefs_noCerebell_onlyAge_svaAdj.csv"
metDatFile <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS_NPC_NEU_MIC.csv"
coefsFile <- "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_chron_age_onSignGenes/modFuncsAlpha1_coefs.csv"
outDir <- "/home/users/gsantamaria/projects/brain_clock/results/model_test_sCell/modAllGenes_integWAllLincs_and_sc_oAge_chronAge_alph1_onSignGenes"
youngThrshld <- 30
oldThrshld <- 70

outDir <- addSlashIfNot(outDir)
createIfNot(outDir)

# Load data
################################################################################
dat <- readCsvFast(datFile)
metDat <- readCsvFast(metDatFile)
coefs <- readCsvFast(coefsFile)

# Keep only ageAnno rows and coefs columns.
metDat <- metDat[metDat$substudy == "ageAnno", ]
ageAnnoSamps <- make.names(metDat$specimenID)
dat_all <- dat[rownames(dat) %in% ageAnnoSamps, ]
dat <- dat[rownames(dat) %in% ageAnnoSamps, coefs$ensembl_gene_id]

write.csv(dat, file = sprintf("%spseuBulkCounts_genesInMod.csv", outDir))

# Divide dat in young and old
young_samps <- metDat$specimenID[metDat$ageDeath <= youngThrshld]
old_samps <- metDat$specimenID[metDat$ageDeath >= oldThrshld]

dat_young <- dat[rownames(dat) %in% make.names(young_samps), ]
dat_old <- dat[rownames(dat) %in% make.names(old_samps), ]

cell_types <- unique(metDat$tissue)

sign_match_mat <- data.frame(matrix(nrow = 0, ncol = ncol(dat), dimnames = list(NULL, colnames(dat))))
pVal_binom_vec <- c()
n_perm <- 10000
pVal_perm_vec <- c()
for(i in seq_along(cell_types)){
        cell <- cell_types[i]
        print(sprintf("Computing: %s...", cell))
        dat_young_cell <- dat_young[grepl(cell, rownames(dat_young)), ]
        dat_old_cell <- dat_old[grepl(cell, rownames(dat_old)), ]
        cell_sign_match <- c()
        for(j in 1:ncol(dat)){
                gene <- colnames(dat)[j]
                gene_diff <- median(dat_old_cell[, gene]) - median(dat_young_cell[, gene])
                gene_coef <- coefs$coefficients[coefs$ensembl_gene_id == gene]
                gene_match <- 0
                if (gene_diff < 0 & gene_coef < 0 | gene_diff > 0 & gene_coef > 0 ){
                        gene_match = 1
                }
                cell_sign_match <- c(cell_sign_match, gene_match)
        }
        perm_match_prop_cell_vec <- c()
        for (p in 1:n_perm){
                dat_perm <- dat_all[, sample(colnames(dat_all), size = ncol(dat))]
                colnames(dat_perm) <- colnames(dat)
                dat_cell_perm_young <- dat_perm[young_samps, ]
                dat_cell_perm_old <- dat_perm[old_samps, ]
                dat_cell_perm_young <- dat_cell_perm_young[grepl(cell, rownames(dat_cell_perm_young)), ]
                dat_cell_perm_old <- dat_cell_perm_old[grepl(cell, rownames(dat_cell_perm_old)), ]
                cell_sign_match_perm <- c()
                for(j in 1:ncol(dat_perm)){
                        gene <- colnames(dat_perm)[j]
                        gene_diff <- median(dat_cell_perm_old[, gene]) - median(dat_cell_perm_young[, gene])
                        gene_coef <- coefs$coefficients[coefs$ensembl_gene_id == gene]
                        gene_match <- 0
                        if (gene_diff < 0 & gene_coef < 0 | gene_diff > 0 & gene_coef > 0 ){
                                gene_match = 1
                        }
                        cell_sign_match_perm <- c(cell_sign_match_perm, gene_match)
                }
                perm_match_prop_cell <- sum(cell_sign_match_perm)/length(cell_sign_match_perm)
                perm_match_prop_cell_vec <- c(perm_match_prop_cell_vec,
                                              perm_match_prop_cell)
        }
        real_cell_sign_match_prop <- sum(cell_sign_match)/length(cell_sign_match)
        pVal_perm <- sum(perm_match_prop_cell_vec >= real_cell_sign_match_prop)/n_perm
        cell_sign_match <- data.frame(matrix(cell_sign_match,
                                             nrow = 1,
                                             ncol = ncol(dat),
                                             dimnames = list(cell,
                                                             colnames(dat))))
        pVal_binom <- binom.test(x = rowSums(cell_sign_match),
                                 n = ncol(cell_sign_match),
                                 alternative = "greater")$p.value
        sign_match_mat <- rbind.data.frame(sign_match_mat, cell_sign_match)
        pVal_binom_vec <- c(pVal_binom_vec, pVal_binom)
        pVal_perm_vec <- c(pVal_perm_vec, pVal_perm)
}

sign_match_mat$prop_same_dir <- apply(sign_match_mat,
                                      1,
                                      function(x) sum(x)/length(x))

sign_match_mat$pVal_binom <- pVal_binom_vec
sign_match_mat$pVal_perm <- pVal_perm_vec

write.csv(sign_match_mat, file = paste0(outDir, "same_dir_sign_genes_prop.csv"))

# Do a plot
df_4Plot <- data.frame(matrix(nrow = 0, ncol = 3,
                              dimnames = list(NULL,
                                              c("cell_type",
                                                "gene",
                                                "is_correlated"))))
for(i in seq_along(cell_types)){
        cell <- cell_types[i]
        cell_sign_vec <- sign_match_mat[cell,
                                        !colnames(sign_match_mat) %in% c("prop_same_dir",
                                                                         "pVal_binom")]
        cell_sign_vec <- unlist(as.vector(cell_sign_vec))
        toBind <- data.frame(cell_type = rep(cell, length(cell_sign_vec)),
                             gene = names(cell_sign_vec),
                             is_correlated = cell_sign_vec)
        df_4Plot <- rbind.data.frame(df_4Plot, toBind)
}
df_4Plot$cell_type <- gsub(".", " ", df_4Plot$cell_type, fixed = T)
df_4Plot$is_correlated <- factor(df_4Plot$is_correlated)
plt <- ggplot(df_4Plot, mapping = aes(x = cell_type, fill = is_correlated)) +
        geom_bar(position = "fill") +
        scale_fill_manual(values = c("red", "chartreuse3"), labels = c("No correlation", "Correlation")) +
        ylab("proportion") +
        xlab("cell type") +
        theme(title = ggtext::element_markdown(),
                      axis.title.x = ggtext::element_markdown(),
                      axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, size = 13),
                      axis.title.y = ggtext::element_markdown(),
                      legend.title = element_blank(),
                      legend.text = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
ggsave(filename = sprintf("%sSC_sameSign_prop.pdf", outDir), plot = plt, height = 10, width = 4)
