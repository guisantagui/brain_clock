################################################################################
# Brain clock: this script gets as input the lincs expression data (NPCs and   #
# NEUs) and the compilation of perturbation studies datasets, merges them and  #
# quantile-normalizes the merged set using the ranked means obtained during    #
# quantile-normalization of the combined clinical dataset (ROSMAP etc.) as     #
# reference.                                                                   #
# It takes as input the human merged metadata file too, and uses it to reshape #
# integrated perturbation metadata to have the same columns as the human one.  #
################################################################################

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}
if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
if (!require("preprocessCore", quietly = T)){
        BiocManager::install("preprocessCore", configure.args=c(preprocessCore = "--disable-threading"),
                             update = F)
}
library(preprocessCore)
library(plotUtils)

# Files
################################################################################

lincsFile_NPC <- "/mnt/lscratch/users/gsantamaria/test_large_files/NPC/parsed_mats/lincs_NPC_concat_expMat.csv"
lincsFile_NEU <- "/mnt/lscratch/users/gsantamaria/test_large_files/NEU/parsed_mats/lincs_NEU_concat_expMat.csv"
lincsFile_MIC <- "/mnt/lscratch/users/gsantamaria/test_large_files/MICROGLIA-PSEN1/parsed_mats/lincs_MICROGLIA-PSEN1_concat_expMat.csv"

pertsFile <- "/home/users/gsantamaria/projects/brain_clock/results/parsed/counts_111.csv"

metDatPerts_f <- "/home/users/gsantamaria/projects/brain_clock/results/parsed/metDat_111.csv"
metDatLincs_NPC_f <- "/mnt/lscratch/users/gsantamaria/test_large_files/NPC/parsed_mats/lincs_NPC_concat_metDat.csv"
metDatLincs_NEU_f <- "/mnt/lscratch/users/gsantamaria/test_large_files/NEU/parsed_mats/lincs_NEU_concat_metDat.csv"
metDatLincs_MIC_f <- "/mnt/lscratch/users/gsantamaria/test_large_files/MICROGLIA-PSEN1/parsed_mats/lincs_MICROGLIA-PSEN1_concat_metDat.csv"

#counts_clin_f <- "/home/users/gsantamaria/projects/brain_clock/results/parsed/merged/merged_counts_log2_qnorm.csv"
counts_clin_f <- "/home/users/gsantamaria/projects/brain_clock/results/preproc/test_no_lincs/merged_counts_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
counts_clin_f <- "/home/users/gsantamaria/projects/brain_clock/results/preproc/second_round/merged_counts_mod_alpha1_coefs_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
metDat_clin_f <- "/home/users/gsantamaria/projects/brain_clock/results/parsed/merged/merged_metdat.csv"
train_test_f <- "/home/users/gsantamaria/projects/brain_clock/results/parsed/merged/train_test.csv"

outDir <- "/home/users/gsantamaria/projects/brain_clock/results/parsed/merged_perts/"
create_dir_if_not(outDir)

outPath <- sprintf("%smerged_perts_exprsn.csv", outDir)
lincsOut <- sprintf("%slincs_exprsn.csv", outDir)

outPathMetDat <- sprintf("%smerged_perts_metdat.csv", outDir)

# Load data
################################################################################
lincs_NPC <- read_table_fast(lincsFile_NPC, row.names = 1)
lincs_NEU <- read_table_fast(lincsFile_NEU, row.names = 1)
lincs_MIC <- read_table_fast(lincsFile_MIC, row.names = 1)
perts <- read_table_fast(pertsFile, row.names = 1)
lincs_NPC <- data.frame(t(lincs_NPC))
lincs_NEU <- data.frame(t(lincs_NEU))
lincs_MIC <- data.frame(t(lincs_MIC))
perts <- data.frame(t(perts))
counts_clin <- read_table_fast(counts_clin_f, row.names = 1)

metDat_lincs_NPC <- read_table_fast(metDatLincs_NPC_f, row.names = 1)
metDat_lincs_NEU <- read_table_fast(metDatLincs_NEU_f, row.names = 1)
metDat_lincs_MIC <- read_table_fast(metDatLincs_MIC_f, row.names = 1)
metDat_perts <- read_table_fast(metDatPerts_f, row.names = 1)
metDat_clin <- read_table_fast(metDat_clin_f, row.names = 1)

# Merge LINCS expression and metadata files
################################################################################

lincs <- rbind.data.frame(lincs_NPC,
                          #lincs_MIC,
                          lincs_NEU)
metDat_lincs <- rbind.data.frame(metDat_lincs_NPC,
                                 #metDat_lincs_MIC,
                                 metDat_lincs_NEU)

# Test: try to include each NPC subtype
#metDat_lincs$cell <- sapply(metDat_lincs$rna_plate,
#                            function(x) strsplit(x, split = "_")[[1]][2])

#lincs <- rbind.data.frame(lincs_NPC,
#                          lincs_NEU,
#                          lincs_MIC)
#metDat_lincs <- rbind.data.frame(metDat_lincs_NPC,
#                                 metDat_lincs_NEU,
#                                 metDat_lincs_MIC)

# There are several samples where the cell type in the metadata of the GCTX file is
# NPC but the ID indicates subtypes of NPCs such as FIBRNPC, NPC.170 etc. Remove then
# all that in the id is not labeled as exactly NPC 
#metDat_lincs <- metDat_lincs[!(metDat_lincs$cell == "NPC" & sapply(metDat_lincs$id,
#                                                                   function(x) strsplit(x, split = "_")[[1]][2]) != "NPC"), ]

# There are other types of controls than DMSO. Keep only DMSO
metDat_lincs <- metDat_lincs[!(metDat_lincs$group == "ctrl" & !grepl("DMSO", metDat_lincs$pertname)), ]

# Keep in expression only what was retained in the metadata
lincs <- lincs[make.names(rownames(lincs)) %in% make.names(metDat_lincs$specimenID), ]

# Merge expression files
################################################################################

commGenes <- intersect(colnames(lincs), colnames(perts))

lincsPerts <- rbind.data.frame(lincs[, commGenes],
                               perts[, commGenes])

# Save common genes as a csv file
commGenes_df <- data.frame(gene = commGenes)
write_table_fast(commGenes_df, f = sprintf("%slincs_genes.csv", outDir))

# Merge metadata files
################################################################################
metDat_lincs$diagn_4BrainClck <- rep("perturbations",
                                     nrow(metDat_lincs))

metDat_lincs$perturbation <- paste(metDat_lincs$pertname,
                                   metDat_lincs$dose,
                                   metDat_lincs$timepoint,
                                   sep = "_")

metDat_lincs$substudy <- rep("LINCS", nrow(metDat_lincs))
#metDat_lincs$RIN <- rep(max(metDat_perts$RIN[!is.na(metDat_perts$RIN)]),
#                        nrow(metDat_lincs))

metDat_lincs <- metDat_lincs[, c("specimenID",
                                 "diagn_4BrainClck",
                                 "perturbation",
                                 "timepoint",
                                 "group",
                                 #"RIN",
                                 "substudy",
                                 "cell",
                                 "rna_plate")]
colnames(metDat_lincs) <- gsub("group", "exper_group", colnames(metDat_lincs))
colnames(metDat_lincs) <- gsub("cell", "tissue", colnames(metDat_lincs))

# Change name of rna_plate to rna_batch to match perturbation dataset
colnames(metDat_lincs) <- gsub("rna_plate", "batch_rna", colnames(metDat_lincs))

# Add columns that are in clinical data as empty columns in perturbation
# datasets
cols2Add_metDatLincs <- colnames(metDat_clin)[!colnames(metDat_clin) %in% colnames(metDat_lincs)]
df2AddLincs <- data.frame(matrix(nrow = nrow(metDat_lincs),
                                 ncol = length(cols2Add_metDatLincs),
                                 dimnames = list(rownames(metDat_lincs),
                                                 cols2Add_metDatLincs)))
metDat_lincs <- cbind.data.frame(metDat_lincs, df2AddLincs)

cols2Add_metDatPerts <- colnames(metDat_clin)[!colnames(metDat_clin) %in% colnames(metDat_perts)]
df2AddPerts <- data.frame(matrix(nrow = nrow(metDat_perts),
                                 ncol = length(cols2Add_metDatPerts),
                                 dimnames = list(rownames(metDat_perts),
                                                 cols2Add_metDatPerts)))
metDat_perts <- cbind.data.frame(metDat_perts, df2AddPerts)

# Add batch info for lincs metadata and copy batch_rna in batch_lib for
# pert111 dataset
metDat_lincs$batch_lib <- gsub("\\_.*", "", metDat_lincs$specimenID)
metDat_perts$batch_lib <- metDat_perts$batch_rna

# Homogenize cell type names in pert111
table(metDat_perts$tissue)
metDat_perts$tissue <- gsub("Cultured cortical neurons,", "NEU", metDat_perts$tissue)
metDat_perts$tissue <- gsub("GABAergic neurons", "NEU", metDat_perts$tissue)
metDat_perts$tissue <- gsub("human motor neurons", "NEU", metDat_perts$tissue)
metDat_perts$tissue <- gsub("human neurons", "NEU", metDat_perts$tissue)
metDat_perts$tissue <- gsub("Neuron", "NEU", metDat_perts$tissue)
metDat_perts$tissue <- gsub("NPCs", "NPC", metDat_perts$tissue)
metDat_perts$tissue <- gsub("NSCs", "NSC", metDat_perts$tissue)
metDat_perts$tissue <- gsub("U5-NPC (G1 Phase)", "NPC", metDat_perts$tissue,
                            fixed = T)

# Add empty column for timepoint for the pert111
metDat_perts$timepoint <- NA
# Get the intersection of the column names
metDat_commCols <- intersect(colnames(metDat_lincs), colnames(metDat_perts))

# Merge metadata dataframes
metDat_lincsPert <- rbind.data.frame(metDat_perts[, metDat_commCols],
                                     metDat_lincs[, metDat_commCols])

# Keep only genes at the intersection with the clinical data
commGenes_clin <- intersect(colnames(lincsPerts), colnames(counts_clin))
lincsPerts <- lincsPerts[, commGenes_clin]
# Save expression and metadata datasets
################################################################################
write_table_fast(lincsPerts, f = outPath)
print(sprintf("%s saved at %s.",
              basename(outPath),
              dirname(outPath)))

#write_table_fast(lincs, f = lincsOut)
#print(sprintf("%s saved at %s.",
#              basename(lincsOut),
#              dirname(lincsOut)))

write_table_fast(metDat_lincsPert, f = outPathMetDat)
print(sprintf("%s saved at %s.",
              basename(outPathMetDat),
              dirname(outPathMetDat)))

# Quantile-normalized merge dataset, using clinical data as reference
################################################################################

# Obtain target distribution from clinical data
Sys.setenv(OMP_NUM_THREADS = "1")

if (!is.null(train_test_f) && file.exists(train_test_f)){
        print(sprintf("Using the training reference extracted from the training samples in %s as indicated in %s...", counts_clin_f, train_test_f))
        train_test <- read_table_fast(train_test_f, row.names = 1)
        train_samps <- make.names(train_test$specimenID[train_test$train_test_split == "train"])
        counts_clin_train <- counts_clin[make.names(rownames(counts_clin)) %in% train_samps, ]
        target_dist <- normalize.quantiles.determine.target(t(counts_clin_train))
}else{
        print(sprintf("Using the training reference extracted from all the samples in %s...", counts_clin_f))
        target_dist <- normalize.quantiles.determine.target(t(counts_clin))
}


lincsPerts <- t(lincsPerts)
lincsPerts_quantNorm <- normalize.quantiles.use.target(as.matrix(lincsPerts),
                                                       target_dist)
# Reassign the dimnames, as normalize.quantiles removes them
dimnames(lincsPerts_quantNorm) <- dimnames(lincsPerts)
lincsPerts_quantNorm <- t(lincsPerts_quantNorm)
outPath_quantNorm <- gsub(".csv", "_qnorm.csv", outPath)
write_table_fast(lincsPerts_quantNorm, f = outPath_quantNorm)
print(sprintf("%s saved at %s.",
              basename(outPath_quantNorm),
              dirname(outPath_quantNorm)))