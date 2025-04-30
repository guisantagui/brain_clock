################################################################################
# Brain clock: Parse perturbation datasets curated from the literature         #
# (described in table S5).                                                     #
################################################################################

if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}

if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
library(plotUtils)

# Directory stuff
################################################################################
setwd("../..")

counts_111_sas_f <- "./data/perturbation/111_sascha_complete_counts_df.csv"
counts_111_jav_f <- "./data/perturbation/javier_111_counts.csv"

metDatSasFile <- "./data/perturbation/Sascha_metadata.csv"
metDatJavFile <- "./data/perturbation/javier_metadata.csv"

out_dir <- "./results/parsed/"
create_dir_if_not(out_dir)

# Load data
################################################################################
counts_111_sas <- read_table_fast(counts_111_sas_f, row.names = 1)
counts_111_jav <- read_table_fast(counts_111_jav_f, row.names = 1)
metDatSas <- read_table_fast(metDatSasFile)
metDatJav <- read_table_fast(metDatJavFile)

# Parsing
################################################################################

# Merge counts files
common_genes <- intersect(rownames(counts_111_sas),
                          rownames(counts_111_jav))

# Merge metadata files
metDatSas$Intervention <- paste(metDatSas$Gene,
                                metDatSas$Method,
                                sep = "_")

colnames(metDatSas) <- gsub("Biosample.Name", "Celltype", colnames(metDatSas))

metDat_111 <- rbind.data.frame(metDatSas[, c("GSE.Series", "Celltype",
                                             "Intervention", "Experimental_GSM",
                                             "Control_GSM")],
                               metDatJav[, c("GSE.Series", "Celltype",
                                             "Intervention", "Experimental_GSM",
                                             "Control_GSM")])

# Filter cell types to keep only the ones we're interested in

# Cell types to be kept
cellTypes2KeepJav <- c("astrocyte",
                       "Cultured cortical neurons,",
                       "GABAergic neurons",
                       "hCMEC/D3",
                       "human brain",
                       "human neurons",
                       "monocyte-derived dendritic cells",
                       "monocytes-derived dendritic cells ",
                       "NSC ")

cellTypes2KeepSas <- c("human motor neurons",
                       "Neural precursor cells",
                       "neuroepithelial stem cells",
                       "Neuron",
                       "NHA",
                       "NPC",
                       "NPCs",
                       "NSC",
                       "NSCs",
                       "U5-NPC (G1 Phase)")

cellTypesKeep <- c(cellTypes2KeepJav,
                   cellTypes2KeepSas)
metDat_111 <- metDat_111[metDat_111$Celltype %in% cellTypesKeep, ]

# Iterate over GSE.Series to select each GSM from the table, in case it exists,
# and create new metadata file that has info for each gsm about gseSeries, cell
# type, perturbation and experimental group (if it's control or experimental).
# When a unique GSM maps to more than one column in counts_111 compute average,
# as it is due to different runs of same sample.
counts_111_filt <- data.frame(matrix(nrow = nrow(counts_111),
                                     ncol = 0,
                                     dimnames = list(rownames(counts_111), NULL)))

metDat_111_parsed <- data.frame(matrix(nrow = 0, ncol = 5,
                                       dimnames = list(NULL,
                                                       c("specimenID",
                                                         "GSE_series",
                                                         "cell_type",
                                                         "perturbation",
                                                         "exper_group"))))

for(i in seq_along(metDat_111$GSE.Series)){
        gse <- metDat_111$GSE.Series[i]
        xprs <- strsplit(metDat_111$Experimental_GSM[i], split = ";")[[1]]
        ctls <- strsplit(metDat_111$Control_GSM[i], split = ";")[[1]]
        cellType <- metDat_111$Celltype[i]
        pert <- metDat_111$Intervention[i]
        gsms_kept <- c()
        gsms <- c(xprs, ctls)
        expGroup <- c()
        for(j in seq_along(gsms)){
                gsm <- gsms[j]
                idx <- grep(gsm, colnames(counts_111))
                if(length(idx) > 0){
                        counts_111_gsm <- counts_111[, idx]
                        if(!is.null(dim(counts_111_gsm))){
                                counts_111_gsm <- apply(counts_111_gsm, 1, mean)
                        }
                        counts_111_gsm <- matrix(counts_111_gsm,
                                                 nrow = length(counts_111_gsm),
                                                 ncol = 1,
                                                 dimnames = list(names(counts_111_gsm),
                                                                 gsm))
                        counts_111_filt <- cbind.data.frame(counts_111_filt, counts_111_gsm)
                        gsms_kept <- c(gsms_kept, gsm)
                        if(gsm %in% xprs){
                                group <- "expr"
                        }else if(gsm %in% ctls){
                                group <- "ctrl"
                        }
                        expGroup <- c(expGroup, group)
                }
        }
        toBindMetDF <- data.frame(specimenID = gsms_kept,
                                  GSE_series = rep(gse, length(gsms_kept)),
                                  cell_type = rep(cellType, length(gsms_kept)),
                                  perturbation = rep(pert, length(gsms_kept)),
                                  exper_group = expGroup)
        metDat_111_parsed <- rbind.data.frame(metDat_111_parsed,
                                              toBindMetDF)
}

# We are dealing with some duplicates due to same GSMs appearing in 
# several rows of the metadata. In Counts are identical columns, so let's
# just get rid of them

metDat_111_parsed$perturbation[metDat_111_parsed$exper_group == "ctrl"] <- "none"
metDat_111_parsed <- metDat_111_parsed[!duplicated(metDat_111_parsed$specimenID), ]
counts_111_filt <- counts_111_filt[, !duplicated(colnames(counts_111_filt))]

metDat_111_parsed$substudy <- rep("perts_111",
                                  nrow(metDat_111_parsed))
metDat_111_parsed$assay <- rep("rnaSeq",
                               nrow(metDat_111_parsed))
metDat_111_parsed$diagn_4BrainClck <- rep("perturbations",
                                          nrow(metDat_111_parsed))
# Change GSE_series to batch_rna, as GSE series correspond to experiments
colnames(metDat_111_parsed) <- gsub("GSE_series", "batch_rna", colnames(metDat_111_parsed))
# Chance cel_type for tissue, for later unifying the category as tissue_cellType
# (once merged)
colnames(metDat_111_parsed) <- gsub("cell_type", "tissue", colnames(metDat_111_parsed))
metDat_111_parsed$RIN <- 10

# Remove version from gene names
rownames(counts_111_filt) <- gsub("\\..*", "", rownames(counts_111_filt))

write_table_fast(counts_111_filt,
                 f = sprintf("%scounts_111.csv", out_dir))
write_table_fast(metDat_111_parsed,
                 f = sprintf("%smetDat_111.csv", out_dir))
