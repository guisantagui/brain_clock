################################################################################
# Brain clock: this script gets as input the lincs expression data (NPCs and   #
# NEUs) and the quantNorm(log2(counts + 1)) integrated dataset                 #
# (Synapse+TBI+111_perFiles) and merges them in a single dataset, which will   #
# be the input of the processing pipeline (sva etc).                           #
################################################################################

lincsFile_NPC <- "/mnt/lscratch/users/gsantamaria/test_large_files/NPC/parsed_mats/lincs_NPC_concat_expMat.csv"
lincsFile_NEU <- "/mnt/lscratch/users/gsantamaria/test_large_files/NEU/parsed_mats/lincs_NEU_concat_expMat.csv"
lincsFile_MIC <- "/mnt/lscratch/users/gsantamaria/test_large_files/MICROGLIA-PSEN1/parsed_mats/lincs_MICROGLIA-PSEN1_concat_expMat.csv"

integFile <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc.csv"

metDatInteg <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC.csv"
metDatLincs_NPC <- "/mnt/lscratch/users/gsantamaria/test_large_files/NPC/parsed_mats/lincs_NPC_concat_metDat.csv"
metDatLincs_NEU <- "/mnt/lscratch/users/gsantamaria/test_large_files/NEU/parsed_mats/lincs_NEU_concat_metDat.csv"
metDatLincs_MIC <- "/mnt/lscratch/users/gsantamaria/test_large_files/MICROGLIA-PSEN1/parsed_mats/lincs_MICROGLIA-PSEN1_concat_metDat.csv"

outDir <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/"
outName <- gsub(".csv", "", basename(integFile))
outPath <- sprintf("%s%s_wLINCS_wLINCS_NPC_NEU_MIC.csv", outDir, outName)

outNameMetDat <- gsub(".csv", "", basename(metDatInteg))
outPathMetDat <- sprintf("%s%s_wLINCS_NPC_NEU_MIC.csv", outDir, outNameMetDat)
# Functions
################################################################################

# read.csv but faster
readCsvFst <- function(pth, header = T){
        df <- data.frame(data.table::fread(pth, header = header))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

# write.csv, but faster
writeCsvFst <- function(df, file, rowNames = T, colNames = T){
        if(rowNames){
                rn <- rownames(df)
                df <- data.table::data.table(df)
                df[, V1 := rn]
                data.table::setcolorder(df, c("V1", setdiff(names(df), "V1")))
        }else{
                df <- data.table::data.table(df)
        }
        data.table::fwrite(df, file, col.names = colNames)
}

# Load data
################################################################################
lincs_NPC <- readCsvFst(lincsFile_NPC)
lincs_NEU <- readCsvFst(lincsFile_NEU)
lincs_MIC <- readCsvFst(lincsFile_MIC)
integ <- readCsvFst(integFile)
lincs_NPC <- data.frame(t(lincs_NPC))
lincs_NEU <- data.frame(t(lincs_NEU))
lincs_MIC <- data.frame(t(lincs_MIC))

metDat_lincs_NPC <- readCsvFst(metDatLincs_NPC)
metDat_lincs_NEU <- readCsvFst(metDatLincs_NEU)
metDat_lincs_MIC <- readCsvFst(metDatLincs_MIC)
metDat_integ <- readCsvFst(metDatInteg)

# Merge LINCS expression and metadata files
################################################################################

lincs <- rbind.data.frame(lincs_NPC, lincs_NEU, lincs_MIC)
metDat_lincs <- rbind.data.frame(metDat_lincs_NPC, metDat_lincs_NEU, metDat_lincs_MIC)

# Merge expression files
################################################################################

commGenes <- intersect(colnames(lincs), colnames(integ))

lincsInteg <- rbind.data.frame(lincs[, commGenes],
                               integ[, commGenes])

# Merge metadata files
################################################################################
metDat_lincs$diagn_4BrainClck <- rep("perturbations",
                                     nrow(metDat_lincs))

metDat_lincs$perturbation <- paste(metDat_lincs$pertname,
                                   metDat_lincs$dose,
                                   metDat_lincs$timepoint,
                                   sep = "_")

metDat_lincs$substudy <- rep("LINCS", nrow(metDat_lincs))
metDat_lincs$RIN <- rep(max(metDat_integ$RIN[!is.na(metDat_integ$RIN)]),
                        nrow(metDat_lincs))

metDat_lincs <- metDat_lincs[, c("specimenID", "diagn_4BrainClck",
                                 "perturbation", "group", "RIN", "substudy", "cell")]
colnames(metDat_lincs) <- gsub("group", "exper_group", colnames(metDat_lincs))
colnames(metDat_lincs) <- gsub("cell", "tissue", colnames(metDat_lincs))

cols2Add_metDatLincs <- colnames(metDat_integ)[!colnames(metDat_integ) %in% colnames(metDat_lincs)]

df2AddLincs <- data.frame(matrix(nrow = nrow(metDat_lincs),
                                 ncol = length(cols2Add_metDatLincs),
                                 dimnames = list(rownames(metDat_lincs),
                                                 cols2Add_metDatLincs)))

metDat_lincs <- cbind.data.frame(metDat_lincs, df2AddLincs)

metDat_lincsInteg <- rbind.data.frame(metDat_integ,
                                      metDat_lincs[, colnames(metDat_integ)])

medianAge <- median(metDat_lincsInteg$ageDeath[!is.na(metDat_lincsInteg$ageDeath)])
metDat_lincsInteg$ageDeath[is.na(metDat_lincsInteg$ageDeath)] <- medianAge
# Save expression and metadata datasets
################################################################################
writeCsvFst(lincsInteg, file = outPath)
print(sprintf("%s saved at %s.",
              basename(outPath),
              dirname(outPath)))

writeCsvFst(metDat_lincsInteg, file = outPathMetDat)
print(sprintf("%s saved at %s.",
              basename(outPathMetDat),
              dirname(outPathMetDat)))