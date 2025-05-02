################################################################################
# Brain clock: intersects merged dataset with LINCS genes, remove lowly        #
# expressed genes and samples with RIN lower than 6, log2 transform and        #
# quantile-normalize.                                                          #
################################################################################
if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}
if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
library(plotUtils)

# Functions
################################################################################
# Finds genes that are expressed at less than an average of 1 per sample 
# in each one of the ages
findLowExpGenes_perAge <- function(df, metDat){
        lowExpGenes_perAge <- list()
        ageVec <- unique(metDat$ageDeath)
        ageVec <- ageVec[!is.na(ageVec)]
        for(age in ageVec){
                sampsAge <- metDat$specimenID[metDat$ageDeath == age]
                sampsAge <- sampsAge[!is.na(sampsAge)]
                sampsAge <- make.names(sampsAge)
                sampsAge <- sampsAge[sampsAge %in% rownames(df)]
                df_age <- df[make.names(sampsAge), ]
                if(length(sampsAge) > 1){
                        lowExpGenes <- colnames(df_age)[colSums(df_age) < nrow(df_age)]
                        lowExpGenes_perAge[[paste("age",
                                                  as.character(age),
                                                  sep = "_")]] <- lowExpGenes
                }
        }
        lowExpGenes_perAge <- lowExpGenes_perAge[lapply(lowExpGenes_perAge,
                                                        length) != 0]
        return(lowExpGenes_perAge)
}

findLowExpGenes_perAgeBin <- function(df, metDat, bin_width = 5){
    lowExpGenes_perBin <- list()
    
    # Bin ages into 5-year intervals (e.g., 60–64, 65–69, etc.)
    metDat$ageBin <- cut(metDat$ageDeath,
                         breaks = seq(floor(min(metDat$ageDeath, na.rm = TRUE)),
                                      ceiling(max(metDat$ageDeath, na.rm = TRUE)),
                                      by = bin_width),
                         include.lowest = TRUE, right = FALSE)

    ageBins <- unique(metDat$ageBin)
    ageBins <- ageBins[!is.na(ageBins)]

    for (bin in ageBins) {
        sampsBin <- metDat$specimenID[metDat$ageBin == bin]
        sampsBin <- sampsBin[!is.na(sampsBin)]
        sampsBin <- make.names(sampsBin)
        sampsBin <- sampsBin[sampsBin %in% rownames(df)]
        df_bin <- df[make.names(sampsBin), ]
        
        if (length(sampsBin) > 1) {
            lowExpGenes <- colnames(df_bin)[colSums(df_bin) < nrow(df_bin)]
            lowExpGenes_perBin[[as.character(bin)]] <- lowExpGenes
        }
    }

    # Remove empty entries
    lowExpGenes_perBin <- lowExpGenes_perBin[lapply(lowExpGenes_perBin, length) != 0]
    return(lowExpGenes_perBin)
}

# Look for the intersection across all elements in a list.
recurs_intersect <- function(lst){
        if (length(lst) == 1) {
                return(lst[[1]])
        }
        
        # Recursive case: find the intersection of the first element with the
        # intersection of the rest of the list
        return(intersect(lst[[1]], recurs_intersect(lst[-1])))
}

# Directory stuff
################################################################################
counts_f <- "../../results/parsed/merged/merged_counts.csv"
lincsGeneInfo <- "../../data/perturbation/lincs/lincs_genes.csv"
metdat_f <- "../../results/parsed/merged/merged_metdat.csv"

propZeroRem <- .8
rinFilt <- 6
outDir <- "../../results/parsed/merged/"

create_dir_if_not(outDir)

outName <- sprintf("%s%s_log2_quantNorm.csv", outDir,
                   gsub(".csv", "", basename(counts_f)))

# Load data
################################################################################

counts <- read_table_fast(counts_f, row.names = 1)
lincs_gi <- read.csv(lincsGeneInfo, row.names = 1)
metdat <- read_table_fast(metdat_f, row.names = 1)

# Preprocessing
################################################################################

# Keep only genes in LINCS
commGenes <- intersect(lincs_gi$gene, colnames(counts))
counts <- counts[, commGenes]
# Ensure that only what is available in the metadata is in counts
counts <- counts[make.names(rownames(counts)) %in% make.names(metdat$specimenID), ]
# Remove genes with a high proportion of zeros:
counts <- counts[, apply(counts,
                         2, function(x) sum(x == 0)/length(x)) < propZeroRem]

# Find genes that at each 5 year interval age, express less than an average of one count per
# gene
lowExp_genes_perAgeBin <- findLowExpGenes_perAgeBin(counts, metdat)
to_rem <- recurs_intersect(lowExp_genes_perAgeBin)
counts <- counts[, !colnames(counts) %in% to_rem]

counts <- counts[make.names(rownames(counts)) %in% make.names(metdat$specimenID[!is.na(metdat$RIN) & metdat$RIN >= rinFilt]), ]

counts <- t(counts)

# Transform
counts_log2 <- log2(counts + 1)


quant_norm <- function(m, axis = 2, train_means = NULL) {
    m_class <- class(m)[1]
    if ((!m_class %in% c("data.frame", "matrix")) | any(dim(m) == 0)){
        stop("Invalid input object.", call. = FALSE)
    }
    if (!axis %in% c(1, 2)){
        stop("Invalid axis.", call. = FALSE)
    }
    if (axis == 1) {
        m <- t(m)
    }

    # Sort the matrix
    m_sort <- apply(m, 2, sort)

    # Calculate row means
    if (!is.null(train_means)) {
        means <- train_means
    } else {
        means <- rowMeans(m_sort)
    }

    # Normalize
    m_norm <- matrix(0, nrow = nrow(m), ncol = ncol(m))
    for (i in 1:ncol(m)) {
        m_norm[, i] <- means[rank(m[, i], ties.method = "average")]
    }

    if (axis == 1) {
        m_norm <- t(m_norm)
    }

    dimnames(m_norm) <- dimnames(m)
    if (m_class == "data.frame") {
        m_norm <- as.data.frame(m_norm)
    }

    return(list(norm = m_norm, ref = means))
}

counts_log2_qNorm <- quant_norm(counts_log2)

means_ref <- counts_log2_qNorm$ref
counts_log2_qNorm <- counts_log2_qNorm$norm

counts_log2_qNorm <- t(counts_log2_qNorm)

# Save
write_table_fast(counts_log2_qNorm, outName)

outName_ref <- gsub(".csv", "_refMeans_4norm.rds", outName)
saveRDS(means_ref, file = outName_ref)