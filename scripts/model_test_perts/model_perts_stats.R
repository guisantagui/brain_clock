################################################################################
# Brain clock: test for statistical significance of the differences in         #
# predicted ages of the perturbations in the pert samples. In particular,      #
# xompute t-test and wilcoxon test p values for both predicted chronological   #
# and transformed ages, and adjust p value with BH method.                     #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Train the GLM on transformed age and assess performance.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--respVar",
                               "--metDat"),
                       help = c("Predicted ages dataframe of perturbations (generated with model_test_perts.R).",
                                "Response variable used to fit the model to (age_chron or age_trans).",
                                "Metadata file"),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
pertFile <- parsed$input
respVar <- parsed$respVar
metDatFile <- parsed$metDat

outName <- sprintf("%s/pred_ages_stats.csv",
                   dirname(pertFile))

# Functions
################################################################################

# Read csv faster
readCsvFast <- function(f){
        df <- data.frame(data.table::fread(f))
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

# Load the data
################################################################################

perts <- readCsvFast(pertFile)
metDat <- readCsvFast(metDatFile)

# Keep only perturbation samples
metDat <- metDat[metDat$diagn_4BrainClck == "perturbations", ]

# Test for statistical significance of unique perturbations
################################################################################
print(sprintf("Running univariate tests of %s (%s)...",
              basename(pertFile),
              basename(dirname(pertFile))))

uniqPerts <- unique(metDat$perturbation)

if(respVar == "age_trans"){
        pert_stats_DF <- data.frame(matrix(nrow = 0,
                                   ncol = 8,
                                   dimnames = list(NULL,
                                                   c("perturbation",
                                                     "cell_type",
                                                     "trans_age_wilcox_pVal",
                                                     "trans_age_tTest_pVal",
                                                     "chron_age_wilcox_pVal",
                                                     "chron_age_tTest_pVal",
                                                     "trans_age_logFC",
                                                     "chron_age_logFC"))))

        pb = txtProgressBar(min = 0, max = length(uniqPerts), initial = 0, style = 3)
        for(i in seq_along(uniqPerts)){
                setTxtProgressBar(pb, i)
                p <- uniqPerts[i]
                metDat_p <- metDat[metDat$perturbation == p & metDat$perturbation != "none", ]
                cTypes <- unique(metDat_p$tissue)
                cTypeVec <- c()
                for(j in seq_along(cTypes)){
                        cTyp <- cTypes[j]
                        ctrls_IDs <- metDat$specimenID[metDat$tissue == cTyp & metDat$exper_group == "ctrl"]
                        exper_IDs <- metDat_p$specimenID[metDat_p$tissue == cTyp & metDat_p$exper_group == "expr"]
                        ctrls_IDs <- make.names(ctrls_IDs)
                        exper_IDs <- make.names(exper_IDs)
                        if(length(ctrls_IDs) >= 2 & length(exper_IDs) >= 2){
                                perts_ctrl_transAgeVec <- perts$trans_age[perts$specimenID %in% ctrls_IDs]
                                perts_ctrl_chronAgeVec <- perts$chron_age[perts$specimenID %in% ctrls_IDs]
                                perts_expr_transAgeVec <- perts$trans_age[perts$specimenID %in% exper_IDs]
                                perts_expr_chronAgeVec <- perts$chron_age[perts$specimenID %in% exper_IDs]

                                transAge_wilcox_pVal <- wilcox.test(perts_ctrl_transAgeVec,
                                                                    perts_expr_transAgeVec)$p.value
                                transAge_tTest_pVal <- t.test(perts_ctrl_transAgeVec,
                                                              perts_expr_transAgeVec)$p.value

                                chronAge_wilcox_pVal <- wilcox.test(perts_ctrl_chronAgeVec,
                                                                    perts_expr_chronAgeVec)$p.value
                                chronAge_tTest_pVal <- t.test(perts_ctrl_chronAgeVec,
                                                              perts_expr_chronAgeVec)$p.value

                                transAge_logFC <- log2(median(perts_expr_transAgeVec)/median(perts_ctrl_transAgeVec))
                                chronAge_logFC <- log2(median(perts_expr_chronAgeVec)/median(perts_ctrl_chronAgeVec))

                                toBindDF <- data.frame(perturbation = p,
                                                       cell_type = cTyp,
                                                       trans_age_wilcox_pVal = transAge_wilcox_pVal,
                                                       trans_age_tTest_pVal = transAge_tTest_pVal,
                                                       chron_age_wilcox_pVal = chronAge_wilcox_pVal,
                                                       chron_age_tTest_pVal = chronAge_tTest_pVal,
                                                       trans_age_logFC = transAge_logFC,
                                                       chron_age_logFC = chronAge_logFC)
                                pert_stats_DF <- rbind.data.frame(pert_stats_DF, toBindDF)
                        }
                }
        }
        close(pb)

        # Adjust p-values
        pert_stats_DF$trans_age_wilcox_pValAdj <- p.adjust(pert_stats_DF$trans_age_wilcox_pVal,
                                                           method = "BH")
        pert_stats_DF$trans_age_tTest_pValAdj <- p.adjust(pert_stats_DF$trans_age_tTest_pVal,
                                                          method = "BH")
        pert_stats_DF$chron_age_wilcox_pValAdj <- p.adjust(pert_stats_DF$chron_age_wilcox_pVal,
                                                           method = "BH")
        pert_stats_DF$chron_age_tTest_pValAdj <- p.adjust(pert_stats_DF$chron_age_tTest_pVal,
                                                          method = "BH")
}else if(respVar == "age_chron"){
        pert_stats_DF <- data.frame(matrix(nrow = 0,
                                           ncol = 5,
                                           dimnames = list(NULL,
                                                           c("perturbation",
                                                             "cell_type",
                                                             "chron_age_wilcox_pVal",
                                                             "chron_age_tTest_pVal",
                                                             "chron_age_logFC"))))

        pb = txtProgressBar(min = 0, max = length(uniqPerts), initial = 0, style = 3)
        for(i in seq_along(uniqPerts)){
                setTxtProgressBar(pb, i)
                p <- uniqPerts[i]
                metDat_p <- metDat[metDat$perturbation == p & metDat$perturbation != "none", ]
                cTypes <- unique(metDat_p$tissue)
                cTypeVec <- c()
                for(j in seq_along(cTypes)){
                        cTyp <- cTypes[j]
                        ctrls_IDs <- metDat$specimenID[metDat$tissue == cTyp & metDat$exper_group == "ctrl"]
                        exper_IDs <- metDat_p$specimenID[metDat_p$tissue == cTyp & metDat_p$exper_group == "expr"]
                        ctrls_IDs <- make.names(ctrls_IDs)
                        exper_IDs <- make.names(exper_IDs)
                        if(length(ctrls_IDs) >= 2 & length(exper_IDs) >= 2){
                                perts_ctrl_chronAgeVec <- perts$chron_age[perts$specimenID %in% ctrls_IDs]
                                perts_expr_chronAgeVec <- perts$chron_age[perts$specimenID %in% exper_IDs]

                                chronAge_wilcox_pVal <- wilcox.test(perts_ctrl_chronAgeVec,
                                                                    perts_expr_chronAgeVec)$p.value
                                chronAge_tTest_pVal <- t.test(perts_ctrl_chronAgeVec,
                                                              perts_expr_chronAgeVec)$p.value

                                chronAge_logFC <- log2(median(perts_expr_chronAgeVec)/median(perts_ctrl_chronAgeVec))

                                toBindDF <- data.frame(perturbation = p,
                                                       cell_type = cTyp,
                                                       chron_age_wilcox_pVal = chronAge_wilcox_pVal,
                                                       chron_age_tTest_pVal = chronAge_tTest_pVal,
                                                       chron_age_logFC = chronAge_logFC)
                                pert_stats_DF <- rbind.data.frame(pert_stats_DF, toBindDF)
                        }
                }
        }
        close(pb)

        # Adjust p-values
        pert_stats_DF$chron_age_wilcox_pValAdj <- p.adjust(pert_stats_DF$chron_age_wilcox_pVal,
                                                           method = "BH")
        pert_stats_DF$chron_age_tTest_pValAdj <- p.adjust(pert_stats_DF$chron_age_tTest_pVal,
                                                          method = "BH")
}

# Save result
################################################################################
writeCsvFst(pert_stats_DF, file = outName)

print(sprintf("%s saved at %s", basename(outName), dirname(outName)))