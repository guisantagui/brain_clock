################################################################################
# Brain clock: test for statistical significance of the differences in         #
# predicted ages of the perturbations in the pert samples. In particular,      #
# compute t-test and wilcoxon test p values for both predicted chronological   #
# and transformed ages, and adjust p value with BH method.                     #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
if(!require(devtools)){
        install.packages("devtools", repos='http://cran.us.r-project.org')
}
if (!require(ggpubr, quietly = T)){
        devtools::install_github("kassambara/ggpubr", upgrade = "never")
}
if (!require(plotUtils, quietly = T)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
if (!require(dplyr, quietly = T)){
        install.packages("dplyr", repos='http://cran.us.r-project.org')
}
if (!require(ggrepel, quietly = T)){
        install.packages("ggrepel", repos='http://cran.us.r-project.org')
}
if (!require(rlang, quietly = T)){
        install.packages("rlang", repos='http://cran.us.r-project.org')
}
library(ggpubr)
library(plotUtils)
library(argparser)
library(dplyr)
library(ggrepel)
library(rlang)

# Terminal argument parser
################################################################################
parser <- arg_parser("Train the GLM on transformed age and assess performance.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat"),
                       help = c("Predicted ages dataframe of perturbations (generated with model_test.R).",
                                "Metadata file"),
                       flag = c(F,
                                F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
pertFile <- parsed$input
metDatFile <- parsed$metDat

outName <- sprintf("%s/pred_ages_stats.csv",
                   dirname(pertFile))

# Functions
################################################################################

# Plot the fold changes of each perturbation, with colors indicating
# significance and indicating the names of the top N perturbations
plotPerts <- function(pertDF, what_stats, cell_type, alph = 0.05,
                      thrshld_labs_low = NULL,
                      thrshld_labs_high = NULL,
                      topN = NULL,
                      comp_vec = NULL,
                      point_size = 0.1,
                      rem_dose = F,
                      y_var = "log2FC"){
        if (y_var == "log2FC"){
                y_var <- "chron_age_logFC"
        }else if (y_var == "delta"){
                y_var <- "chron_age_delta"
        }else{
                stop(wprintf("%s is not an allowed y variable", y_var), call. = F)
        }
        pertDF <- pertDF[pertDF$cell_type == cell_type, ]
        pertDF <- pertDF[order(pertDF[, y_var]), ]
        pertDF$perturbation <- factor(pertDF$perturbation,
                                      levels = pertDF$perturbation)
        pertDF$sign <- rep(NA, nrow(pertDF))
        pertDF$labs <- rep(NA, nrow(pertDF))
        if(what_stats == "t_test"){
                pertDF$sign[pertDF$chron_age_tTest_pValAdj <= alph] <- "significant"
        }else if (what_stats == "wilcox"){
                pertDF$sign[pertDF$chron_age_wilcox_pValAdj <= alph] <- "significant"
        }
        pertDF$sign[is.na(pertDF$sign)] <- "not significant"
        if (!is.null(thrshld_labs_low) &!is.null(thrshld_labs_high)){
                boolLab <- pertDF$sign == "significant" & (pertDF[, y_var] < thrshld_labs_low | pertDF[, y_var] > thrshld_labs_high)
        }
        if(!is.null(topN) & is.null(comp_vec)){
                pertDF_sign <- pertDF[pertDF$sign == "significant", ]
                topN_sign_idxs <-  setNames(order(abs(pertDF_sign[, y_var]),
                                                  decreasing = T)[1:topN],
                                            pertDF_sign$perturbation[order(abs(pertDF_sign[, y_var]),
                                                                           decreasing = T)[1:topN]])
                boolLab <- rep(F, nrow(pertDF))
                boolLab[pertDF$perturbation %in% names(topN_sign_idxs)] <- T
        }else if (!is.null(comp_vec) & is.null(topN)){c
                pertDF_sign <- pertDF[pertDF$sign == "significant" & pertDF[, y_var] < 0, ]
                pert_names_val <- pertDF_sign$perturbation[grepl(paste(tolower(comp_vec),
                                                                       collapse = "|"),
                                                                 tolower(pertDF_sign$perturbation))] 
                boolLab <- rep(F, nrow(pertDF))
                boolLab[pertDF$perturbation %in% pert_names_val] <- T
        }else if (!is.null(comp_vec) & !is.null(topN)){
                pertDF_sign <- pertDF[pertDF$sign == "significant" & pertDF[, y_var] < 0, ]
                topN_sign_idxs <-  setNames(order(abs(pertDF_sign[, y_var]),
                                                  decreasing = T),
                                            pertDF_sign$perturbation[order(abs(pertDF_sign[, y_var]),
                                                                           decreasing = T)])
                pert_names_val <- topN_sign_idxs[grepl(paste(tolower(comp_vec),
                                                             collapse = "|"),
                                                       tolower(names(topN_sign_idxs)))]
                pert_names_val <- pert_names_val[1:topN]
                boolLab <- rep(F, nrow(pertDF))
                boolLab[pertDF$perturbation %in% names(pert_names_val)] <- T
        }
        pertDF$labs[boolLab] <- as.character(pertDF$perturbation[boolLab])
        pertDF$labs <- gsub(sprintf("_%s", cell_type), "", pertDF$labs)
        if (rem_dose){
                pertDF$labs <- gsub("\\_.*", "", pertDF$labs)
        }
        if (y_var == "chron_age_logFC"){
                y_label <- "predicted age log2 fold change"
        }else if(y_var == "chron_age_delta"){
                y_label <- "Delta predicted age (years)"
        }
        pertDF$labs <- gsub("_", " ", pertDF$labs)
        pertDF <- pertDF %>% arrange(sign)
        plt <- ggplot(pertDF,
                      mapping = aes(x = perturbation,
                                    y = !!sym(y_var),
                                    col = sign,
                                    label = labs)) +
                geom_hline(yintercept = 0) +
                geom_point(size = point_size) +
                scale_color_manual(values = c("red", "darkgreen")) +
                geom_text_repel(data = subset(pertDF, !is.na(labs)),
                                aes(label = labs),
                                max.overlaps = 100) +
                labs(x = "perturbations", y = y_label) +
                scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
                theme(title = element_text(size = 20),
                      axis.text.y = element_text(size=15),
                      axis.text.x = element_blank(),
                      axis.title = element_text(size=20),
                      axis.ticks.x = element_blank(),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      legend.text = element_text(size=12),
                      legend.title = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  linewidth = 1))
        return(plt)
}

# Load the data
################################################################################

perts <- read_table_fast(pertFile, row.names = 1)
metDat <- read_table_fast(metDatFile, row.names = 1)

# Keep only perturbation samples
metDat <- metDat[metDat$diagn_4BrainClck == "perturbations", ]

# Test for statistical significance of unique perturbations
################################################################################
print(sprintf("Running univariate tests of %s (%s)...",
              basename(pertFile),
              basename(dirname(pertFile))))

uniqPerts <- unique(paste(metDat$perturbation,
                          metDat$tissue, sep = "_"))

uniqPerts <- uniqPerts[!grepl("DMSO", uniqPerts)]

metDat$pert_cell <- paste(metDat$perturbation,
                          metDat$tissue, sep = "_")

pert_stats_DF <- data.frame(matrix(nrow = 0,
                                   ncol = 6,
                                   dimnames = list(NULL,
                                                   c("perturbation",
                                                     "cell_type",
                                                     "chron_age_wilcox_pVal",
                                                     "chron_age_tTest_pVal",
                                                     "chron_age_logFC",
                                                     "chron_age_delta"))))

pb = txtProgressBar(min = 0, max = length(uniqPerts), initial = 0, style = 3)
for(i in seq_along(uniqPerts)){
        #i <- 200
        setTxtProgressBar(pb, i)
        p <- uniqPerts[i]
        metDat_p <- metDat[metDat$pert_cell == p & metDat$perturbation != "none", ]
        cTypes <- unique(metDat_p$tissue)
        cTypeVec <- c()
        for(j in seq_along(cTypes)){
                cTyp <- cTypes[j]
                # If it's not lincs, do comparison by same accession. If it's, do comparison
                # taking timempoints into account.

                if (all(metDat_p$timepoint == "")){
                        p_batch <- unique(metDat_p$batch_rna)
                        ctrls_IDs <- metDat$specimenID[metDat$tissue == cTyp & metDat$exper_group == "ctrl" & metDat$batch_rna == p_batch]
                }else{
                        timepoint <- unique(metDat_p$timepoint)
                        ctrls_IDs <- metDat$specimenID[metDat$tissue == cTyp & metDat$exper_group == "ctrl" & metDat$timepoint == timepoint]
                         #ctrls_IDs <- metDat$specimenID[metDat$tissue == cTyp & metDat$exper_group == "ctrl"]
                }
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

                        chronAge_delta <- median(perts_expr_chronAgeVec) - median(perts_ctrl_chronAgeVec)

                        toBindDF <- data.frame(perturbation = p,
                                               cell_type = cTyp,
                                               chron_age_wilcox_pVal = chronAge_wilcox_pVal,
                                               chron_age_tTest_pVal = chronAge_tTest_pVal,
                                               chron_age_logFC = chronAge_logFC,
                                               chron_age_delta = chronAge_delta)
                        pert_stats_DF <- rbind.data.frame(pert_stats_DF, toBindDF)
                }
        }
}
close(pb)

# Adjust p-values per cell type
pert_stats_DF$chron_age_wilcox_pValAdj <- NA
pert_stats_DF$chron_age_tTest_pValAdj <- NA

uniq_cTypes <- unique(pert_stats_DF$cell_type)
for (cTyp in uniq_cTypes){
        pert_stats_DF$chron_age_wilcox_pValAdj[pert_stats_DF$cell_type == cTyp] <- p.adjust(pert_stats_DF$chron_age_wilcox_pVal[pert_stats_DF$cell_type == cTyp],
                                                                                            method = "BH")
        
        pert_stats_DF$chron_age_tTest_pValAdj[pert_stats_DF$cell_type == cTyp] <- p.adjust(pert_stats_DF$chron_age_tTest_pVal[pert_stats_DF$cell_type == cTyp],
                                                                                           method = "BH")
}
#pert_stats_DF$chron_age_wilcox_pValAdj <- p.adjust(pert_stats_DF$chron_age_wilcox_pVal,
#                                                   method = "BH")
#pert_stats_DF$chron_age_tTest_pValAdj <- p.adjust(pert_stats_DF$chron_age_tTest_pVal,
#                                                  method = "BH")


# Save result
################################################################################
write_table_fast(pert_stats_DF, f = outName)

print(sprintf("%s saved at %s", basename(outName), dirname(outName)))

# Plot transcriptional age log2FC and delta per perturbation
################################################################################
npc_neu_pertPlt_vert <- ggarrange(plotPerts(pert_stats_DF,
                                            "t_test",
                                            "NPC",
                                            topN = 10,
                                            rem_dose = T),
                                  plotPerts(pert_stats_DF,
                                            "t_test",
                                            "NEU",
                                            topN = 10,
                                            rem_dose = T),
                                  common.legend = T,
                                  legend = "bottom",
                                  nrow = 2)

ggsave(filename = sprintf("%s/pert_stats_plot.pdf", dirname(pertFile)),
       plot = npc_neu_pertPlt_vert, height = 12, width = 6)

npc_neu_pertPlt_delta_vert <- ggarrange(plotPerts(pert_stats_DF,
                                                  "t_test",
                                                  "NPC",
                                                  topN = NULL,
                                                  thrshld_labs_low = -24,
                                                  thrshld_labs_high = 31.5,
                                                  rem_dose = T,
                                                  y_var = "delta"),
                                        plotPerts(pert_stats_DF,
                                                  "t_test",
                                                  "NEU",
                                                  topN = NULL,
                                                  thrshld_labs_low = -8.7,
                                                  thrshld_labs_high = 19.5,
                                                  rem_dose = T,
                                                  y_var = "delta"),
                                        common.legend = T,
                                        legend = "bottom",
                                        nrow = 2)

ggsave(filename = sprintf("%s/pert_stats_plot_delta.pdf", dirname(pertFile)),
       plot = npc_neu_pertPlt_delta_vert, height = 12, width = 6)