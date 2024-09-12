################################################################################
# Brain clock: test for statistical significance of the relationship between   #
# chronological age and predicted age in single cell data from ageAnno. In     #
# particular, fit linear model for pred_age ~ chron_age for each cell type     #
# and assess significance via F-test.                                          #
################################################################################
library(ggplot2)
library(ggpubr)
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
                               "--metDat",
                               "--excludeYoung"),
                       help = c("Predicted ages dataframe of sc samples (generated with model_test.R).",
                                "Response variable used to fit the model to (age_chron or age_trans).",
                                "Metadata file",
                                "If samples below 18 YO should be removed prior of running the analysis."),
                       flag = c(F, F, F, T))

parsed <- parse_args(parser)

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

# Fit a linear model using desired response variable to age_death
# (actual chronological age) andn return the model, the pValue and 
# a plot.
getLMod <- function(df, respVar, cellType = NULL){
        form <- as.formula(sprintf("%s ~ age_death", respVar))
        lMod <- lm(form, df)
        fStat <- summary(lMod)$f
        pVal <- pf(fStat[1], fStat[2], fStat[3], lower.tail = F)
        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p-value)~"="~pVal, 
                         list(a = format(unname(coef(lMod)[1]), digits = 2),
                              b = format(unname(coef(lMod)[2]), digits = 2),
                              r2 = format(summary(lMod)$r.squared, digits = 3),
                              pVal = format(unname(pVal), digits = 3)))
        eq <- as.character(as.expression(eq))
        eq <- paste0("atop(",
                     gsub("(x) * \",\"", "(x),", eq, fixed = T),
                     ")")
        plt <- ggplot(data = df, mapping = aes(y = !!sym(respVar), x = age_death)) +
                geom_point() +
                geom_smooth(method = "lm", se = F) +
                labs(x = "chronological age",
                     y = sprintf("predicted %s",
                                 gsub("_", " ", respVar))) +
                theme(title = element_text(size = 20),
                      axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      axis.title = element_text(size=20),
                      panel.background = element_blank(),
                      panel.grid.major = element_line(colour = "gray"), 
                      panel.grid.minor = element_blank(),
                      legend.text = element_text(size=12),
                      legend.title = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  linewidth = 1))
        xRange <- layer_scales(plt)$x$range$range
        yRange <- layer_scales(plt)$y$range$range
        if(lMod$coefficients[2] > 0.3){
                xPos <- min(xRange) + round((abs(min(xRange)) + 40) * .8)
                yPos <- max(yRange) - round(abs(max(yRange)) * .1)
        }else if(lMod$coefficients[2] <= 0.3 & lMod$coefficients[2] > 0){
                xPos <- round(mean(xRange))
                yPos <- max(yRange) - round(abs(max(yRange)) * .051)
        }else{
                xPos <- max(xRange) - round((abs(max(xRange)) + 10) * .8)
                yPos <- max(yRange) - round(abs(max(yRange)) * .1)
        }
        plt <- plt +
                geom_label(x = xPos,
                           y = yPos,
                           label = eq, parse = TRUE)
        if(!is.null(cellType)){
                plt <- plt +
                        ggtitle(cellType)
        }
        out <- list(lMod = lMod, pVal = pVal, plot = plt)
        return(out)
}

# Run young vs old univariate test and plot it
run_young_vs_old <- function(df, respVar, youngThrshld = 30, oldThrshld = 70,
                             test = "wilcox.test", cell){
        df$group <- rep(NA, nrow(df))
        df$group[df$age_death <= youngThrshld] <- "young"
        df$group[df$age_death >= oldThrshld] <- "old"
        df <- df[!is.na(df$group), ]
        df$group <- factor(df$group, levels = c("young", "old"))
        if(test == "wilcox.test"){
                pVal <- wilcox.test(df[df$group == "young", respVar],
                                    df[df$group == "old", respVar])$p.val
        }else if(test == "t.test"){
                pVal <- t.test(df[df$group == "young", respVar],
                               df[df$group == "old", respVar])$p.val
        }
        
        statSymbDF <- data.frame(value = c(1, 0.1, 0.05, 0.01, 0.001),
                                 symbol = c("ns", ".", "*", "**", "***"))
        symb <- statSymbDF$symbol[max(which(pVal <= statSymbDF$value))]
        maxY <- max(df[, respVar])
        maxY <- maxY + max(maxY) * .025
        signifDF <- data.frame(x = 1,
                               xend = 2,
                               y = maxY,
                               annotation = symb)
        plt <- ggplot(df, aes(x = group, y = !!sym(respVar))) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter() +
                labs(title = gsub(".", " ", cell, fixed = T),
                     y = "predicted age") +
                theme(title = ggtext::element_markdown(size = 20),
                      axis.text.x = ggtext::element_markdown(size = 15),
                      axis.title.y = ggtext::element_markdown(size = 15),
                      axis.title.x = element_blank(),
                      legend.text = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4")) +
                geom_signif(stat = "identity",
                            data = signifDF,
                            aes(x=x, xend=xend,
                                y=y,
                                yend=y,
                                fill = NULL,
                                annotation = annotation,
                                #group = group,
                                col = NULL))
        out <- list(p_value = pVal, plot = plt)
        return(out)
}

# Directory stuff
################################################################################
#predFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/model_test_sCell/modAllGenes_integWLincs_and_sc_oAge_chronAge_alph1/pred_ages.csv"
#respVar <- "age_chron"
#metDatFile<-"/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast/combined_metDat_wTBI_wPert111_wSC_wLINCS.csv"

predFile <- parsed$input
respVar <- parsed$respVar
metDatFile <- parsed$metDat
excludeYoung <- parsed$excludeYoung

outName <- sprintf("%s/pred_ages_stats.csv",
                   dirname(predFile))
outDir <- addSlashIfNot(dirname(outName))

# Load the data
################################################################################

preds <- readCsvFast(predFile)
metDat <- readCsvFast(metDatFile)

# Keep only ageAnno samples
metDat <- metDat[metDat$substudy == "ageAnno", ]

if(excludeYoung){
        keepSamps <- metDat$specimenID[metDat$ageDeath >= 18]
        preds <- preds[preds$specimenID %in% keepSamps, ]
        metDat <- metDat[metDat$specimenID %in% preds$specimenID, ]
        #metDat <- metDat[!grepl("mid7", metDat$specimenID), ] #Try to remove mid7
        #metDat <- metDat[!grepl(paste(ols, collapse = "|"), metDat$specimenID), ]
        preds <- preds[preds$specimenID %in% metDat$specimenID, ]
}
# Test for statistical significance of association between chronological age
# and predicted age for each cell type
################################################################################
print(sprintf("Running univariate tests for each cell type in %s (%s)...",
              basename(predFile),
              basename(dirname(predFile))))

if(respVar == "age_chron"){
        pltChronLst <- list()
        uniqCells <- unique(unique(metDat$tissue))
        cellStatsDF <- data.frame(matrix(nrow = 0,
                                         ncol = 3,
                                         dimnames = list(NULL,
                                                         c("cell_type",
                                                           "pVal_chron_age",
                                                           "young_vs_old_chron_pVal"))))
        univPltChrnLst <- list()
        for(cell in uniqCells){
                df_cell <- data.frame(specimenID = metDat$specimenID[metDat$tissue == cell],
                                      age_death = metDat$ageDeath[metDat$tissue == cell])
                df_cell$chron_age <- preds$chron_age[match(make.names(df_cell$specimenID),
                                                           make.names(preds$specimenID))]
                lm_chron <- getLMod(df_cell, "chron_age", cellType = cell)
                plotChronName <- sprintf("%slm_chron_age_vs_pred_chron_age_%s.pdf",
                                         outDir,
                                         cell)
                pValChron <- lm_chron$pVal
                plotChron <- lm_chron$plot
                ggsave(plotChronName, plot = plotChron, height = 5, width = 6)
                pltChronLst[[cell]] <- plotChron
                toBind <- data.frame(cell_type = cell,
                                     pVal_chron_age = pValChron)
                univ_chron <- run_young_vs_old(df_cell,
                                               "chron_age",
                                               cell = cell,
                                               test = "wilcox.test",
                                               youngThrshld = 30,
                                               oldThrshld = 70)
                toBind$young_vs_old_chron_pVal <- univ_chron$p_value
                cellStatsDF <- rbind.data.frame(cellStatsDF, toBind)
                univPltChrnLst[[cell]] <- univ_chron$plot
        }
        cellStatsDF$pAdj_chron_age <- p.adjust(cellStatsDF$pVal_chron_age,
                                               method = "BH")
        nrowCombPlot <- round(sqrt(length(uniqCells)))
        ncolCombPlot <- length(uniqCells) / nrowCombPlot
        chronPlotsComb <- ggarrange(plotlist = pltChronLst,
                                    ncol = ncolCombPlot,
                                    nrow = nrowCombPlot)
        plotChronAllCellsName <- sprintf("%slm_chron_age_vs_pred_chron_age_allCellTyp.pdf",
                                         outDir)
        ggsave(plotChronAllCellsName, chronPlotsComb, height = 8, width = 12)
        
        chronUnivPlotsComb <- ggarrange(plotlist = univPltChrnLst,
                                        ncol = ncolCombPlot,
                                        nrow = nrowCombPlot)
        chronUnivPlotsName <- sprintf("%suniv_chron_age_young_vs_old.pdf",
                                      outDir)
        ggsave(chronUnivPlotsName, chronUnivPlotsComb, height = 10, width = 10)
}else if(respVar == "age_trans"){
        pltChronLst <- list()
        pltTransLst <- list()
        univPltChrnLst <- list()
        univPltTrnsLst <- list()
        uniqCells <- unique(unique(metDat$tissue))
        cellStatsDF <- data.frame(matrix(nrow = 0,
                                         ncol = 5,
                                         dimnames = list(NULL,
                                                         c("cell_type",
                                                           "pVal_chron_age",
                                                           "pVal_trans_age",
                                                           "young_vs_old_chron_pVal",
                                                           "young_vs_old_trans_pVal"))))
        for(cell in uniqCells){
                df_cell <- data.frame(specimenID = metDat$specimenID[metDat$tissue == cell],
                                      age_death = metDat$ageDeath[metDat$tissue == cell])
                df_cell$chron_age <- preds$chron_age[match(make.names(df_cell$specimenID),
                                                           make.names(preds$specimenID))]
                df_cell$trans_age <- preds$trans_age[match(make.names(df_cell$specimenID),
                                                           make.names(preds$specimenID))]
                lm_chron <- getLMod(df_cell, "chron_age", cellType = cell)
                lm_trans <- getLMod(df_cell, "trans_age", cellType = cell)
                plotChronName <- sprintf("%slm_chron_age_vs_pred_chron_age_%s.pdf",
                                         outDir,
                                         cell)
                plotTransName <- sprintf("%slm_chron_age_vs_pred_trans_age_%s.pdf",
                                         outDir,
                                         cell)
                pValChron <- lm_chron$pVal
                pValTrans <- lm_trans$pVal
                plotChron <- lm_chron$plot
                plotTrans <- lm_trans$plot
                ggsave(plotChronName, plot = plotChron, height = 5, width = 6)
                ggsave(plotTransName, plot = plotTrans, height = 5, width = 6)
                pltChronLst[[cell]] <- plotChron
                pltTransLst[[cell]] <- plotTrans
                toBind <- data.frame(cell_type = cell,
                                     pVal_chron_age = pValChron,
                                     pVal_trans_age = pValTrans)
                univ_chron <- run_young_vs_old(df_cell,
                                               "chron_age",
                                               cell = cell,
                                               test = "wilcox.test",
                                               youngThrshld = 30,
                                               oldThrshld = 70)
                univ_trans <- run_young_vs_old(df_cell,
                                               "trans_age",
                                               cell = cell,
                                               test = "wilcox.test",
                                               youngThrshld = 30,
                                               oldThrshld = 70)
                toBind$young_vs_old_chron_pVal <- univ_chron$p_value
                toBind$young_vs_old_trans_pVal <- univ_trans$p_value
                cellStatsDF <- rbind.data.frame(cellStatsDF, toBind)
                univPltChrnLst[[cell]] <- univ_chron$plot
                univPltTrnsLst[[cell]] <- univ_trans$plot
        }
        cellStatsDF$pAdj_chron_age <- p.adjust(cellStatsDF$pVal_chron_age,
                                               method = "BH")
        cellStatsDF$pAdj_trans_age <- p.adjust(cellStatsDF$pVal_trans_age,
                                               method = "BH")
        nrowCombPlot <- round(sqrt(length(uniqCells)))
        ncolCombPlot <- length(uniqCells) / nrowCombPlot
        chronPlotsComb <- ggarrange(plotlist = pltChronLst,
                                    ncol = ncolCombPlot,
                                    nrow = nrowCombPlot)
        transPlotsComb <- ggarrange(plotlist = pltTransLst,
                                    ncol = ncolCombPlot,
                                    nrow = nrowCombPlot)
        plotChronAllCellsName <- sprintf("%slm_chron_age_vs_pred_chron_age_allCellTyp.pdf",
                                         outDir)
        plotTransAllCellsName <- sprintf("%slm_chron_age_vs_pred_trans_age_allCellTyp.pdf",
                                         outDir)
        ggsave(plotChronAllCellsName, chronPlotsComb, height = 8, width = 12)
        ggsave(plotTransAllCellsName, transPlotsComb, height = 8, width = 12)
        
        chronUnivPlotsComb <- ggarrange(plotlist = univPltChrnLst,
                                        ncol = ncolCombPlot,
                                        nrow = nrowCombPlot)
        transUnivPlotsComb <- ggarrange(plotlist = univPltTrnsLst,
                                        ncol = ncolCombPlot,
                                        nrow = nrowCombPlot)
        transUnivPlotsName <- sprintf("%suniv_trans_age_young_vs_old.pdf",
                                      outDir)
        ggsave(transUnivPlotsName, transUnivPlotsComb, height = 10, width = 10)
}

writeCsvFst(cellStatsDF, file = outName)
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))