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
                               "--metDat"),
                       help = c("Predicted ages dataframe of sc samples (generated with model_test.R).",
                                "Response variable used to fit the model to (age_chron or age_trans).",
                                "Metadata file"),
                       flag = c(F, F, F))

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

# Directory stuff
################################################################################
#predFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/model_test_sCell/modAllGenes_integWLincs_and_sc_oAge_chronAge_alph1/pred_ages.csv"
#respVar <- "age_chron"
#metDatFile<-"/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast/combined_metDat_wTBI_wPert111_wSC_wLINCS.csv"

predFile <- parsed$input
respVar <- parsed$respVar
metDatFile <- parsed$metDat

outName <- sprintf("%s/pred_ages_stats.csv",
                   dirname(predFile))
outDir <- addSlashIfNot(dirname(outName))

# Load the data
################################################################################

preds <- readCsvFast(predFile)
metDat <- readCsvFast(metDatFile)

# Keep only ageAnno samples
metDat <- metDat[metDat$substudy == "ageAnno", ]

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
                                         ncol = 2,
                                         dimnames = list(NULL,
                                                         c("cell_type",
                                                           "pVal_chron_age"))))
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
                cellStatsDF <- rbind.data.frame(cellStatsDF, toBind)
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
}else if(respVar == "age_trans"){
        pltChronLst <- list()
        pltTransLst <- list()
        uniqCells <- unique(unique(metDat$tissue))
        cellStatsDF <- data.frame(matrix(nrow = 0,
                                         ncol = 3,
                                         dimnames = list(NULL,
                                                         c("cell_type",
                                                           "pVal_chron_age",
                                                           "pVal_trans_age"))))
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
                cellStatsDF <- rbind.data.frame(cellStatsDF, toBind)
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
}

writeCsvFst(cellStatsDF, file = outName)
print(sprintf("%s saved at %s.", basename(outName), dirname(outName)))