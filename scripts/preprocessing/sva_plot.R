################################################################################
# Plot SVA surrogate variables colored by known batches to assess if there are #
# batch effects left and if the surrogate variables are capturing known.       #
# batches.                                                                     #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require("ggplot2", quietly = T)){
        install.packages("ggplot2", repos='http://cran.us.r-project.org')
}
library(ggplot2)
if(!require("ggpubr", quietly = T)){
        install.packages("ggpubr", repos='http://cran.us.r-project.org')
}
library(ggpubr)
if(!require(sva, quietly = T)) BiocManager::install("sva", update = F)
library(sva)

# Terminal argument parser
################################################################################
parser <- arg_parser("Plot surrogate variables from SVA")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--sva",
                               "--metDat",
                               "--outDir"),
                       help = c("Input matrix, features in columns, samples rows. CSV.",
                                "SVA object (in RDS format)",
                                "Metadata CSV file including batch info.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F))

parsed <- parse_args(parser)


# Directory stuff
################################################################################
datFile <- parsed$input
svaFile <- parsed$sva
metDatFile <- parsed$metDat
outDir <- parsed$outDir

outBaseName <- basename(svaFile)
outBaseName <- gsub(".rds|.RDS", "", outBaseName)
outName <- sprintf("%s%s_plt.pdf", outDir, outBaseName)

# Load data
################################################################################
svaObj <- readRDS(svaFile)
dat <- data.frame(data.table::fread(datFile))
rownames(dat) <- dat$V1
dat <- dat[, colnames(dat) != "V1"]
metDat <- read.csv(metDatFile, row.names = 1)


# Plot the surrogate variables
################################################################################
plotSVs <- function(df, x, y){
        plt <- ggplot(data = df, mapping = aes_string(x = x,
                                                      y = y,
                                                      col = "batch")) +
                geom_point() +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        return(plt)
}

plotAllSVs <- function(svObj, dat, metDat){
        svMat <- svObj$sv
        rownames(svMat) <- rownames(dat)
        colnames(svMat) <- paste("SV", as.character(1:ncol(svMat)), sep = "_")
        svMat <- data.frame(svMat)
        
        svMat$batch <- metDat$substudy[match(rownames(svMat), make.names(metDat$specimenID))]
        plotList <- list()
        for(j in 2:(ncol(svMat) - 1)){
                for(i in 1:(ncol(svMat) - 2)){
                        if(j > i){
                                svPlot <- plotSVs(svMat,
                                                  x = sprintf("SV_%s", i),
                                                  y = sprintf("SV_%s", j))
                                if(j < (ncol(svMat) - 1)){
                                        svPlot <- svPlot +
                                                theme(axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                                }
                                if(i > 1){
                                        svPlot <- svPlot +
                                                theme(axis.title.y = element_blank(),
                                                      axis.text.y = element_blank())
                                }
                        }else{
                                svPlot <- NA
                        }
                        plotList[[sprintf("SV_%s_SV_%s", i, j)]] <- svPlot
                }
        }
        multPlot <- ggarrange(plotlist = plotList,
                              common.legend = T,
                              ncol = ncol(svMat) - 2,
                              nrow = ncol(svMat) - 2,
                              widths = c(1, rep(.8, ncol(svMat)-3)),
                              heights = c(rep(.8, ncol(svMat)-3), 1))
        return(multPlot)
}

allSVs <- plotAllSVs(svObj = svaObj, dat = dat, metDat = metDat)
ggsave(plot = allSVs, filename = outName, width = 10, height = 10)