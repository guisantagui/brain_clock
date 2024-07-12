################################################################################
# Plot score/biplots given RDS object created with bigPCA_exe.R.               #
################################################################################

if(!require(ggplot2, quietly = T)){
        install.packages("ggplot2",
                         repos='http://cran.us.r-project.org')
}
if(!require(ggtext, quietly = T)){
        install.packages("ggtext",
                         repos='http://cran.us.r-project.org')
}
library(ggtext)
if(!require("ggpubr", quietly = T)){
        install.packages("ggpubr", repos='http://cran.us.r-project.org')
}
library(ggpubr)
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Run PCA for big matrices.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--x",
                               "--y",
                               "--outDir"),
                       help = c("RDS file created with bigPCA_exe.R.",
                                "Metadata CSV file.",
                                "PC in the X axis.",
                                "PC in the Y axis",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

bigPCAFile <- parsed$input
metDatFile <- parsed$metDat
outDir <- parsed$outDir
x <- parsed$x
y <- parsed$y

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

outName <- gsub(".rds|.RDS", "", basename(bigPCAFile))

# Functions
################################################################################
pcBiplot <- function(PC, x="PC1", y="PC2", varPlotFilt = NULL, biPlot = F,
                     colBy = "tissue", colROSMAPBatch = F){
        #metDat <- metDat
        if(colROSMAPBatch){
                metDat$substudy[metDat$substudy == "ROSMAP"] <- paste("ROSMAP",
                                                                      metDat$batch_seq[metDat$substudy == "ROSMAP"],
                                                                      sep = "_")
        }
        data <- data.frame(obsnames=row.names(PC$x), PC$x)
        data <- data[, c("obsnames", x, y)]
        
        data$tissue <- metDat$tissue[match(rownames(data),
                                           make.names(metDat$specimenID))]
        
        data$substudy <- metDat$substudy[match(rownames(data),
                                               make.names(metDat$specimenID))]
        data$pmi <- metDat$pmi[match(rownames(data),
                                     make.names(metDat$specimenID))]
        data$individualID <- metDat$individualID[match(rownames(data),
                                                       make.names(metDat$specimenID))]
        data$age_death <- metDat$ageDeath[match(rownames(data),
                                                make.names(metDat$specimenID))]
        data$diagnosis <- metDat$diagn_4BrainClck[match(rownames(data),
                                                  make.names(metDat$specimenID))]
        propVar <- PC$summary[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        plot <- ggplot(data, aes_string(x = x, 
                                        y = y, 
                                        color = colBy)) + 
                #scale_discrete_manual("Time",
                #                      aesthetics = "colour",
                #                      values = c("cyan",
                #                                 "dodgerblue",
                #                                 "dodgerblue4")) +
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() + 
                xlab(sprintf("%s (%s%%)", x, propX)) +
                ylab(sprintf("%s (%s%%)", y, propY)) +
                #geom_text_repel() +
                #theme_minimal()
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if(biPlot){
                datapc <- data.frame(varnames=rownames(PC$rotation), 
                                     PC$rotation)
                mult <- min(
                        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
                        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
                )
                datapc <- transform(datapc,
                                    v1 = .7 * mult * (get(x)),
                                    v2 = .7 * mult * (get(y))
                )
                datapc$x0 <- rep(0, nrow(datapc))
                datapc$y0 <- rep(0, nrow(datapc))
                if(!is.null(varPlotFilt)){
                        datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
                }
                plot <- plot +
                        geom_text_repel(data=datapc, 
                                        aes(x=v1, y=v2, label=varnames), 
                                        color = "black", 
                                        size = 3) + 
                        geom_segment(data = datapc, aes(x=x0, 
                                                        y=y0, 
                                                        xend=v1, 
                                                        yend=v2, 
                                                        label = varnames),
                                     arrow = arrow(length=unit(0.2,"cm"),
                                                   type = "closed",
                                                   angle = 20), 
                                     alpha=0.75, 
                                     color="black", 
                                     size = 0.5)
        }
        return(plot)
}

# Plots Multi PC score/biplot
doPCAMultiPlot <- function(PC,
                           nComps,
                           varPlotFilt = NULL,
                           biPlot = F,
                           colBy = "tissue",
                           colROSMAPBatch = F){
        plotList <- list()
        for(j in 2:(nComps)){
                for(i in 1:(nComps - 1)){
                        if(j > i){
                                scPlot <- pcBiplot(PC,
                                                   biPlot = F,
                                                   x = sprintf("PC%s", i),
                                                   y = sprintf("PC%s", j),
                                                   varPlotFilt = varPlotFilt,
                                                   colBy = colBy,
                                                   colROSMAPBatch = colROSMAPBatch)
                                if(j < nComps){
                                        scPlot <- scPlot +
                                                theme(axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                                }
                                if(i > 1){
                                        scPlot <- scPlot +
                                                theme(axis.title.y = element_blank(),
                                                      axis.text.y = element_blank())
                                }
                        }else{
                                scPlot <- NA
                        }
                        plotList[[sprintf("PC%s_PC%s", i, j)]] <- scPlot
                }
        }
        multPlot <- ggarrange(plotlist = plotList,
                              common.legend = T,
                              ncol = nComps - 1,
                              nrow = nComps - 1,
                              widths = c(1, rep(.8, nComps-2)),
                              heights = c(rep(.8, nComps-2), 1))
        return(multPlot)
}

# Load the data
################################################################################
bigPCA <- readRDS(bigPCAFile)
metDat <- read.csv(metDatFile, row.names = 1)


# Do the plots and save them
################################################################################
pca_substudy <- pcBiplot(bigPCA,
                         x = x,
                         y = y,
                         colBy = "substudy",
                         biPlot = F,
                         colROSMAPBatch = T)
ggsave(plot = pca_substudy,
       filename = sprintf("%s%s_substudy.pdf",
                          outDir,
                          outName))


pca_tissue <- pcBiplot(bigPCA,
                       x = x,
                       y = y,
                       colBy = "tissue",
                       biPlot = F)

ggsave(plot = pca_tissue,
       filename = sprintf("%s%s_tissue.pdf",
                          outDir,
                          outName))


pca_age_death <- pcBiplot(bigPCA,
                          x = x,
                          y = y,
                          colBy = "age_death",
                          biPlot = F)

ggsave(plot = pca_age_death,
       filename = sprintf("%s%s_ageDeath.pdf",
                          outDir,
                          outName))


pca_pmi <- pcBiplot(bigPCA,
                    x = x,
                    y = y,
                    colBy = "pmi",
                    biPlot = F)

ggsave(plot = pca_pmi,
       filename = sprintf("%s%s_pmi.pdf",
                          outDir,
                          outName))

pca_diagn <- pcBiplot(bigPCA,
                      x = x,
                      y = y,
                      colBy = "diagnosis",
                      biPlot = F)

ggsave(plot = pca_diagn,
       filename = sprintf("%s%s_diagn.pdf",
                          outDir,
                          outName))


pcaMult_subst <- doPCAMultiPlot(bigPCA,
                                nComps = 5,
                                varPlotFilt = NULL,
                                biPlot = F,
                                colBy = "substudy",
                                colROSMAPBatch = T)

ggsave(plot = pcaMult_subst,
       filename = sprintf("%s%s_mult_substudy.pdf",
                          outDir,
                          outName))


pcaMult_tiss <- doPCAMultiPlot(bigPCA,
                               nComps = 5,
                               varPlotFilt = NULL,
                               biPlot = F,
                               colBy = "tissue")

ggsave(plot = pcaMult_tiss,
       filename = sprintf("%s%s_mult_tissue.pdf",
                          outDir,
                          outName))


pcaMult_age <- doPCAMultiPlot(bigPCA,
                              nComps = 5,
                              varPlotFilt = NULL,
                              biPlot = F,
                              colBy = "age_death")

ggsave(plot = pcaMult_age,
       filename = sprintf("%s%s_mult_ageDeath.pdf",
                          outDir,
                          outName))


pcaMult_pmi <- doPCAMultiPlot(bigPCA,
                              nComps = 5,
                              varPlotFilt = NULL,
                              biPlot = F,
                              colBy = "pmi")

ggsave(plot = pcaMult_pmi,
       filename = sprintf("%s%s_mult_pmi.pdf",
                          outDir,
                          outName))


pcaMult_diagn <- doPCAMultiPlot(bigPCA,
                                nComps = 5,
                                varPlotFilt = NULL,
                                biPlot = F,
                                colBy = "diagnosis")

ggsave(plot = pcaMult_diagn,
       filename = sprintf("%s%s_mult_diagn.pdf",
                          outDir,
                          outName))