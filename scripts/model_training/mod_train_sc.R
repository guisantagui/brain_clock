################################################################################
#                                                                              #
# Brain Clock: fit a GLM model to the processed SC data from ageAnno to find   #
# which are the features that allow for better prediction in individual cell   #
# These features will be added to a second round of fitting to ensure that the #
# model captures differences in individual cell types. Use chronological age.  #
#                                                                              #
################################################################################

if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require(limma, quietly = T)) BiocManager::install("limma", update = F)
if(!require("h2o", quietly = T)){
        install.packages("h2o",
                         type="source",
                         repos="https://h2o-release.s3.amazonaws.com/h2o/rel-3.46.0/2/R")
}
library(h2o)
if(!require("ggplot2")) install.packages("ggplot2",
                                         repos='https://pbil.univ-lyon1.fr/CRAN/')
library(ggplot2)
if(!require("ggpubr")) install.packages("ggpubr",
                                        repos='https://pbil.univ-lyon1.fr/CRAN/')
library(ggpubr)

# Terminal argument parser
################################################################################
parser <- arg_parser("Train the GLM on transformed age and assess performance.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--alpha",
                               "--mem",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features columns, samples rows. CSV.",
                                "Metadata file (including ages).",
                                "Alpha value for the elastic net.",
                                "Memory to allocate to h2o instance.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# Read csv faster
readCsvFast <- function(f){
        df <- data.frame(data.table::fread(f))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if(!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Predict the age viven the expression matrix and the model
predictAge <- function(model,
                       predExpr
){
        predExpr <- predExpr[model@parameters$x,]
        predAges <- c()
        maxiter <- ceiling(ncol(predExpr)/1000)
        for(i in 1:maxiter){
                #i <- 1
                from <- ((i-1)*1000+1)
                to <- min(i*1000,ncol(predExpr))
                
                inp_h2o <- as.h2o(t(predExpr[,from:to]))
                
                pred <- h2o.predict(model,inp_h2o)
                pred <- as.data.frame(pred)$predict
                names(pred) <- colnames(predExpr)[from:to]
                
                predAges <- c(predAges,pred)
        }
        return(predAges)
}

# Given a vector of values, and the predicted values according to the model,
# compute some metrics.
getMetrics <- function(vals, pred){
        mae <- sum(abs(vals - pred))/length(vals)
        mse <- sum((vals - pred)^2)/length(vals)
        rmse <- sqrt(mse)
        mape <- (sum(abs((vals - pred)/vals))/length(vals))*100
        r2 <- 1 - sum((vals - pred)^2)/sum((vals - mean(vals))^2)
        out <- data.frame(metric = c("mae", "mse", "rmse", "mape", "r2"),
                          value = c(mae, mse, rmse, mape, r2))
        return(out)
}

# Plot ages vs predicted ages in test control dataset and in the high braak
# dataset
comp_predAges_test_ND <- function(test_ages, nd_ages, test_pred, nd_pred,
                                  color = "group", filtCommAgeRank = T){
        df <- data.frame(age_death = c(test_ages, nd_ages),
                         specimenID = names(c(test_ages, nd_ages)),
                         pred_age = c(test_pred, nd_pred),
                         group = c(rep("control", length(test_ages)),
                                   rep("nd", length(nd_ages))))
        df$braak <- metDat$Braak[match(make.names(df$specimenID),
                                       make.names(metDat$specimenID))]
        df$substudy <- metDat$substudy[match(make.names(df$specimenID),
                                             make.names(metDat$specimenID))]
        if(filtCommAgeRank){
                minAgeND <- min(df$age_death[df$group == "nd"])
                df <- df[df$age_death >= minAgeND, ]
        }
        plt <- ggplot(data = df, mapping = aes(x = age_death,
                                               y = pred_age,
                                               col = .data[[color]])) +
                geom_point() + 
                xlab("chronological age") +
                ylab("predicted age") +
                theme(axis.text.y = element_text(size=15),
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
                                                  fill=NA, linewidth = 1))
        out <- list(data = df, plot = plt)
        return(out)
}

# Directory stuff
################################################################################
dataFile <- "/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast_allLINCSBrain/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_NPC_NEU_MIC_noCerebell_onlyAge_svaAdj.csv"
metDatFile <- "/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wSC_wLINCS_NPC_NEU_MIC.csv"
alph <- 1
mem  <- "24G"
excludeYoung <- T
outDir <- "/home/users/gsantamaria/projects/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_SCmod_chron_age/"

dataFile <- parsed$input
metDatFile <- parsed$metDat
alph <- as.numeric(parsed$alpha)
mem <- parsed$mem
outDir <- addSlashIfNot(parsed$outDir)

createIfNot(outDir)

# Load and parse data
################################################################################

metDat <- readCsvFast(metDatFile)

dat <- readCsvFast(dataFile)

metDat <- metDat[metDat$substudy == "ageAnno", ]
if(excludeYoung){
        metDat <- metDat[metDat$ageDeath >= 18, ]
}

dat <- dat[make.names(rownames(dat)) %in% make.names(metDat$specimenID), ]


# Split in train and test
################################################################################

# Ensure that all cell types are represented, and that the 3 age stages (young,
# 18 to 40, mid, 40 to 70, and old 60 onwards) are represented.

trainProp <- .66

cellTypes <- names(table(metDat$tissue))

ageInts <- c("18-40", "40-70", "70-100")

seed <- 123

sampsInTrain <- c()
for(i in seq_along(cellTypes)){
        cell <- cellTypes[i]
        samps_cell <- metDat$specimenID[metDat$tissue == cell]
        for(j in seq_along(ageInts)){
                agInt <- ageInts[j]
                ag1 <- as.numeric(strsplit(agInt, split = "-")[[1]][1])
                ag2 <- as.numeric(strsplit(agInt, split = "-")[[1]][2])
                samps_agInt <- metDat$specimenID[metDat$ageDeath >= ag1 & metDat$ageDeath < ag2]
                samps_cell_agInt <- samps_cell[samps_cell %in% samps_agInt]
                nToSamp <- round(length(samps_cell_agInt) * trainProp)
                set.seed(seed + i + j)
                sampsInTrain <- c(sampsInTrain,
                                  sample(samps_cell_agInt, nToSamp))
        }
}

dat_train <- dat[rownames(dat) %in% make.names(sampsInTrain), ]
dat_test <- dat[!rownames(dat) %in% make.names(sampsInTrain), ]


ageChronTrain <- metDat$ageDeath[match(rownames(dat_train),
                                       make.names(metDat$specimenID))]
train4Mod <- cbind.data.frame(ageChronTrain, dat_train)

colnames(train4Mod)[1] <- "age_chron"

# Model training
################################################################################

# Initialize h2o and fit the model
conn <- h2o.init(max_mem_size=mem)
trainData_h2o <- as.h2o(train4Mod)

lambda <- NULL
lambda_search <- T

model <- h2o.glm(y = "age_chron",
                 training_frame = trainData_h2o,
                 nfolds = 10,
                 fold_assignment = "Random",
                 family = "AUTO",
                 link = "family_default",
                 lambda_search = lambda_search,
                 lambda = lambda,
                 standardize = T,
                 alpha = alph,
                 #alpha = 0,
                 seed = 100,
                 max_active_predictors = ncol(trainData_h2o),
                 solver = "IRLSM")

h2o.saveModel(model, path = sprintf("%smodFuncsAlpha%s", outDir, alph), force = T)
saveRDS(model, file = sprintf("%smodFuncsAlpha%s.rds", outDir, alph))

mod_coefs <- model@model$coefficients_table

mod_coefs <- mod_coefs[mod_coefs$coefficients != 0, ]

rownames(mod_coefs) <- mod_coefs$names
mod_coefs <- mod_coefs[!grepl("Intercept", mod_coefs$names), ]
colnames(mod_coefs) <- gsub("names", "ensembl_gene_id", colnames(mod_coefs))

coefsName <- sprintf("%smodFuncsAlpha%s_coefs.csv", outDir, as.character(alph))
write.csv(mod_coefs, coefsName)
print(sprintf("%s saved at %s", basename(coefsName), dirname(coefsName)))

# Model evaluation in test dataset
################################################################################

# Plot metrics
cv_r2 <- unlist(as.vector(model@model$cross_validation_metrics_summary["r2", -c(1, 2)]))
cv_mae <- unlist(as.vector(model@model$cross_validation_metrics_summary["mae", -c(1, 2)]))
cv_rmse <- unlist(as.vector(model@model$cross_validation_metrics_summary["rmse", -c(1, 2)]))

plotDF <- data.frame(metric = c(rep("r2", length(cv_r2)),
                                rep("mae", length(cv_mae)),
                                rep("rmse", length(cv_rmse))),
                     value = c(cv_r2, cv_mae, cv_rmse))

                     ggplot(data = plotDF, mapping = aes(x = metric, y = value)) +
        geom_boxplot(outlier.shape = NA) + geom_jitter() +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              axis.title = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, linewidth = 1))

ggsave(filename = sprintf("%smodFuncsAlpha%s_cvMetrics.pdf", outDir, alph))

# Predict the age in the test dataset
predicted <- predictAge(model, t(dat_test))

testAges <- metDat$ageDeath[match(names(predicted), make.names(metDat$specimenID))]
metrics_chronAge <- getMetrics(testAges, predicted)
print("Metrics for test dataset on chronological ages:")
print(metrics_chronAge)
write.csv(metrics_chronAge,
          file = sprintf("%smetricsTest_chronAge_alpha%s.csv", outDir, alph))

# Plot R2 of training, CV and testing
dfr2 <- data.frame(r2_type = factor(c("r2_training",
                                      rep("r2_cv",
                                          length(model@model$cross_validation_metrics_summary["r2", -c(1, 2)])), "r2_test"),
                                    levels = c("r2_training",
                                               "r2_cv",
                                               "r2_test")),
                   value = c(h2o.r2(model, train = T),
                             unlist(model@model$cross_validation_metrics_summary["r2", -c(1, 2)]),
                             metrics_chronAge$value[metrics_chronAge$metric == "r2"]))

dfr2_bxplt <- dfr2[dfr2$r2_type == "r2_cv", ]
dfr2_bxplt$r2_type <- factor(dfr2_bxplt$r2_type, levels = c("r2_training",
                                                            "r2_cv",
                                                            "r2_test"))
dfr2_point <- dfr2[dfr2$r2_type != "r2_cv", ]
dfr2_point$r2_type <- factor(dfr2_point$r2_type, levels = c("r2_training",
                                                            "r2_cv",
                                                            "r2_test"))

ggplot(data = dfr2, mapping = aes(x = r2_type, y = value)) +
        geom_point(data = dfr2_point) + ylim(0, 1) +
        geom_boxplot(data = dfr2_bxplt, outlier.shape = NA) +
        geom_jitter(data = dfr2_bxplt) +
        scale_x_discrete(limits = c("r2_training",
                                    "r2_cv",
                                    "r2_test"),
                         labels = c("r2_training" = "R2 training",
                                    "r2_cv" = "R2 CV",
                                    "r2_test" = "R2 test")) +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              axis.title = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, linewidth = 1))

ggsave(sprintf("%smodelAlph%sR2.pdf", outDir, as.character(alph)),
       height = 3, width = 4)

# Plot the chronological age vs the predicted age for the test set, coloring 
# according to the cell type.

df_chron_pred <- data.frame(sample = names(predicted),
                            cell_type = gsub("\\_.*", "", names(predicted)),
                            chron_age = testAges,
                            pred_age = predicted)
plt <- ggplot(df_chron_pred, mapping = aes(x = chron_age, y = pred_age, col = cell_type)) +
        geom_point() +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              axis.title = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, linewidth = 1))

ggsave(filename = sprintf("%schron_vs_pred.pdf", outDir), plot = plt)

h2o.shutdown(prompt = F)