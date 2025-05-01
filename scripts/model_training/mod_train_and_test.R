################################################################################
#                                                                              #
# Brain Clock: fit a GLM model to the processed brain expression dataset to    #
# predict age, and assess if there are differences in predicted ages in the    #
# neurodegenerated subjects                                                    #
#                                                                              #
################################################################################

if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}
if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
if(!require(rstatix, quietly = T)){
        install.packages("rstatix", repos='http://cran.us.r-project.org')
}
if(!require(limma, quietly = T)) BiocManager::install("limma", update = F)
if(!require("h2o", quietly = T)){
        install.packages("h2o",
                         type="source",
                         repos="https://h2o-release.s3.amazonaws.com/h2o/rel-3.46.0/2/R")
}
if(!require("ggplot2")) install.packages("ggplot2",
                                         repos='https://pbil.univ-lyon1.fr/CRAN/')
if(!require("ggpubr")) install.packages("ggpubr",
                                        repos='https://pbil.univ-lyon1.fr/CRAN/')

library(argparser)
library(rstatix)
library(limma)
library(h2o)
library(ggplot2)
library(ggpubr)
library(plotUtils)

# Terminal argument parser
################################################################################
parser <- arg_parser("Train the GLM on transformed age and assess performance.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--alpha",
                               "--mem",
                               "--preFiltGenes",
                               "--lambda",
                               "--braakThrshld",
                               "--refit",
                               "--overwrite_mod",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features columns, samples rows. CSV.",
                                "Metadata file (including ages).",
                                "Alpha value for the elastic net.",
                                "Memory to allocate to h2o instance.",
                                "CSV file with genes to be prefiltered previous to the fitting. Has to have a column called 'ensembl_gene_id'. If set to 'none' all the genes will be used",
                                "If set, lambda search will be done. Otherwise it will be set to zero and no regularization will be done.",
                                "Braak score threshold for considering samples neurodegenerated.",
                                "If, even when a model is available in the output path, it should be refit",
                                "Force overwrite",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F, T, F, T, T, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

dataFile <- "../../results/preproc/test_no_lincs/merged_counts_log2_quantNorm_noCerebell_combat_onlyAge_svaAdj.csv"
metDatFile <- "../../results/parsed/merged/merged_metdat.csv"
alph <- 1
mem <- "24G"
preFiltGenes <- "none"
braakThrshld <- 4
outDir <- "../../results/models/first_round/"
lambdaFlag <- T
refit <- F
overwrite_mod <- T

dataFile <- parsed$input
metDatFile <- parsed$metDat
alph <- as.numeric(parsed$alpha)
mem <- parsed$mem
preFiltGenes <- parsed$preFiltGenes
lambdaFlag <- parsed$lambda
braakThrshld <- as.numeric(parsed$braakThrshld)
refit <- parsed$refit
overwrite_mod <- parsed$overwrite_mod
outDir <- parsed$outDir

create_dir_if_not(outDir)

respVar <- "age_chron"
# Functions
################################################################################

# Predict the age given the expression matrix and the model
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
        pearson2 <- cor(vals, pred)^2
        bias_mean <- mean(pred - vals)
        bias_median <- median(pred - vals)
        out <- data.frame(metric = c("mae", "mse", "rmse", "mape", "r2",
                                     "pearson2", "bias_mean", "bias_median"),
                          value = c(mae, mse, rmse, mape, r2,
                                    pearson2, bias_mean, bias_median))
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

# Given a dataframe with the predicted ages in each sample, its age interval and
# its braak index, plot predicted age per braak index.
plotAgeInt <- function(dfBraak, ageInt, color = NULL){
        plt <- ggplot(data = dfBraak[dfBraak$ageInt == ageInt, ],
                      mapping = aes(x = braak,
                                    y = pred_age)) +
                geom_boxplot(outlier.shape = NA) +
                ylab("predicted age") +
                xlab("Braak index") +
                ggtitle(gsub("_", " to ", ageInt)) +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_line(colour = "gray"), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA, linewidth = 1))
        if(is.null(color)){
                plt <- plt +
                        geom_jitter()
        }else{
                plt <- plt +
                        geom_jitter(mapping = aes(col = .data[[color]]))
        }
        return(plt)
}

# Given an age interval, dfBraak and outName, fits a linear model on braak
# index and pred_age for individuals not classified and controls, and 
# writes out the summary.
getModPredAgeVsBraak <- function(df,
                                 age_int = NULL,
                                 outName = NULL,
                                 doPlot = T,
                                 savePlt = F,
                                 plotEq = T,
                                 accountChronAge = F,
                                 xPos = NULL,
                                 yPos = NULL){
        #df <- dfBraak
        #age_int <- NULL
        if(!is.null(age_int)){
                df <- df[df$ageInt == age_int, ]
        }
        df <- df[df$braak != "control", ]
        df$braak <- as.numeric(df$braak)
        
        if(accountChronAge){
                lMod <- lm(pred_age ~ braak + age_death, df)
        }else{
                lMod <- lm(pred_age ~ braak, df)
        }
        print(summary(lMod))
        if(!is.null(outName)){
                capture.output(summary(lMod),
                               file = outName)
        }
        if(plotEq){
                if(!accountChronAge){
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
                }else{
                        summ <- summary(lMod)
                        braak_pVal <- summ$coefficients["braak", 4]
                        eq <- substitute("Braak index's" ~~italic(p-value)~"="~braak_pVal,
                                         list(braak_pVal = format(unname(braak_pVal), digits = 3)))
                        eq <- as.character(as.expression(eq))
                }
                
        }
        if(doPlot){
                plt <- ggplot(data = df, mapping = aes(x = braak, y = pred_age)) +
                        geom_point() +
                        geom_smooth(method = "lm", se = FALSE) +
                        labs(x = "Braak index", y = "predicted age") +
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
                if(!is.null(age_int)){
                        tit <- gsub("_", " to ", age_int)
                        plt <- plt +
                                ggtitle(tit)
                }
                if(!is.null(outName) & savePlt){
                        plotName <- gsub(".txt", ".pdf", outName)
                        ggsave(plotName,
                               plot = plt,
                               height = 8,
                               width = 9)
                }
                if(plotEq){
                        if(is.null(xPos) & is.null(yPos)){
                                xRange <- layer_scales(plt)$x$range$range
                                yRange <- layer_scales(plt)$y$range$range
                                if(lMod$coefficients[2] > 0.3){
                                        xPos <- min(xRange) + abs(min(xRange) + 1.4)
                                        yPos <- max(yRange) - abs(max(yRange) * 0.01)
                                }else if(lMod$coefficients[2] <= 0.3 & lMod$coefficients[2] > 0){
                                        xPos <- max(xRange) - abs(max(xRange) * .6)
                                        yPos <- min(yRange) + abs(min(yRange) * 0.05)
                                }else{
                                        xPos <- max(xRange) - abs(max(xRange) * .6)
                                        yPos <- max(yRange) - abs(max(yRange) * 0.05)
                                }
                        }
                        plt <- plt +
                                geom_label(x = xPos,
                                           y = yPos,
                                           label = eq, parse = TRUE)
                }
                return(plt)
        }
}

# Load and parse data
################################################################################

metDat <- read.csv(metDatFile, row.names = 1)

dat <- read_table_fast(dataFile, row.names = 1)

if (preFiltGenes != "none"){
        preFilt <- read.csv(preFiltGenes, row.names = 1)
        keepGenes <- preFilt$ensembl_gene_id
        keepGenes <- keepGenes[!grepl("Intercept", keepGenes)]
        dat <- dat[, colnames(dat) %in% keepGenes]
}

# Data training/testing partition
################################################################################

# Split controls and test
# Split controls and the external dataset we will use as validation,
# (besides the 1/3 that of the controls that we are going to leave out of
# training)
dat_ctrls <- dat[rownames(dat) %in% make.names(metDat$specimenID[metDat$diagn_4BrainClck == "Control" & metDat$substudy != "brainSeq_pII" & metDat$substudy != "ageAnno"]), ]
dat_extrn <- dat[rownames(dat) %in% make.names(metDat$specimenID[metDat$diagn_4BrainClck == "Control" & metDat$substudy == "brainSeq_pII"]), ]
dat_rest <- dat[rownames(dat) %in% make.names(metDat$specimenID[metDat$diagn_4BrainClck == "Rest"]), ]

# Get a subdataset with individuals with a high braak index. We will use it 
# later to check the predictions.
dat_rest_highBraak <- dat_rest[rownames(dat_rest) %in% make.names(metDat$specimenID[metDat$Braak >= braakThrshld]), ]

# Obtain the ages of the controls and transform them
ageVec <- metDat$ageDeath[match(rownames(dat_ctrls),
                                make.names(metDat$specimenID))]
names(ageVec) <- rownames(dat_ctrls)

ageVecExtn <- metDat$ageDeath[match(rownames(dat_extrn),
                                    make.names(metDat$specimenID))]
names(ageVecExtn) <- rownames(dat_extrn)

# Split controls in training and testing and create the training DF,
# including the transformed ages
trainProp <- .66
sizeTrain <- round(nrow(dat_ctrls) * trainProp)
set.seed(123)
inTrain <- sample(rownames(dat_ctrls), size = sizeTrain)


# Do train/test splitting ensuring that all the 5 year intervals are represented
# at the desired proportion.
agesBins <- 20:100
agesBins <- agesBins[20:100 %% 5 == 0]

inTrain <- c()
seed <- 111
for(i in 1:(length(agesBins) - 1)){
        #i <- 1
        intLow <- agesBins[i]
        intHigh <- agesBins[i + 1]
        toSamp <- metDat$specimenID[metDat$ageDeath >= intLow & metDat$ageDeath < intHigh]
        toSamp <- toSamp[make.names(toSamp) %in% rownames(dat_ctrls)]
        sizeTrain <- round(length(toSamp) * trainProp)
        set.seed(seed)
        samped <- sample(toSamp, size = sizeTrain)
        inTrain <- c(inTrain, samped)
        seed <- seed + i
}

testSet_metDat <- metDat[make.names(metDat$specimenID) %in% rownames(dat_ctrls) & !(make.names(metDat$specimenID) %in% make.names(inTrain)), ]
extnSet_metDat <- metDat[make.names(metDat$specimenID) %in% rownames(dat_extrn), ]

write.csv(testSet_metDat, file = sprintf("%smetDat_testSet.csv", outDir))
write.csv(extnSet_metDat, file = sprintf("%smetDat_extnSet.csv", outDir))

dat_ctrls_train <- dat_ctrls[make.names(inTrain), ]
dat_ctrls_test <- dat_ctrls[!rownames(dat_ctrls) %in% make.names(inTrain), ]

ageChronTrain <- ageVec[match(rownames(dat_ctrls_train),
                              names(ageVec))]
        
train4Mod <- cbind.data.frame(ageChronTrain, dat_ctrls_train)


colnames(train4Mod)[1] <- respVar

# Model training
################################################################################

# Initialize h2o and fit the model
conn <- h2o.init(max_mem_size=mem)
trainData_h2o <- as.h2o(train4Mod)

if(lambdaFlag){
        lambda <- NULL
        lambda_search <- T
}else{
        lambda <- 0
        lambda_search <- F
}

model_dir <- sprintf("%smod_alpha%s", outDir, alph)
model_files <- list.files(model_dir)

if (length(model_files) == 0 | refit){
        # Train if model doesn't exist or if refit is selected to TRUE
        print("Training the model...")
        model <- h2o.glm(y = respVar,
                 training_frame = trainData_h2o,
                 nfolds = 10,
                 fold_assignment = "Random",
                 family = "AUTO",
                 link = "family_default",
                 lambda_search = lambda_search,
                 lambda = lambda,
                 standardize = T,
                 alpha = alph,
                 seed = 100,
                 max_active_predictors = ncol(trainData_h2o),
                 solver = "IRLSM")
        h2o.saveModel(model, path = sprintf("%smod_alpha%s", outDir, alph),
                      force = overwrite_mod)
}else{
        # If model exists and refit == FALSE, load latest model in the directory
        model_file <- sprintf("%s/%s",
                              model_dir,
                              model_files[length(model_files)])
        print(sprintf("Loading %s from %s...", basename(model_file), dirname(model_file)))
        model <- h2o.loadModel(model_file)
}

# Model evaluation in test dataset
################################################################################

# Plot metrics
cv_r2 <- unlist(as.vector(model@model$cross_validation_metrics_summary["r2", -c(1, 2)]))
cv_mae <- unlist(as.vector(model@model$cross_validation_metrics_summary["mae", -c(1, 2)]))
cv_rmse <- unlist(as.vector(model@model$cross_validation_metrics_summary["rmse", -c(1, 2)]))

# Obtain CV pearson2
#cv_preds <- h2o.cross_validation_predictions(model)
#actuals <- h2o.get_frame("your_training_frame_name")[["response_column_name"]]

plotDF <- data.frame(metric = c(rep("r2", length(cv_r2)),
                                rep("mae", length(cv_mae)),
                                rep("rmse", length(cv_rmse))),
                     value = c(cv_r2, cv_mae, cv_rmse))


mod_cvMetrics_plt <- ggplot(data = plotDF, mapping = aes(x = metric, y = value)) +
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

ggsave(filename = sprintf("%smod_alpha%s_cvMetrics.pdf", outDir, alph),
       plot = mod_cvMetrics_plt)

# Predict the age in the both test datasets (1:3 leave out and external dataset i.e.
# BrainSeq Phase II)
predicted <- predictAge(model, t(dat_ctrls_test))
predicted_extn <- predictAge(model, t(dat_extrn))

# Back-transform to chronological age
#if(respVar == "age_trans"){
#        predicted2Age <- back2Age(predicted)
#        testAges <- ageVec[match(names(predicted2Age), names(ageVec))]
#        testAgesTrns <- ageVecTrans[match(rownames(dat_ctrls_test),
#                                          names(ageVecTrans))]
#        
#        # Get metrics for test set
#        metrics_chronAge <- getMetrics(testAges, predicted2Age)
#        metrics_transAge <- getMetrics(testAgesTrns, predicted)
#        print("Metrics for test dataset on chronological ages:")
#        print(metrics_chronAge)
#        
#        print("Metrics for test dataset on transformed ages:")
#        print(metrics_transAge)
#       write.csv(metrics_chronAge,
#                 file = sprintf("%smetricsTest_chronAge_alpha%s.csv", outDir, alph))
#        write.csv(metrics_transAge,
#                  file = sprintf("%smetricsTest_transAge_alpha%s.csv", outDir, alph))
#        
#        # Plot R2 of training, CV and testing
#        dfr2 <- data.frame(r2_type = factor(c("r2_training",
#                                              rep("r2_cv",
#                                                  length(model@model$cross_validation_metrics_summary["r2", -c(1, 2)])),
#                                              "r2_test"),
#                                            levels = c("r2_training",
#                                                       "r2_cv",
#                                                       "r2_test")),
#                           value = c(h2o.r2(model, train = T),
#                                     unlist(model@model$cross_validation_metrics_summary["r2", -c(1, 2)]),
#                                     metrics_transAge$value[metrics_transAge$metric == "r2"]))
#}else if(respVar == "age_chron"){

# 1:3 leave out test dataset
testAges <- ageVec[match(names(predicted), names(ageVec))]
metrics_chronAge <- getMetrics(testAges, predicted)
print("Metrics for 1:3 leave out test dataset on chronological ages:")
print(metrics_chronAge)
write.csv(metrics_chronAge,
          file = sprintf("%smetricsTest_chronAge_alpha%s.csv",
                         outDir, alph))

# Obtain the residuals to check heteroskedasticity
test_data_preds <- data.frame(specimenID = names(testAges),
                              real_age = testAges,
                              pred_age = predicted)
test_realvspred_lm <- lm(formula = pred_age ~ real_age, data = test_data_preds)
test_data_preds$residuals <- residuals(test_realvspred_lm)
test_resid_plt <- ggplot(test_data_preds,
                         mapping = aes(x = real_age, y = residuals)) +
        geom_point()
ggsave(filename = sprintf("%stest_realvspred_resids.pdf", outDir),
       plot = test_resid_plt)
# Plot real age vs predicted age to see the relationship
tst_predage_plt <- ggplot(test_data_preds, mapping = aes(x = real_age, y = pred_age)) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_regline_equation()
ggsave(filename = sprintf("%stest_dat_predvsreal.pdf", outDir), plot = tst_predage_plt)



# External test dataset
metrics_chronAgeExtn <- getMetrics(ageVecExtn, predicted_extn)
print("Metrics for the external test dataset on chronological ages:")
print(metrics_chronAgeExtn)
write.csv(metrics_chronAgeExtn,
          file = sprintf("%smetricsExtn_chronAge_alpha%s.csv",
                         outDir, alph))

# Obtain the residuals to check heteroskedasticity
extn_data_preds <- data.frame(specimenID = names(ageVecExtn),
                              real_age = ageVecExtn,
                              pred_age = predicted_extn)
extn_realvspred_lm <- lm(formula = pred_age ~ real_age, data = extn_data_preds)
extn_data_preds$residuals <- residuals(extn_realvspred_lm)
extn_resid_plt <- ggplot(extn_data_preds, mapping = aes(x = real_age, y = residuals)) + geom_point()
ggsave(filename = sprintf("%sextn_realvspred_resids.pdf", outDir),
       plot = extn_resid_plt)
# Plot real age vs predicted age to see the relationship
ext_predage_plt <- ggplot(extn_data_preds, mapping = aes(x = real_age, y = pred_age)) +
        geom_point() +
        geom_smooth(method='lm') +
        stat_regline_equation()
ggsave(filename = sprintf("%sext_dat_predvsreal.pdf", outDir), plot = ext_predage_plt)

        # Plot R2 of training, CV and testing
        dfr2 <- data.frame(r2_type = factor(c("r2_training",
                                              rep("r2_cv",
                                                  length(model@model$cross_validation_metrics_summary["r2", -c(1, 2)])),
                                                  "r2_test_leave_out",
                                                  "r2_test_external"),
                                            levels = c("r2_training",
                                                       "r2_cv",
                                                       "r2_test_leave_out",
                                                       "r2_test_external")),
                           value = c(h2o.r2(model, train = T),
                                     unlist(model@model$cross_validation_metrics_summary["r2", -c(1, 2)]),
                                     metrics_chronAge$value[metrics_chronAge$metric == "r2"],
                                     metrics_chronAgeExtn$value[metrics_chronAgeExtn$metric == "r2"]))
        dfmae <- data.frame(mae_type = factor(c("mae_training",
                                              rep("mae_cv",
                                                  length(model@model$cross_validation_metrics_summary["mae", -c(1, 2)])),
                                                  "mae_test_leave_out",
                                                  "mae_test_external"),
                                            levels = c("mae_training",
                                                       "mae_cv",
                                                       "mae_test_leave_out",
                                                       "mae_test_external")),
                           value = c(h2o.mae(model, train = T),
                                     unlist(model@model$cross_validation_metrics_summary["mae", -c(1, 2)]),
                                     metrics_chronAge$value[metrics_chronAge$metric == "mae"],
                                     metrics_chronAgeExtn$value[metrics_chronAgeExtn$metric == "mae"]))
#}


dfr2_bxplt <- dfr2[dfr2$r2_type == "r2_cv", ]
dfr2_bxplt$r2_type <- factor(dfr2_bxplt$r2_type, levels = c("r2_training",
                                                            "r2_cv",
                                                            "r2_test_leave_out",
                                                            "r2_test_external"))
dfr2_point <- dfr2[dfr2$r2_type != "r2_cv", ]
dfr2_point$r2_type <- factor(dfr2_point$r2_type, levels = c("r2_training",
                                                            "r2_cv",
                                                            "r2_test_leave_out",
                                                            "r2_test_external"))

r2_plt <- ggplot(data = dfr2, mapping = aes(x = r2_type, y = value)) +
        geom_point(data = dfr2_point) + ylim(0, 1) +
        geom_boxplot(data = dfr2_bxplt, outlier.shape = NA) +
        geom_jitter(data = dfr2_bxplt) +
        scale_x_discrete(limits = c("r2_training",
                                    "r2_cv",
                                    "r2_test_leave_out",
                                    "r2_test_external"),
                         labels = c("r2_training" = "train set",
                                    "r2_cv" = "cross-validation",
                                    "r2_test_leave_out" = "test set (1:3 leave-out)",
                                    "r2_test_external" = "test set (external)")) +
        ylab(expression(R^2)) +
        ylim(c(floor(min(dfr2$value * 10)) / 10, 1)) + 
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15, angle = 90, hjust = 1),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, linewidth = 1))

ggsave(sprintf("%smodelAlph%sR2.pdf", outDir, as.character(alph)),
       height = 5, width = 2,
       plot = r2_plt)

dfmae_bxplt <- dfmae[dfmae$mae_type == "mae_cv", ]
dfmae_bxplt$mae_type <- factor(dfmae_bxplt$mae_type, levels = c("mae_training",
                                                                "mae_cv",
                                                                "mae_test_leave_out",
                                                                "mae_test_external"))
dfmae_point <- dfmae[dfmae$mae_type != "mae_cv", ]
dfmae_point$mae_type <- factor(dfmae_point$mae_type, levels = c("mae_training",
                                                                "mae_cv",
                                                                "mae_test_leave_out",
                                                                "mae_test_external"))

mae_plt <- ggplot(data = dfmae, mapping = aes(x = mae_type, y = value)) +
        geom_point(data = dfmae_point) + ylim(0, 1) +
        geom_boxplot(data = dfmae_bxplt, outlier.shape = NA) +
        geom_jitter(data = dfmae_bxplt) +
        scale_x_discrete(limits = c("mae_training",
                                    "mae_cv",
                                    "mae_test_leave_out",
                                    "mae_test_external"),
                         labels = c("mae_training" = "train set",
                                    "mae_cv" = "cross-validation",
                                    "mae_test" = "test set (1:3 leave-out)",
                                    "mae_test" = "test set (external)")) +
        ylab("MAE") +
        ylim(c(floor(min(dfmae$value)), ceiling(max(dfmae$value)))) + 
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15, angle = 90, hjust = 1),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, linewidth = 1))

ggsave(sprintf("%smodelAlph%sMAE.pdf", outDir, as.character(alph)),
       height = 5, width = 2)

# Assessment in neurodegenerated individuals
################################################################################

# Predict age in neurodegenerated individuals
predicted_hBraak <- predictAge(model, t(dat_rest_highBraak))

agesHighBRaak <- metDat$ageDeath[match(rownames(dat_rest_highBraak), make.names(metDat$specimenID))]
names(agesHighBRaak) <- rownames(dat_rest_highBraak)



#if(respVar == "age_trans"){
#        # Transform actual chronological ages to transformed ages to compute
#       # R2
#        agesHighBRaakTrns <- test(agesHighBRaak)
#       # Back-transform to chronological age
#        predicted_hBraak2Age <- back2Age(predicted_hBraak)
#        # Get metrics for the high braak 
#        metrics_hBraak_chronAge <- getMetrics(agesHighBRaak, predicted_hBraak2Age)
#        metrics_hBraak_transAge <- getMetrics(agesHighBRaakTrns, predicted_hBraak)
#        
#        print("Metrics for high-Braak samples on chronological ages:")
#        print(metrics_hBraak_chronAge)
#        
#        print("Metrics for high-Braak samples on transformed ages:")
#        print(metrics_hBraak_transAge)
#        
#        write.csv(metrics_hBraak_chronAge,
#                  file = sprintf("%smetricsHBraak_chronAge_alpha%s.csv", outDir, alph))
#        write.csv(metrics_hBraak_transAge,
#                  file = sprintf("%smetricsHBraak_transAge_alpha%s.csv", outDir, alph))
        
#        # Plot chronological age vs predicted age in test control dataset, and 
#        # in ND dataset. Consider only ages from the minimum neurodegenerated onwards.
#        comp_predAges_test_ND(testAges,
#                              agesHighBRaak,
#                              predicted2Age,
#                              predicted_hBraak2Age,
#                              color = "group")$plot
#        
#        ggsave(sprintf("%sctrls_highBraak_diagn_alpha%s.pdf", outDir, alph), 
#               height = 8,
#               width = 9)
#        
#        comp_predAges_test_ND(testAges,
#                              agesHighBRaak,
#                              predicted2Age,
#                              predicted_hBraak2Age,
#                              color = "braak")$plot
#        
#        ggsave(sprintf("%sctrls_highBraak_braak_alpha%s.pdf", outDir, alph),
#               height = 8, width = 9)
#       
#        comp_predAges_test_ND(testAges,
#                              agesHighBRaak,
#                              predicted2Age,
#                              predicted_hBraak2Age,
#                              color = "substudy")$plot
#        
#        ggsave(sprintf("%sctrls_highBraak_batch_alpha%s.pdf", outDir, alph),
#               height = 8, width = 9)
#        
#        
#        # Plot chronological age vs predicted age in test control dataset, and 
#        # in ND dataset. Consider all ages..
#        predAgeDF_tst_vs_nd <- comp_predAges_test_ND(testAges,
#                                                     agesHighBRaak,
#                                                     predicted2Age,
#                                                     predicted_hBraak2Age,
#                                                     color = "group",
#                                                     filtCommAgeRank = F)
#        
#        predAgeDF_tst_vs_nd$plot
#        
#        ggsave(sprintf("%sctrls_highBraak_all_diagn_alpha%s.pdf", outDir, alph),
#               height = 8, width = 9)
#        
#        comp_predAges_test_ND(testAges,
#                              agesHighBRaak,
#                              predicted2Age,
#                              predicted_hBraak2Age,
#                              color = "braak",
#                              filtCommAgeRank = F)$plot
#        
#        ggsave(sprintf("%sctrls_highBraak_all_braak_alpha%s.pdf", outDir, alph),
#               height = 8, width = 9)
#        
#        comp_predAges_test_ND(testAges,
#                              agesHighBRaak,
#                              predicted2Age,
#                              predicted_hBraak2Age,
#                              color = "substudy",
#                              filtCommAgeRank = F)$plot
#        
#        ggsave(sprintf("%sctrls_highBraak_all_batch_alpha%s.pdf", outDir, alph),
#               height = 8, width = 9)
#        
#}else if(respVar == "age_chron"){
        metrics_hBraak_chronAge <- getMetrics(agesHighBRaak, predicted_hBraak)
        write.csv(metrics_hBraak_chronAge,
                  file = sprintf("%smetricsHBraak_chronAge_alpha%s.csv", outDir, alph))
        
        # Plot chronological age vs predicted age in test control dataset, and 
        # in ND dataset. Consider only ages from the minimum neurodegenerated onwards.
        comp_predAges_test_ND(testAges,
                              agesHighBRaak,
                              predicted,
                              predicted_hBraak,
                              color = "group")$plot
        
        ggsave(sprintf("%sctrls_highBraak_diagn_alpha%s.pdf", outDir, alph), 
               height = 8,
               width = 9)
        
        ctrls_hBraak_plt <- comp_predAges_test_ND(testAges,
                              agesHighBRaak,
                              predicted,
                              predicted_hBraak,
                              color = "braak")$plot
        
        ggsave(sprintf("%sctrls_highBraak_braak_alpha%s.pdf", outDir, alph),
               height = 8, width = 9, plot = ctrls_hBraak_plt)
        
        ctrls_hBraak_batchCol_plt <- comp_predAges_test_ND(testAges,
                              agesHighBRaak,
                              predicted,
                              predicted_hBraak,
                              color = "substudy")$plot
        
        ggsave(sprintf("%sctrls_highBraak_batch_alpha%s.pdf", outDir, alph),
               height = 8, width = 9,
               plot = ctrls_hBraak_batchCol_plt)
        
        
        # Plot chronological age vs predicted age in test control dataset, and 
        # in ND dataset. Consider all ages..
        predAgeDF_tst_vs_nd <- comp_predAges_test_ND(testAges,
                                                     agesHighBRaak,
                                                     predicted,
                                                     predicted_hBraak,
                                                     color = "group",
                                                     filtCommAgeRank = F)
        
        cntrls_hBraak_allDiagn_plt <- predAgeDF_tst_vs_nd$plot
        
        ggsave(sprintf("%sctrls_highBraak_all_diagn_alpha%s.pdf", outDir, alph),
               height = 8, width = 9,
               plot = cntrls_hBraak_allDiagn_plt)
        
        cntrls_hBraak_allBraak_plt <- comp_predAges_test_ND(testAges,
                              agesHighBRaak,
                              predicted,
                              predicted_hBraak,
                              color = "braak",
                              filtCommAgeRank = F)$plot
        
        ggsave(sprintf("%sctrls_highBraak_all_braak_alpha%s.pdf", outDir, alph),
               height = 8, width = 9,
               plot = cntrls_hBraak_allBraak_plt)
        
        cntrls_hBraak_allBatch_plt <- comp_predAges_test_ND(testAges,
                              agesHighBRaak,
                              predicted,
                              predicted_hBraak,
                              color = "substudy",
                              filtCommAgeRank = F)$plot
        
        ggsave(sprintf("%sctrls_highBraak_all_batch_alpha%s.pdf", outDir, alph),
               height = 8, width = 9,
               plot = cntrls_hBraak_allBatch_plt)
#}



# Do Mann-Whitney test between high-braak and controls in the 60 to 70 age
# range
df4DisContComp_60to70 <- predAgeDF_tst_vs_nd$data[predAgeDF_tst_vs_nd$data$age_death >= 60 & predAgeDF_tst_vs_nd$data$age_death <= 70, ]


#plot(density(df4DisContComp_60to70$pred_age[df4DisContComp_60to70$group == "nd"]))
#plot(density(df4DisContComp_60to70$pred_age[df4DisContComp_60to70$group == "control"]))

dis_vs_cont_60to70bxplt <-ggplot(df4DisContComp_60to70, mapping = aes(x = group, y = pred_age)) +
        geom_boxplot() +
        stat_compare_means(method = "wilcox.test") +
        stat_compare_means(method = "wilcox.test") +
        ylab("predicted age") +
        
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              axis.title.y = element_text(size=20),
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              legend.text = element_text(size=12),
              legend.title = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, linewidth = 1))

ggsave(filename = sprintf("%spredAge60to70ctrlVsHBraak.pdf", outDir),
       width = 4,
       height = 5,
       plot = dis_vs_cont_60to70bxplt)

# Do a LM and an ANCOVA to assess if the influence of neurodegeneration on
# predicted age is significant.

predAgeDF_tst_vs_nd_ancova <- predAgeDF_tst_vs_nd$data[predAgeDF_tst_vs_nd$data$age_death >= min(predAgeDF_tst_vs_nd$data$age_death[predAgeDF_tst_vs_nd$data$group == "nd"]), ]

mod_4Anc <- lm(pred_age ~ age_death * group, data = predAgeDF_tst_vs_nd_ancova)

print("LM")
print(summary(mod_4Anc))
capture.output(summary(mod_4Anc),
               file = sprintf("%sND_lm_pred_vs_chronND_alph%s.txt",
                              outDir, as.character(alph)))

ancova_mod <- aov(pred_age ~ age_death + group + age_death:group,
                  data = predAgeDF_tst_vs_nd_ancova)

print("ANCOVA")
print(summary(ancova_mod))
capture.output(summary(ancova_mod),
               file = sprintf("%sND_ancova_alph%s.txt",
                              outDir, as.character(alph)))

ancova_plot <- ggplot(predAgeDF_tst_vs_nd_ancova, aes(x = age_death, y = pred_age, color = group)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        labs(x = "chronological age", y = "predicted age") +
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
ggsave(sprintf("%sctrls_highBraak_diagn_alpha%s_wLines.pdf", outDir, alph), 
       height = 8,
       width = 9,
       ancova_plot)

ancova_mod_2 <- rstatix::anova_test(data = predAgeDF_tst_vs_nd_ancova,
                                    formula = pred_age ~ age_death + group + age_death:group)

print(ancova_mod_2)
capture.output(ancova_mod_2,
               file = sprintf("%sND_ancovaFitRstatix_alph%s.txt",
                              outDir, as.character(alph)))

# Obtain boxplots of predicted age in each braak stage in each decade, from
# 60 to 100

# Predict ages in all the samples that are in rest and have braak index
dat_rest_wBraak <- dat_rest[rownames(dat_rest) %in% make.names(metDat$specimenID[!is.na(metDat$Braak)]), ]
wBraak_ages <- metDat$ageDeath[match(rownames(dat_rest_wBraak), make.names(metDat$specimenID))]
predicted_wBraak <- predictAge(model, t(dat_rest_wBraak))


#if(respVar == "age_trans"){
#        predicted_wBraak2Age <- back2Age(predicted_wBraak)
#        
#        dfBraak <- data.frame(specimenID = metDat$specimenID[match(names(predicted_wBraak2Age),
#                                                                   make.names(metDat$specimenID))],
#                              age_death = wBraak_ages,
#                              pred_age = predicted_wBraak2Age,
#                              braak = metDat$Braak[match(names(predicted_wBraak2Age),
#                                                         make.names(metDat$specimenID))],
#                              group = "Rest",
#                              substudy = metDat$substudy[match(names(predicted_wBraak2Age),
#                                                               make.names(metDat$specimenID))])
#}else if (respVar == "age_chron"){
        dfBraak <- data.frame(specimenID = metDat$specimenID[match(names(predicted_wBraak),
                                                                   make.names(metDat$specimenID))],
                              age_death = wBraak_ages,
                              pred_age = predicted_wBraak,
                              braak = metDat$Braak[match(names(predicted_wBraak),
                                                         make.names(metDat$specimenID))],
                              group = "Rest",
                              substudy = metDat$substudy[match(names(predicted_wBraak),
                                                               make.names(metDat$specimenID))])
#}

df4DisContComp <- predAgeDF_tst_vs_nd$data
df4DisContComp$braak <- rep(NA, nrow(df4DisContComp))
dfBraak <- rbind.data.frame(dfBraak,
                            df4DisContComp[df4DisContComp$group == "control",
                                           c("specimenID", "age_death", "pred_age", "braak", "group", "substudy")])

dfBraak$ageInt <- rep(NA, nrow(dfBraak))
ageInt <- 20:100
ageInt <- ageInt[ageInt %% 10 == 0]

for(i in 1:(length(ageInt) - 1)){
        age_1 <- ageInt[i]
        age_2 <- ageInt[i + 1]
        dfBraak$ageInt[dfBraak$age_death >= age_1 & dfBraak$age_death < age_2] <- sprintf("%s_%s", age_1, age_2)
        
}

dfBraak$braak[dfBraak$group == "control"] <- "control"


age60_70_batch <- plotAgeInt(dfBraak = dfBraak,
                             ageInt = "60_70",
                             color = "substudy")
age70_80_batch <- plotAgeInt(dfBraak = dfBraak,
                             ageInt = "70_80",
                             color = "substudy")
age80_90_batch <- plotAgeInt(dfBraak = dfBraak,
                             ageInt = "80_90",
                             color = "substudy")
age90_100_batch <- plotAgeInt(dfBraak = dfBraak,
                              ageInt = "90_100",
                              color = "substudy")

ggarrange(plotlist = list(age60_70_batch,
                          age70_80_batch,
                          age80_90_batch,
                          age90_100_batch))
ggsave(filename = sprintf("%spredAg_byBraakAndAge_batch_alph%s.pdf",
                          outDir, alph),
       width = 10,
       height = 6)

age60_70 <- plotAgeInt(dfBraak = dfBraak,
                       ageInt = "60_70")
age70_80 <- plotAgeInt(dfBraak = dfBraak,
                       ageInt = "70_80")
age80_90 <- plotAgeInt(dfBraak = dfBraak,
                       ageInt = "80_90")
age90_100 <- plotAgeInt(dfBraak = dfBraak,
                        ageInt = "90_100")

ggarrange(plotlist = list(age60_70,
                          age70_80,
                          age80_90,
                          age90_100))
ggsave(filename = sprintf("%spredAg_byBraakAndAge_alph%s.pdf",
                          outDir, alph),
       width = 10,
       height = 6)

# Assesss if the relationship between braak index and predicted age is
# significant
getModPredAgeVsBraak(dfBraak,
                     outName = sprintf("%sND_lmFit_alph%s_all.txt",
                                       outDir, as.character(alph)),
                     savePlt = T,
                     plotEq = T)

nd_lmPlt_60_70 <- getModPredAgeVsBraak(dfBraak,
                                       age_int = "60_70",
                                       outName = sprintf("%sND_lmFit_alph%s_60_70.txt",
                                                         outDir,
                                                         as.character(alph)),
                                       plotEq = T)

nd_lmPlt_70_80 <- getModPredAgeVsBraak(dfBraak,
                                       age_int = "70_80",
                                       outName = sprintf("%sND_lmFit_alph%s_70_80.txt",
                                                         outDir,
                                                         as.character(alph)),
                                       plotEq = T)

nd_lmPlt_80_90 <- getModPredAgeVsBraak(dfBraak,
                                       age_int = "80_90",
                                       outName = sprintf("%sND_lmFit_alph%s_80_90.txt",
                                                         outDir,
                                                         as.character(alph)),
                                       plotEq = T)

nd_lmPlt_90_100 <- getModPredAgeVsBraak(dfBraak,
                                        age_int = "90_100",
                                        outName = sprintf("%sND_lmFit_alph%s_90_100.txt",
                                                          outDir,
                                                          as.character(alph)),
                                        plotEq = T)
ggarrange(plotlist = list(nd_lmPlt_60_70,
                          nd_lmPlt_70_80,
                          nd_lmPlt_80_90,
                          nd_lmPlt_90_100))

ggsave(filename = sprintf("%sND_lmFit_alph%s_byAge.pdf",
                          outDir, as.character(alph)),
       width = 12,
       height = 10)

# Assess relationship between braak index and predicted age, accounting
# for chronoligical age (pred age ~ braak + chron_age)

nd_lmPlt_all_accCrhonAge <- getModPredAgeVsBraak(dfBraak,
                                                 outName = sprintf("%sND_lmFit_alph%s_all_accChronAge.txt",
                                                                   outDir, as.character(alph)),
                                                 savePlt = T,
                                                 plotEq = T,
                                                 accountChronAge = T)

ggsave(sprintf("%sND_lmFit_alph%s_all_accChronAge.pdf",
               outDir, as.character(alph)),
       plot = nd_lmPlt_all_accCrhonAge,
       height = 8,
       width = 9)

nd_lmPlt_60_70_accCrhonAge <- getModPredAgeVsBraak(dfBraak,
                                                   age_int = "60_70",
                                                   outName = sprintf("%sND_lmFit_alph%s_60_70_accChronAge.txt",
                                                                     outDir,
                                                                     as.character(alph)),
                                                   plotEq = T,
                                                   accountChronAge = T)
nd_lmPlt_70_80_accCrhonAge <- getModPredAgeVsBraak(dfBraak,
                                                   age_int = "70_80",
                                                   outName = sprintf("%sND_lmFit_alph%s_70_80_accChronAge.txt",
                                                                     outDir,
                                                                     as.character(alph)),
                                                   plotEq = T,
                                                   accountChronAge = T)
nd_lmPlt_80_90_accCrhonAge <- getModPredAgeVsBraak(dfBraak,
                                                   age_int = "80_90",
                                                   outName = sprintf("%sND_lmFit_alph%s_80_90_accChronAge.txt",
                                                                     outDir,
                                                                     as.character(alph)),
                                                   plotEq = T,
                                                   accountChronAge = T)
nd_lmPlt_90_100_accCrhonAge <- getModPredAgeVsBraak(dfBraak,
                                                    age_int = "90_100",
                                                    outName = sprintf("%sND_lmFit_alph%s_90_100_accChronAge.txt",
                                                                      outDir,
                                                                      as.character(alph)),
                                                    plotEq = T,
                                                    accountChronAge = T)


ggarrange(plotlist = list(nd_lmPlt_60_70_accCrhonAge,
                          nd_lmPlt_70_80_accCrhonAge,
                          nd_lmPlt_80_90_accCrhonAge,
                          nd_lmPlt_90_100_accCrhonAge))

ggsave(filename = sprintf("%sND_lmFit_alph%s_byAge_accCrhonAge.pdf",
                          outDir, as.character(alph)),
       width = 12,
       height = 10)

h2o.shutdown(prompt = F)