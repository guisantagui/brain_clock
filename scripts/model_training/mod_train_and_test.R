################################################################################
#                                                                              #
# Brain Clock: fit a GLM model to the processed brain expression dataset to    #
# predict age, and assess if there are differences in predicted ages in the    #
# neurodegenerated subjects                                                    #
#                                                                              #
################################################################################
if(!require(limma, quietly = T)) BiocManager::install("limma", update = F)
if(!require("h2o", quietly = T)){
        install.packages("h2o",
                         type="source",
                         repos="https://h2o-release.s3.amazonaws.com/h2o/rel-3.46.0/2/R")
}
library(h2o)
if(!require("ggplot2")) install.packages("ggplot2", repos='https://pbil.univ-lyon1.fr/CRAN/')
library(ggplot2)
if(!require("ggpubr")) install.packages("ggpubr", repos='https://pbil.univ-lyon1.fr/CRAN/')
library(ggpubr)

# Directory stuff
################################################################################

geneInfoFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_geneInfo_wTBI.csv" # gene information file
metDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat_wTBI.csv" # Integrated metadata file
ageTransParFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/GompertzMakehamParameters.rds" # rds object with the age transformation parameters
dataFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/combined_TPMs_wTBI_preproc_log10_noCerebell_combat_pmiAgeAdj_onlyAge_svaAdj.csv" # Preprocessed RNAseq dataset
dataFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/combined_TPMs_wTBI_preproc_log10_noCerebell_onlyAge_svaAdj.csv" # Preprocessed RNAseq dataset
#dataFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/combined_postMergeTPMs_preproc_log10_noCerebell_combat_pmiAgeAdj_onlyAge_svaAdj.csv"
modelFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/wTBI_allGenes/modFuncsAlpha0.5.rds"
badamModFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/badam_procs_wTBI/modFuncsAlpha0.5.rds"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/modAllGenes_modBadamSign_wTBI/" # output directory.
alph <- 0.5 # Alpha parameter for the glm


if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Functions for doing age transformation
################################################################################

vals <- readRDS(ageTransParFile)
vals <- round(vals,6) #Makes the parameters more readable with little loss of accuracy
test <- function(x){
        exp(-vals[3]*x-(vals[1]/vals[2])*(exp(vals[2]*x)-1))
}

inverse = function (f, lower = 0, upper = 110) {
        function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}
agetrans_inverse = function(lower,upper) inverse(test, lower, upper)


# Load and parse data
################################################################################
badamMod <- readRDS(badamModFile)
modOld <- readRDS(modelFile)
geneInfo <- read.csv(geneInfoFile, row.names = 1)
metDat <- read.csv(metDatFile, row.names = 1)

dat <- data.frame(data.table::fread(dataFile))
rownames(dat) <- dat$V1
dat <- dat[, colnames(dat) != "V1"]

# Change names of genes to symbols and remove genes that don't have a symbol.
symbNames <- geneInfo$symbol[match(colnames(dat), geneInfo$ensembl_gene_id_orig)]

dat <- dat[, symbNames != ""]

symbNames <- symbNames[symbNames != ""]

colnames(dat) <- symbNames

# Get list of genes from old model that don't have a coefficient equal to zero
modSign <- modOld@model$coefficients_table$names[modOld@model$coefficients_table$coefficients != 0]
modSign <- modSign[modSign != "Intercept"]

# Get list of genes in badam model
badamModSign <- badamMod@model$coefficients_table$names[badamMod@model$coefficients_table$coefficients != 0]

badamModSign <- badamModSign[make.names(badamModSign) %in% colnames(dat)]

intersect(badamModSign, modSign)

# Keep only the genes that are either significant between impaired and 
# not impaired, or that are hallmark genes in aging
dat <- dat[, colnames(dat) %in% c(badamModSign, modSign)]

# Analysis
################################################################################

# Split controls and rest
dat_ctrls <- dat[rownames(dat) %in% make.names(metDat$specimenID[metDat$diagn_4BrainClck == "Control"]), ]
dat_rest <- dat[rownames(dat) %in% make.names(metDat$specimenID[metDat$diagn_4BrainClck == "Rest"]), ]

# Get a subdataset with individuals with a high braak index. We will use it 
# later to check the predictions.
dat_rest_highBraak <- dat_rest[rownames(dat_rest) %in% make.names(metDat$specimenID[metDat$Braak >= 4]), ]

# Obtain the ages of the controls and transform them
ageVec <- metDat$ageDeath[match(rownames(dat_ctrls),
                                make.names(metDat$specimenID))]

names(ageVec) <- rownames(dat_ctrls)

ageVecTrans <- test(ageVec)

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

dat_ctrls_train <- dat_ctrls[make.names(inTrain), ]
dat_ctrls_test <- dat_ctrls[!rownames(dat_ctrls) %in% make.names(inTrain), ]

ageTransTrain <- ageVecTrans[match(rownames(dat_ctrls_train),
                                   names(ageVecTrans))]

train4Mod <- cbind.data.frame(ageTransTrain, dat_ctrls_train)

colnames(train4Mod)[1] <- "age_trans"

# Initialize h2o and fit the model
conn <- h2o.init(max_mem_size="16G")
trainData_h2o <- as.h2o(train4Mod)

lambda <- NULL
lambda_search <- T
#alph <- 0.5
#alph <- 0

model <- h2o.glm(y = "age_trans",
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

#sum(procGenes %in% model@model$coefficients_table$names[model@model$coefficients_table$coefficients != 0])

# Plot metrics
cv_r2 <- unlist(as.vector(model@model$cross_validation_metrics_summary["r2", -c(1, 2)]))
cv_mae <- unlist(as.vector(model@model$cross_validation_metrics_summary["mae", -c(1, 2)]))
cv_rmse <- unlist(as.vector(model@model$cross_validation_metrics_summary["rmse", -c(1, 2)]))

plotDF <- data.frame(metric = c(rep("r2", length(cv_r2)),
                                rep("mae", length(cv_mae)),
                                rep("rmse", length(cv_rmse))),
                     value = c(cv_r2, cv_mae, cv_rmse))


ggplot(data = plotDF, mapping = aes(x = metric, y = value)) +
        geom_boxplot(outlier.shape = NA) + geom_jitter()


ggsave(filename = sprintf("%smodFuncsAlpha%s_cvMetrics.pdf", outDir, alph))

# Predict the age in the test dataset
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

predicted <- predictAge(model, t(dat_ctrls_test))

# Clip values higher than 1 to .99 to avoid errors, and lower than 0 to
# 0.01
predicted[predicted > 1] <- max(train4Mod$age_trans)
predicted[predicted < 0] <- .000012

# Re-convert predicted back to age
predicted2Age <- sapply(predicted, function(x) agetrans_inverse(0, 110)(x))

predicted2Age <- unlist(predicted2Age)
names(predicted2Age) <- gsub(".root", "", names(predicted2Age))

testAges <- ageVec[match(names(predicted2Age), names(ageVec))]

# Get metrics
mae_test <- sum(abs(testAges - predicted2Age))/length(testAges)
mse_test <- sum((testAges - predicted2Age)^2)/length(testAges)
rmse_test <- sqrt(mse_test)
mape_test <- (sum(abs((testAges - predicted2Age)/testAges))/length(testAges))*100
r2_test <- 1 - sum((testAges - predicted2Age)^2)/sum((testAges - mean(testAges))^2)

print(sprintf("MAE in test dataset: %s", round(mae_test, digits = 3)))
print(sprintf("MSE in test dataset: %s", round(mse_test, digits = 3)))
print(sprintf("RMSE in test dataset: %s", round(rmse_test, digits = 3)))
print(sprintf("MAPE in test dataset: %s", round(mape_test, digits = 3)))
print(sprintf("R2 in test dataset: %s", round(r2_test, digits = 3)))

# Plot R2 of training, CV and testing
dfr2 <- data.frame(r2_type = factor(c("r2_training",
                                      rep("r2_cv",
                                          length(model@model$cross_validation_metrics_summary["r2", -c(1, 2)])), "r2_test"),
                                    levels = c("r2_training",
                                               "r2_cv",
                                               "r2_test")),
                   value = c(h2o.r2(model, train = T),
                             unlist(model@model$cross_validation_metrics_summary["r2", -c(1, 2)]),
                             r2_test))

dfr2_bxplt <- dfr2[dfr2$r2_type == "r2_cv", ]
dfr2_bxplt$r2_type <- factor(dfr2_bxplt$r2_type, levels = c("r2_training",
                                                            "r2_cv",
                                                            "r2_test"))
dfr2_point <- dfr2[dfr2$r2_type != "r2_cv", ]
dfr2_point$r2_type <- factor(dfr2_point$r2_type, levels = c("r2_training",
                                                            "r2_cv",
                                                            "r2_test"))
dfr2_point$r2_type

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
                                          fill=NA, size=1))

ggsave(sprintf("%smodelAlph%sR2.pdf",
               outDir, as.character(alph)),
       height = 3, width = 4)

metrics_test <- data.frame(metric = c("mae", "mse", "rmse", "mape", "r2"),
                           value = c(mae_test,
                                     mse_test,
                                     rmse_test,
                                     mape_test,
                                     r2_test))

write.csv(metrics_test,
          file = sprintf("%smetricsTest_alpha%s.csv", outDir, alph))

# Predict age in neurodegenerated individuals
predicted_hBraak <- predictAge(model, t(dat_rest_highBraak))
# Clip values higher than 1 to .99 to avoid errors, and lower than 0 to
# .000012
predicted_hBraak[predicted_hBraak > 1] <- max(train4Mod$age_trans)
predicted_hBraak[predicted_hBraak < 0] <- .000012


predicted_hBraak2Age <- sapply(predicted_hBraak, function(x) agetrans_inverse(0, 110)(x))
predicted_hBraak2Age <- unlist(predicted_hBraak2Age)
names(predicted_hBraak2Age) <- gsub(".root", "", names(predicted_hBraak2Age))



agesHighBRaak <- metDat$ageDeath[match(rownames(dat_rest_highBraak),
                                       make.names(metDat$specimenID))]
names(agesHighBRaak) <- rownames(dat_rest_highBraak)

# Get metrics for the high braak 
mae_hBraak <- sum(abs(agesHighBRaak - predicted_hBraak2Age))/length(agesHighBRaak)
mse_hBraak <- sum((agesHighBRaak - predicted_hBraak2Age)^2)/length(agesHighBRaak)
rmse_hBraak <- sqrt(mse_hBraak)
mape_hBraak <- (sum(abs((agesHighBRaak - predicted_hBraak2Age)/agesHighBRaak))/length(agesHighBRaak))*100
r2_hBraak <- 1 - sum((agesHighBRaak - predicted_hBraak2Age)^2)/sum((agesHighBRaak - mean(agesHighBRaak))^2)

print(sprintf("MAE in high-braak dataset: %s", round(mae_hBraak, digits = 3)))
print(sprintf("MSE in high-braak dataset: %s", round(mse_hBraak, digits = 3)))
print(sprintf("RMSE in high-braak dataset: %s", round(rmse_hBraak, digits = 3)))
print(sprintf("MAPE in high-braak dataset: %s", round(mape_hBraak, digits = 3)))
print(sprintf("R2 in high-braak dataset: %s", round(r2_hBraak, digits = 3)))

metrics_hBraak <- data.frame(metric = c("mae", "mse", "rmse", "mape", "r2"),
                             value = c(mae_hBraak,
                                       mse_hBraak,
                                       rmse_hBraak,
                                       mape_hBraak,
                                       r2_hBraak))

write.csv(metrics_hBraak,
          file = sprintf("%smetricsHBraak_alpha%s.csv", outDir, alph))

# Plot ages vs predicted ages in test control dataset and in the high braak
# dataset
testAges_subsetHigh <- testAges[testAges >= min(as.numeric(names(table(agesHighBRaak)))) & testAges <= max(as.numeric(names(table(agesHighBRaak))))]
predicted2Age_subsetHigh <- predicted2Age[testAges >= min(as.numeric(names(table(agesHighBRaak)))) & testAges <= max(as.numeric(names(table(agesHighBRaak))))]

df4DisContComp <- data.frame(age_death = c(testAges_subsetHigh,
                                           agesHighBRaak),
                             specimenID = names(c(testAges_subsetHigh,
                                                  agesHighBRaak)),
                             pred_age = c(predicted2Age_subsetHigh,
                                          predicted_hBraak2Age),
                             group = c(rep("control", length(testAges_subsetHigh)),
                                       rep("high_braak", length(agesHighBRaak))))

df4DisContComp$braak <- metDat$Braak[match(df4DisContComp$specimenID,
                                           metDat$specimenID)]

df4DisContComp$substudy <- metDat$substudy[match(df4DisContComp$specimenID,
                                                 make.names(metDat$specimenID))]


ggplot(data = df4DisContComp,
       mapping = aes(x = age_death, y = pred_age, col = group)) +
        geom_point()

ggsave(sprintf("%sctrls_highBraak_diagn_alpha%s.pdf", outDir, alph))


ggplot(data = df4DisContComp,
       mapping = aes(x = age_death, y = pred_age, col = braak)) +
        geom_point()

ggsave(sprintf("%sctrls_highBraak_braak_alpha%s.pdf", outDir, alph))


ggplot(data = df4DisContComp,
       mapping = aes(x = age_death, y = pred_age, col = substudy)) +
        geom_point()
ggsave(sprintf("%sctrls_highBraak_batch_alpha%s.pdf", outDir, alph))


df4DisContComp <- data.frame(age_death = c(testAges,
                                           agesHighBRaak),
                             specimenID = names(c(testAges,
                                                  agesHighBRaak)),
                             pred_age = c(predicted2Age,
                                          predicted_hBraak2Age),
                             group = c(rep("control", length(testAges)),
                                       rep("high_braak", length(agesHighBRaak))))

df4DisContComp$group <- gsub("high_braak",
                             "neurodegenerated",
                             df4DisContComp$group)

df4DisContComp$substudy <- metDat$substudy[match(df4DisContComp$specimenID,
                                                 make.names(metDat$specimenID))]

ggplot(data = df4DisContComp,
       mapping = aes(x = age_death, y = pred_age, col = group)) +
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
                                          fill=NA, size=1))

ggsave(sprintf("%sctrls_highBraak_all_diagn_alpha%s.pdf", outDir, alph), height = 8, width = 9)

ggplot(data = df4DisContComp,
       mapping = aes(x = age_death, y = pred_age, col = substudy)) +
        geom_point()

ggsave(sprintf("%sctrls_highBraak_all_batch_alpha%s.pdf", outDir, alph))


df4DisContComp_60to70 <- df4DisContComp[df4DisContComp$age_death >= 60 & df4DisContComp$age_death <= 70, ]


plot(density(df4DisContComp_60to70$pred_age[df4DisContComp_60to70$group == "high_braak"]))
plot(density(df4DisContComp_60to70$pred_age[df4DisContComp_60to70$group == "control"]))

ggplot(df4DisContComp_60to70, mapping = aes(x = group, y = pred_age)) +
        geom_boxplot() +
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
                                          fill=NA, size=1))

ggsave(filename = sprintf("%spredAge60to60ctrlVsHBraak.pdf",
                          outDir),
       width = 4,
       height = 5)

wilcox.test(df4DisContComp_60to70$pred_age[df4DisContComp_60to70$group == "control"],
            df4DisContComp_60to70$pred_age[df4DisContComp_60to70$group == "high_braak"])

# Predict ages in all the samples that are in rest and have braak index
dat_rest_wBraak <- dat_rest[rownames(dat_rest) %in% make.names(metDat$specimenID[!is.na(metDat$Braak)]), ]
wBraak_ages <- metDat$ageDeath[match(rownames(dat_rest_wBraak),
                                     make.names(metDat$specimenID))]
predicted_wBraak <- predictAge(model, t(dat_rest_wBraak))

predicted_wBraak[predicted_wBraak > 1] <- max(train4Mod$age_trans)
predicted_wBraak[predicted_wBraak < 0] <- .000012


predicted_wBraak2Age <- sapply(predicted_wBraak, function(x) agetrans_inverse(0, 110)(x))
predicted_wBraak2Age <- unlist(predicted_wBraak2Age)
names(predicted_wBraak2Age) <- gsub(".root", "", names(predicted_wBraak2Age))

dfBraak <- data.frame(specimenID = metDat$specimenID[match(names(predicted_wBraak2Age),
                                                           make.names(metDat$specimenID))],
                      age_death = wBraak_ages,
                      pred_age = predicted_wBraak2Age,
                      braak = metDat$Braak[match(names(predicted_wBraak2Age),
                                                 make.names(metDat$specimenID))],
                      group = "Rest",
                      substudy = metDat$substudy[match(names(predicted_wBraak2Age),
                                                       make.names(metDat$specimenID))])


df4DisContComp$braak <- rep(NA, nrow(df4DisContComp))
dfBraak <- rbind.data.frame(dfBraak,
                            df4DisContComp[df4DisContComp$group == "control",
                                           c("specimenID",
                                             "age_death",
                                             "pred_age",
                                             "braak",
                                             "group",
                                             "substudy")])

dfBraak$ageInt <- rep(NA, nrow(dfBraak))
ageInt <- 20:100
ageInt <- ageInt[ageInt %% 10 == 0]

for(i in 1:(length(ageInt) - 1)){
        age_1 <- ageInt[i]
        age_2 <- ageInt[i + 1]
        dfBraak$ageInt[dfBraak$age_death >= age_1 & dfBraak$age_death < age_2] <- sprintf("%s_%s", age_1, age_2)
        
}

dfBraak$braak[dfBraak$group == "control"] <- "control"


age60_70 <- ggplot(data = dfBraak[dfBraak$ageInt == "60_70", ],
                   mapping = aes(x = braak,
                                 y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(mapping = aes(col = substudy)) +
        ggtitle("60 to 70")

age70_80 <- ggplot(data = dfBraak[dfBraak$ageInt == "70_80", ],
                   mapping = aes(x = braak,
                                 y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(mapping = aes(col = substudy)) +
        ggtitle("70 to 80")

age80_90 <- ggplot(data = dfBraak[dfBraak$ageInt == "80_90", ],
                   mapping = aes(x = braak,
                                 y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(mapping = aes(col = substudy)) +
        ggtitle("80 to 90")

age90_100 <- ggplot(data = dfBraak[dfBraak$ageInt == "90_100", ],
                    mapping = aes(x = braak,
                                  y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(mapping = aes(col = substudy)) +
        ggtitle("90 to 100")

ggarrange(plotlist = list(age60_70, age70_80, age80_90, age90_100))
ggsave(filename = sprintf("%spredAg_byBraakAndAge_alph%s.pdf",
                          outDir,
                          alph),
       width = 10, height = 6)


age60_70 <- ggplot(data = dfBraak[dfBraak$ageInt == "60_70", ],
                   mapping = aes(x = braak,
                                 y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        ylab("predicted age") +
        xlab("Braak index") +
        ggtitle("60 to 70") +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, size=1))

age70_80 <- ggplot(data = dfBraak[dfBraak$ageInt == "70_80", ],
                   mapping = aes(x = braak,
                                 y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        ylab("predicted age") +
        xlab("Braak index") +
        ggtitle("70 to 80") +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, size=1))

age80_90 <- ggplot(data = dfBraak[dfBraak$ageInt == "80_90", ],
                   mapping = aes(x = braak,
                                 y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        ylab("predicted age") +
        xlab("Braak index") +
        ggtitle("80 to 90") +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, size=1))

age90_100 <- ggplot(data = dfBraak[dfBraak$ageInt == "90_100", ],
                    mapping = aes(x = braak,
                                  y = pred_age)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        ylab("predicted age") +
        xlab("Braak index") +
        ggtitle("90 to 100") +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=15),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray"), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black",
                                          fill=NA, size=1))

ggarrange(plotlist = list(age60_70, age70_80, age80_90, age90_100))
ggsave(filename = sprintf("%spredAg_byBraakAndAge_noBatch_alph%s.pdf",
                          outDir,
                          alph),
       width = 10, height = 6)
