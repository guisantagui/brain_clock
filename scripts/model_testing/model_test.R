################################################################################
# Brain clock: test predicted age of a pre-filtered set of samples, either     #
# perturbation or single cell.                                                 #
################################################################################
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
if(!require("h2o", quietly = T)){
        install.packages("h2o",
                         type="source",
                         repos="https://h2o-release.s3.amazonaws.com/h2o/rel-3.46.0/2/R")
}
if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}
if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
library(h2o)
library(plotUtils)
# Terminal argument parser
################################################################################
parser <- arg_parser("Train the GLM on transformed age and assess performance.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--modDir",
                               "--metDat",
                               "--sizeBatch",
                               "--whatSampsTest",
                               "--mem",
                               "--outDir"),
                       help = c("Input transcriptomic dataset, features columns, samples rows. CSV.",
                                "Directory where the h2o model is (generated with mod_train_and_test.R)",
                                "Metadata file (including ages).",
                                "Size of each batch for predictions.", 
                                "What samples should be tested. Possible values are perturbation or single_cell.",
                                "Memory to allocate to h2o instance.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

datFile <- parsed$input
modDir <- parsed$modDir
metDatFile <- parsed$metDat
sizeBatch <- as.numeric(parsed$sizeBatch)
whatSampsTest <- parsed$whatSampsTest
mem <- parsed$mem
outDir <- parsed$outDir


outName <- sprintf("%spred_ages.csv", outDir)

create_dir_if_not(outDir)

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

# Load the data
################################################################################
dat <- read_table_fast(datFile, row.names = 1)
metDat <- read_table_fast(metDatFile, row.names = 1)

# Keep only either perturbated samples or single cell samples
if(whatSampsTest == "perturbation"){
        filtSamps <- metDat$specimenID[metDat$diagn_4BrainClck == "perturbations"]
}else if(whatSampsTest == "single_cell"){
        filtSamps <- metDat$specimenID[metDat$substudy == "ageAnno"]
}

dat <- dat[make.names(filtSamps), ]


# Initialize h2o and load model
conn <- h2o.init(max_mem_size=mem)

# The script looks inside of the selected model directory, and loads
# the model file that is inside. It is coded like this because h2o gives a
# numerical number to the model, that each time gets larger. This way the
# launch script doesn't need to be changed in case we do a refit. If in the
# directory there are different models, it will load the latest (the one with
# the larger number).
modFiles <- list.files(modDir)

modFile <- sprintf("%s/%s",
                   modDir,
                   modFiles[length(modFiles)])

mod <- h2o.loadModel(modFile)

# Do predictions
################################################################################

# Load only desired samples based on indexes
nBatches <- ceiling(nrow(dat)/sizeBatch)

predVec <- c()
for(i in 1:nBatches){
        
        idx1 <- (i - 1) * sizeBatch + 1
        idx2 <- idx1 + sizeBatch - 1
        if(idx2 > nrow(dat)){
                idx2 <- nrow(dat)
        }
        pred <- predictAge(mod, t(dat[idx1:idx2, ]))
        predVec <- c(predVec, pred)
}


predsDF <- data.frame(specimenID = names(predVec),
                      chron_age = predVec)

writeCsvFst(predsDF, file = outName)
print(sprintf("%s saved in %s.", basename(outName), dirname(outName)))
h2o.shutdown(prompt = F)

################################################################################