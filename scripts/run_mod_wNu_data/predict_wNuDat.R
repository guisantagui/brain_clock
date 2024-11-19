library(h2o)
library(biomaRt)
library("GEOquery")
library(ggplot2)
if(!require(BiocManager, quietly = T)){
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
}
if(!require("limma", quietly = T)){
        BiocManager::install("limma", update = F)
}
library(limma)
library(DESeq2)

################################################################################
# Brain clock: predict new data.                                               #
################################################################################

# Functions 
################################################################################

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if (!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

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

# Predict the age given the expression matrix and the model
predictAge <- function(model,
                       predExpr){
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

get_gene_info <- function(id_list,
                          mart,
                          removeNAs = T,
                          filt = "entrezgene_id"){
        #filt = "entrezgene_id"
        #id_list <- dat$GeneID
        #removeNAs <- T
        #mart <- human
        if(filt == "ensembl_gene_id_version"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "ensembl_gene_id_version",
                          "entrezgene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }else{
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "entrezgene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }
        gene_coords <- getBM(attributes = atts,
                             filters = filt,
                             values = gsub("_PAR_Y", "", id_list),
                             mart = mart)
        
        gene_coords$size <- gene_coords$end_position - gene_coords$start_position
        
        # For cases where the same oid might map to several IDs, prioritize the
        # ones where the chromosome name is not some weird scaffold (i.e., the
        # number of characters of the chromosome name is less than 3).
        dupFilts <- gene_coords[, filt]
        dupFilts <- unique(dupFilts[duplicated(dupFilts)])
        for(i in seq_along(dupFilts)){
                f <- dupFilts[i]
                if(sum(nchar(gene_coords$chromosome_name) < 3 & gene_coords[, filt] == f) > 0){
                        gene_coords <- gene_coords[!(nchar(gene_coords$chromosome_name) > 2 & gene_coords[, filt] %in% dupFilts), ]
                }
        }
        
        mapped_matchVec <- match(gsub("_PAR_Y", "", id_list),
                                 gene_coords[, filt])
        outDF <- data.frame(matrix(nrow = length(id_list),
                                   ncol = ncol(gene_coords) + 1,
                                   dimnames = list(NULL,
                                                   c(sprintf("%s_orig", filt),
                                                     colnames(gene_coords)))))
        outDF[, sprintf("%s_orig", filt)] <- id_list
        for(n in colnames(gene_coords)){
                outDF[, n] <- gene_coords[mapped_matchVec, n]
        }
        unMapped <- outDF[, sprintf("%s_orig", filt)][is.na(outDF$ensembl_gene_id)]
        naVec <- is.na(outDF$ensembl_gene_id)
        if(removeNAs){
                nNAs <- sum(naVec)
                propNAs <- round(nNAs/nrow(outDF) * 100, digits = 2)
                outDF <- outDF[!naVec, ]
                print(sprintf("%s ensembl_gene_ids were removed due to inability to map to ENSEMBL database.", as.character(nNAs)))
                print(sprintf("This represents the %s%% of the total number of ensembl_gene_ids submitted (n = %s)",
                              as.character(propNAs),
                              as.character(length(id_list))))
        }
        return(outDF)
}

# This function adds up columns or rows that ended up with same name 
# after transformation to EMSEMBL.
sumDupsMat <- function(m, dups_in){
        if(dups_in == "cols"){
                m <- t(m)
        }
        featsOrd <- unique(rownames(m))
        dups <- unique(rownames(m)[duplicated(rownames(m))])
        
        for(i in seq_along(dups)){
                dup <- dups[i]
                dupMat <- m[rownames(m) == dup, ]
                sumsMat <- matrix(colSums(dupMat),
                                  nrow = 1,
                                  ncol = ncol(dupMat),
                                  dimnames = list(dup,
                                                  colnames(dupMat)))
                m <- m[rownames(m) != dup, ]
                m <- rbind(m, sumsMat)
        }
        m <- m[match(featsOrd, rownames(m)), ]
        if(dups_in == "cols"){
                m <- t(m)
        }
        return(m)
}

quantNorm <- function(m, axis = 2) {
        if (axis == 1) {
                m <- t(m)
        }
        
        # Get ranks of the matrix
        m_rank <- apply(m, 2, rank, ties.method = "average")
        
        # Sort the matrix
        m_sort <- apply(m, 2, sort)
        
        # Calculate row means of the sorted matrix
        means <- rowMeans(m_sort)
        
        # Create a normalized matrix with the same dimensions
        m_norm <- matrix(0, nrow = nrow(m), ncol = ncol(m))
        
        for (i in 1:ncol(m)) {
                m_norm[, i] <- means[rank(m[, i], ties.method = "average")]
        }
        
        if (axis == 1) {
                m_norm <- t(m_norm)
        }
        
        dimnames(m_norm) <- dimnames(m)
        return(m_norm)
}

parse_soft <- function(sft){
        gsms <- c()
        samps <- c()
        for(g in names(sft@gsms)){
                gsms <- c(gsms, g)
                samps <- c(samps, sft@gsms[[g]]@header$title)
        }
        out_df <- data.frame(GSM = gsms, sample = samps)
        return(out_df)
}

# Directory stuff
################################################################################
mod_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/modAllGenes_ingegWAllLincsBrain_and_sc_sva_chron_age_onSignGenes/modFuncsAlpha1/GLM_model_R_1727080896797_1"
sva_mod_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/svaMod/all_signChron/sva_pred_mod.rds"
mem <- "16G"
datFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/science_all.csv"
softFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/andrew_yoo/science_soft_all.csv"
bg_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_wLINCS_NPC_NEU_MIC_genes.csv"
dat_id_type <- "entrezgene_id"
do_frozen_sva <- F

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/yoo_data_science/"

createIfNot(outDir)


# Load data
################################################################################

# Initialize h2o and load model
conn <- h2o.init(max_mem_size=mem)
model <- h2o.loadModel(mod_file)
sva_mod <- readRDS(sva_mod_file)

dat <- readCsvFast(datFile)
#soft <- getGEO(filename = softFile)
#soft <- parse_soft(soft)
soft <- readCsvFast(softFile)
dfGenes <- read.csv(bg_file, row.names = 1)$ensembl_gene_id


#soft$pheno <- gsub(" ", "", gsub("[0-9]", "", soft$sample))
#soft$age <- as.numeric(gsub("[^0-9]", "", gsub("\\ .*", "", soft$sample)))

# Parse user data
################################################################################

# Change ID, if different than ensembl_gene_id
if(dat_id_type != "ensembl_gene_id"){
        gene_info_file <- sprintf("%s/gene_info.csv", dirname(datFile))
        if(!file.exists(gene_info_file)){
                human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
                gene_info <- get_gene_info(dat$GeneID, human, filt = dat_id_type)
                write.csv(gene_info, file = gene_info_file)
        }else{
                gene_info <- read.csv(gene_info_file, row.names = 1)
        }
        nu_IDs <- gene_info$ensembl_gene_id[match(dat$GeneID,
                                                  gene_info[, sprintf("%s_orig",
                                                                      dat_id_type)])]
        dat$GeneID <- nu_IDs
        dat <- dat[!is.na(dat$GeneID), ]
}

# Remove GeneID column, assign ids to rownames and, if there was a conversion
# to ENSEMBL, sum the counts of the genes that end up with the same ID.
rNams_dat <- dat$GeneID
dat <- as.matrix(dat[, colnames(dat) != "GeneID"])
rownames(dat) <- rNams_dat
dat <- sumDupsMat(dat, "rows")

# Keep only genes in LINCS, to make quantile normalization as close as possible
dat <- dat[rownames(dat) %in% dfGenes, ]

# Log2 transform and quantile normalize
dat <- log2(dat + 1)
dat <- quantNorm(dat)


mod_genes <- model@model$coefficients_table$names
mod_genes <- mod_genes[mod_genes != "Intercept"]

dat <- dat[rownames(dat) %in% mod_genes, ]

rownames(dat)[!rownames(dat) %in% rownames(sva_mod$coefficients)]
rownames(sva_mod$coefficients)[!rownames(sva_mod$coefficients) %in% rownames(dat)]

dat <- dat[match(mod_genes, rownames(dat)), ]
rownames(dat) <- mod_genes

which(is.na(dat), arr.ind = T)

# Predict SVs
if(do_frozen_sva){
        SVs <- predict(sva_mod, data.frame(t(dat)))
        
        # Regress out SVs
        dat <- removeBatchEffect(dat, covariates = SVs)
}

# Predict ages
pred_ages <- predictAge(model = model, dat)

soft$pred_age <- pred_ages[match(soft$GSM, names(pred_ages))]
soft <- soft[!is.na(soft$pred_age), ]

h2o.shutdown(prompt = F)

# Assess significance of relationship between predicted age and chronological
# age.

summary(lm(pred_age ~ age, data = soft[soft$treatment == "H2O", ]))
summary(lm(pred_age ~ age, data = soft[soft$treatment == "H2O" & soft$pheno == "HC", ]))

summary(lm(pred_age ~ age + treatment, data = soft))
summary(lm(pred_age ~ age + treatment + pheno, data = soft))

t.test(soft$pred_age[soft$treatment == "3TC" & soft$pheno == "LOAD"],
       soft$pred_age[soft$treatment == "H2O" & soft$pheno == "LOAD"])

boxplot(soft$pred_age[soft$treatment == "3TC" & soft$pheno == "LOAD"],
        soft$pred_age[soft$treatment == "H2O" & soft$pheno == "LOAD"])

#plot(soft$age[soft$GSE != "GSE253174" & soft$pheno == "HC"], soft$pred_age[soft$GSE != "GSE253174" & soft$pheno == "HC"])

ggplot(data = soft[soft$treatment == "H2O", ], mapping = aes(x = age, y = pred_age, col = pheno)) +
        geom_point()

ggsave(filename = sprintf("%sscatter_allPhen_noTreat.pdf", outDir),
       height = 5, width = 5)
ggsave(filename = sprintf("%sscatter_allPhen_noTreat.png", outDir),
       height = 5, width = 5)

ggplot(data = soft[soft$treatment == "H2O" & soft$pheno == "HC", ], mapping = aes(x = age, y = pred_age, col = pheno)) +
        geom_point()

ggsave(filename = sprintf("%sscatter_HC_noTreat.pdf", outDir),
       height = 5, width = 5)
ggsave(filename = sprintf("%sscatter_HC_noTreat.png", outDir),
       height = 5, width = 5)

ggplot(data = soft, mapping = aes(x = pheno, y = pred_age)) +
        geom_boxplot()

ggsave(filename = sprintf("%spheno_bxplt.pdf", outDir),
       height = 5, width = 5)
ggsave(filename = sprintf("%spheno_bxplt.png", outDir),
       height = 5, width = 5)

# Do a LM and an ANCOVA to assess if the influence of treatment on
# predicted age is significant.
mod_4Anc <- lm(pred_age ~ age * treatment, data = soft)

print("LM")
print(summary(mod_4Anc))

ancova_mod <- aov(pred_age ~ age + treatment + age:treatment + pheno + age:pheno,
                  data = soft)

print("ANCOVA")
print(summary(ancova_mod))

summary(aov(pred_age ~ age + treatment + pheno ,
            data = soft))

summary(aov(pred_age ~ age,
            data = soft[soft$treatment == "H2O" & soft$pheno == "HC", ]))



ggplot(soft, aes(x = age, y = pred_age, color = treatment)) +
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

ggsave(filename = sprintf("%sscatter_treat.pdf", outDir),
       height = 5, width = 5)
ggsave(filename = sprintf("%sscatter_treat.png", outDir),
       height = 5, width = 5)

ggplot(soft, aes(x = age, y = pred_age, color = pheno)) +
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
#ggsave(sprintf("%sctrls_highBraak_diagn_alpha%s_wLines.pdf", outDir, alph), 
#       height = 8,
#       width = 9)

ancova_mod_2 <- rstatix::anova_test(data = soft,
                                    formula = pred_age ~ age + treatment + age:treatment)

print(ancova_mod_2)

mod_coefs <- as.data.frame(model@model$coefficients_table)


mod_coefs <- mod_coefs[order(abs(mod_coefs$standardized_coefficients), decreasing = T), ]
mod_coefs <- mod_coefs[mod_coefs$standardized_coefficients != 0, ]

dat_4assessGenes <- data.frame(t(dat[rownames(dat) %in% mod_coefs$names, ]))

dat_4assessGenes$pred_age <- soft$pred_age[match(rownames(dat_4assessGenes),
                                                 soft$GSM)]

modList <- list()
pVals <- c()
for(gene in colnames(dat_4assessGenes)[1:(ncol(dat_4assessGenes) - 1)]){
        form <- as.formula(sprintf("pred_age ~ %s", gene))
        lmMod <- lm(form, data = dat_4assessGenes)
        modList[[gene]] <- lmMod
        pVals <- c(pVals, summary(lmMod)$coefficients[gene, "Pr(>|t|)"])
}
geneAssPVals <- data.frame(gene = colnames(dat_4assessGenes)[1:(ncol(dat_4assessGenes) - 1)],
                           p_val = pVals,
                           p_adj = p.adjust(pVals, method = "BH"))


geneAssPVals <- geneAssPVals[order(geneAssPVals$p_adj), ]
rownames(geneAssPVals) <- 1:nrow(geneAssPVals)
write.csv(geneAssPVals, file = sprintf("%sscience_genes_sign.csv", outDir))
write.csv(mod_coefs, file = sprintf("%smod_coefs.csv", outDir))
