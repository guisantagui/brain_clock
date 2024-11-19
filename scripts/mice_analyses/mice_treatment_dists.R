if(!require(TabulaMurisSenisData, quietly = T)){
        devtools::install_github("fmicompbio/TabulaMurisSenisData")
}
library(TabulaMurisSenisData)
if(!require(org.Mm.eg.db, quietly = T)){
        BiocManager::install("org.Mm.eg.db", update = F)
}
if(!require(biomaRt, quietly = T)){
        BiocManager::install("biomaRt", update = F)
}

library(DESeq2)

library(org.Mm.eg.db)
library(biomaRt)
library(ggplot2)
library(sva)
library(SmartSVA)
library(ropls)
library(ggpubr)

stand <- function(m, axis = 1){
        if (axis == 1){
                m <- t(m)
        }
        m_stand <- m[, apply(m,
                             2,
                             function(x) sd(x) != 0)]
        m_stand <- apply(m, 2, function(x) (x - mean(x))/sd(x))
        if (axis ==1){
                m_stand <- t(m_stand)
        }
        return(m_stand)
}

getTopContrib <- function(PC, topN = 12, x = "PC1", y = "PC2"){
        contrib <- facto_summarize(PC, 
                                   "var", 
                                   axes = c(as.numeric(gsub("PC", "", x)),
                                            as.numeric(gsub("PC", "", y))))
        contrib <- contrib[order(contrib$contrib, decreasing = T), 
                           c("name", "contrib")]
        topContrib <- as.character(contrib$name[1:topN])
        return(topContrib)
}

plot_pca <- function(PC, sample_info, col = NULL, shape = NULL,
                     x = "PC1", y = "PC2",
                     labs = NULL, topNFeats = NULL, biplot = F,
                     coord_fix = T){
        dat <- data.frame(obsnames=row.names(PC$x), PC$x)
        dat <- dat[, c("obsnames", x, y)]
        
        dat <- cbind.data.frame(dat,
                                sample_info[match(make.names(dat$obsnames),
                                                  make.names(sample_info$sample)), ])
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        
        datapc <- data.frame(varnames=rownames(PC$rotation), 
                             PC$rotation)
        mult <- min(
                (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
                (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
        datapc <- transform(datapc,
                            v1 = .7 * mult * (get(x)),
                            v2 = .7 * mult * (get(y))
        )
        datapc$x0 <- rep(0, nrow(datapc))
        datapc$y0 <- rep(0, nrow(datapc))
        if(!is.null(topNFeats)){
                varPlotFilt <- getTopContrib(PC, topN = topNFeats, x = x, y = y)
                datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
        }
        
        if(!is.null(col)){
                col <- sym(col)
        }
        
        if(!is.null(shape)){
                shape <- sym(shape)
        }
        
        if(!is.null(labs)){
                labs <- sym(labs)
        }
        pcaPlt <- ggplot(dat, aes(!!sym(x), !!sym(y),
                                  color = !!col,
                                  label = !!labs,
                                  shape = !!shape)) +
                geom_point(size=3) +
                xlab(sprintf("%s (%s %%)", x, propX)) +
                ylab(sprintf("%s (%s %%)", y, propY)) + 
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if(coord_fix){
                pcaPlt <- pcaPlt +
                        coord_fixed()
        }
        if (!is.null(labs)){
                pcaPlt <- pcaPlt +
                        geom_text_repel()
        }
        if (biplot){
                pcaPlt <- pcaPlt +
                        geom_text_repel(data=datapc, 
                                        aes(x=v1, y=v2, label=varnames), 
                                        color = "black", 
                                        size = 3,
                                        max.overlaps = 100,
                                        inherit.aes = F) + 
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
                                     size = 0.5,
                                     inherit.aes = F)
        }
        return(pcaPlt)
}

doPCAMultiPlot <- function(PC,
                           sample_info,
                           nComps,
                           col = "age",
                           shape = "tissue",
                           labs = NULL,
                           topNFeats = NULL,
                           biplot = F,
                           coord_fix = T){
        plotList <- list()
        for(j in 2:(nComps)){
                for(i in 1:(nComps - 1)){
                        if(j > i){
                                scPlot <- plot_pca(PC,
                                                   biplot = biplot,
                                                   x = sprintf("PC%s", i),
                                                   y = sprintf("PC%s", j),
                                                   sample_info = sample_info,
                                                   col = col,
                                                   shape = shape,
                                                   labs = labs,
                                                   coord_fix = coord_fix)
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

expFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/mouse_experiments/RNAseq_Hippocampus_Cortex.RData"
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/mice_experiments/"

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

load(expFile)

tms_bulk <- TabulaMurisSenisBulk()

inp_df[1:10, 1:10]

table(tms_bulk$organ)

counts <- tms_bulk@assays@data$counts

counts[1:10, 1:10]

metDat <- tms_bulk@colData
metDat <- metDat[metDat$organ == "Brain", ]

counts <- counts[, colnames(counts) %in% metDat$`Sample name`]

# Convert names to ENSEMBL
################################################################################

# Create mart for conversion of IDs
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id",
                                "external_gene_name",
                                "hgnc_symbol",
                                "entrezgene_id",
                                "refseq_mrna",
                                "description",
                                "chromosome_name",
                                "start_position",
                                "end_position"),
                 filters = "external_gene_name",
                 values = rownames(counts),     
                 mart = mouse)

unmaped <- rownames(counts)[!rownames(counts) %in% results$external_gene_name]

length(unmaped)

res_unmapped <- getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name",
                                     "hgnc_symbol",
                                     "entrezgene_id",
                                     "refseq_mrna",
                                     "description",
                                     "chromosome_name",
                                     "start_position",
                                     "end_position"),
                 filters = "ensembl_transcript_id",
                 values = gsub("\\..*", "", unmaped),
                 mart = mouse)

ensembl_names <- results$ensembl_gene_id[match(rownames(counts),
                                               results$external_gene_name)]
counts <- counts[!is.na(ensembl_names), ]
ensembl_names <- ensembl_names[!is.na(ensembl_names)]
rownames(counts) <- ensembl_names

# Edit rownames to remove version and keep only cortex samples
rownames(inp_df) <- gsub("\\..*", "", rownames(inp_df))
inp_df <- inp_df[, meta_df$Sample[grep("Cortex", meta_df$Condition)]]

# Merge datasets
comm_genes <- intersect(rownames(inp_df), rownames(counts))
merged <- cbind.data.frame(inp_df[comm_genes, ], counts[comm_genes, ])

# Remove not expressed genes
merged <- merged[rowSums(merged) != 0, ]

merged_log2 <- log2(merged + 1)


merged_log2_stand <- stand(merged_log2, axis = 1)

merged_log2_stand_pca <- prcomp(t(merged_log2_stand),
                                scale. = F,
                                center = F)


metDatComb <- meta_df
metDatComb <- metDatComb[grep("Cortex", metDatComb$Condition), ]
metDatComb$treatment <- gsub(".*_", "", metDatComb$Condition)
metDatComb$tissue <- gsub("\\_.*", "", metDatComb$Condition)
metDatComb <- metDatComb[, !colnames(metDatComb) %in% c("Condition",
                                                       "Group")]
metDatComb$age <- rep(18, nrow(metDatComb))

toBind <- data.frame(Sample = metDat$`Sample name`,
                     treatment = "none",
                     tissue = metDat$organ,
                     age = as.numeric(metDat$`characteristics: age`))
metDatComb <- rbind.data.frame(metDatComb, toBind)

colnames(metDatComb) <- tolower(colnames(metDatComb))

plot_pca(merged_log2_stand_pca, metDatComb, col = "age", shape = "tissue",
         x = "PC1", y = "PC2",
         labs = NULL, topNFeats = NULL, biplot = F)

ggsave(filename = sprintf("%spca_merged_log2_stand.pdf", outDir),
       height = 5,
       width = 5)

merged_log2_stand_tms_pca <- prcomp(t(merged_log2[, metDatComb$sample[metDatComb$treatment == "none"]]),
                                    scale. = F,
                                    center = T)

plot_pca(merged_log2_stand_tms_pca, metDatComb, col = "age", shape = "tissue",
         x = "PC1", y = "PC2",
         labs = NULL, topNFeats = NULL, biplot = F)

ggsave(filename = sprintf("%spca_tms_log2_stand.pdf", outDir),
       height = 5,
       width = 5)

doPCAMultiPlot(merged_log2_stand_tms_pca,
               metDatComb,
               6,
               col = "age",
               shape = "tissue",
               labs = NULL,
               topNFeats = NULL,
               biplot = F,
               coord_fix = F)

ggsave(filename = sprintf("%spcaMult_tms_log2_stand.pdf", outDir),
       height = 10,
       width = 10)


# Run SVA to remove batch effect
pheno <- data.frame(sample = colnames(merged_log2),
                    age = metDatComb$age[match(colnames(merged_log2),
                                               metDatComb$sample)])

mod <- model.matrix(~age, data = pheno)
mod0 <- model.matrix(~1, data = pheno)

rownames(mod) <- pheno$sample
rownames(mod0) <- pheno$sample

n.sv <- num.sv(merged_log2, mod, method = "leek")

svobj <- smartsva.cpp(as.matrix(merged_log2),
                      mod,
                      mod0,
                      n.sv = n.sv)


merged_log2_svAdj <- removeBatchEffect(merged_log2, covariates = svobj$sv)

merged_log2_svAdj_stand <- stand(merged_log2_svAdj, axis = 1)
merged_log2_svAdj_stand_pca <- prcomp(t(merged_log2_svAdj),
                                      scale. = F,
                                      center = T)

plot_pca(merged_log2_svAdj_stand_pca,
         metDatComb, col = "age", shape = "tissue",
         x = "PC1", y = "PC2",
         labs = NULL, topNFeats = NULL, biplot = F)

ggsave(filename = sprintf("%spca_merged_log2_svaAdj_stand.pdf", outDir),
       height = 5,
       width = 5)

doPCAMultiPlot(merged_log2_svAdj_stand_pca,
               metDatComb,
               6,
               col = "age",
               shape = "tissue",
               labs = NULL,
               topNFeats = NULL,
               biplot = F,
               coord_fix = F)

ggsave(filename = sprintf("%spcaMult_merged_log2_svaAdj_stand.pdf", outDir),
       height = 10,
       width = 10)

plot_pca(merged_log2_svAdj_stand_pca,
         metDatComb, col = "age", shape = "tissue",
         x = "PC5", y = "PC6",
         labs = NULL, topNFeats = NULL, biplot = F)


# Keep only the ones that are 18 months old, and the young ones

metDatSubset <- metDatComb[metDatComb$age >= 18 | metDatComb$age == 3 | metDatComb$age == 6, ]

metDatSubset <- metDatSubset[!(metDatSubset$treatment == "none" & metDatSubset$age == 18), ]

metDatSubset$group <- metDatSubset$treatment

metDatSubset$group[metDatSubset$treatment == "none" & metDatSubset$age > 18] <- "old"
metDatSubset$group[metDatSubset$treatment == "none" & metDatSubset$age <= 6] <- "young"


merged_log2_svAdj_subset <- merged_log2_svAdj[, metDatSubset$sample]
merged_log2_svAdj_subset_stand <- stand(merged_log2_svAdj_subset, axis = 1)
merged_log2_svAdj_subset_stand_pca <- prcomp(t(merged_log2_svAdj_subset_stand),
                                             scale. = F,
                                             center = F)

plot_pca(merged_log2_svAdj_subset_stand_pca,
         metDatSubset, col = "group", shape = "tissue",
         x = "PC1", y = "PC2",
         labs = NULL, topNFeats = NULL, biplot = F)

ggsave(filename = sprintf("%spca_subset_log2_svaAdj_stand.pdf", outDir),
       height = 5,
       width = 5)

doPCAMultiPlot(merged_log2_svAdj_subset_stand_pca,
               metDatSubset,
               6,
               col = "group",
               shape = "tissue",
               labs = NULL,
               topNFeats = NULL,
               biplot = F,
               coord_fix = F)

ggsave(filename = sprintf("%spcaMult_subset_log2_svaAdj_stand.pdf", outDir),
       height = 10,
       width = 10)

control <- metDatSubset$sample[metDatSubset$group == "Control"]
treated <- metDatSubset$sample[metDatSubset$group == "Treated"]
young <- metDatSubset$sample[metDatSubset$group == "young"]
old <- metDatSubset$sample[metDatSubset$group == "old"]

subset_dist <- as.matrix(dist(merged_log2_svAdj_subset_stand_pca$x[, 1]))

ctrl_vs_young <- c()
ctrl_vs_old <- c()
for(ctrl in control){
        for(y in young){
                d <- subset_dist[ctrl, y]
                names(d) <- sprintf("%s_vs_%s", ctrl, y)
                ctrl_vs_young <- c(ctrl_vs_young, d)
        }
        for(o in old){
                d <- subset_dist[ctrl, o]
                names(d) <- sprintf("%s_vs_%s", ctrl, o)
                ctrl_vs_old <- c(ctrl_vs_old, d)
        }
}

treat_vs_young <- c()
treat_vs_old <- c()
for(treat in treated){
        for(y in young){
                d <- subset_dist[treat, y]
                names(d) <- sprintf("%s_vs_%s", treat, y)
                treat_vs_young <- c(treat_vs_young, d)
        }
        for(o in old){
                d <- subset_dist[treat, o]
                names(d) <- sprintf("%s_vs_%s", treat, o)
                treat_vs_old <- c(treat_vs_old, d)
        }
}

old_bxplt_df <- data.frame(comp = c(names(ctrl_vs_old),
                                    names(treat_vs_old)),
                           dist = c(ctrl_vs_old,
                                    treat_vs_old),
                           group = c(rep("control", length(ctrl_vs_old)),
                                     rep("treated", length(treat_vs_old))))

ggplot(old_bxplt_df, aes(x = group, y = dist)) +
        geom_boxplot(outliers = F) +
        geom_jitter() +
        ylab("distance to old") +
        theme(title = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA,
                                          linewidth = 1),
              panel.grid.major = element_line(colour = "#d4d4d4"),
              legend.position = "right")

ggsave(filename = sprintf("%sdist_bxplt_old.pdf", outDir), height = 5,
       width = 5)

young_bxplt_df <- data.frame(comp = c(names(ctrl_vs_young),
                                      names(treat_vs_young)),
                             dist = c(ctrl_vs_young,
                                      treat_vs_young),
                             group = c(rep("control", length(ctrl_vs_young)),
                                       rep("treated", length(treat_vs_young))))

ggplot(young_bxplt_df, aes(x = group, y = dist)) +
        geom_boxplot(outliers = F) +
        geom_jitter() +
        ylab("Distance to young") +
        theme(title = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA,
                                          linewidth = 1),
              panel.grid.major = element_line(colour = "#d4d4d4"),
              legend.position = "right")

ggsave(filename = sprintf("%sdist_bxplt_young.pdf", outDir), height = 5,
       width = 5)

t.test(ctrl_vs_old, treat_vs_old)
t.test(ctrl_vs_young, treat_vs_young)

wilcox.test(ctrl_vs_old, treat_vs_old)
wilcox.test(ctrl_vs_young, treat_vs_young)

#ropls::opls(t(merged_log2_svAdj[, metDatComb$sample]),
#            y = metDatComb$age,
#            predI = 1,
#            orthoI = 3)

#ropls::opls(t(merged_log2_svAdj_subset[, metDatSubset$sample[metDatSubset$group == "old" | metDatSubset$group == "young"]]),
#            y = metDatSubset$group[metDatSubset$group == "old" | metDatSubset$group == "young"],
#            predI = 1,
#            orthoI = 3)
merged_log2_svAdj_subset

# Do a DE analysis with DESeq between young and old to see if the distances
# in these genes decrease

metdat_4DESeq <- metDatSubset[metDatSubset$group == "old" | metDatSubset$group == "young", ]
rownames(metdat_4DESeq) <- metdat_4DESeq$sample
metdat_4DESeq$group <- factor(metdat_4DESeq$group, levels = c("young", "old"))

counts_old_vs_young <- counts[, metdat_4DESeq$sample]
counts_old_vs_young <- counts_old_vs_young[rowSums(counts_old_vs_young) != 0, ]
dds <- DESeqDataSetFromMatrix(countData = counts_old_vs_young,
                              colData = metdat_4DESeq,
                              design = ~group)

DE <- DESeq(dds,
            test = "Wald",
            fitType = "glmGamPoi",
            sfType = "ratio")
res <- data.frame(results(DE))
res <- res[!is.na(res$padj), ]
res_sign <- res[res$padj < 0.05, ]


metdat_4DESeq_treatment <- metDatSubset[metDatSubset$group == "Treated" | metDatSubset$group == "Control", ]
rownames(metdat_4DESeq_treatment) <- metdat_4DESeq_treatment$sample
metdat_4DESeq_treatment$group <- factor(metdat_4DESeq_treatment$group,
                                        levels = c("Control", "Treated"))

counts_treat <- merged[, metdat_4DESeq_treatment$sample]
counts_treat <- counts_treat[rowSums(counts_treat) != 0, ]
dds_treat <- DESeqDataSetFromMatrix(countData = counts_treat,
                                    colData = metdat_4DESeq_treatment,
                                    design = ~group)

DE_treat <- DESeq(dds_treat,
                  test = "Wald",
                  fitType = "glmGamPoi",
                  sfType = "ratio")
res_treat <- data.frame(results(DE_treat))
res_treat_sign <- res_treat[!is.na(res_treat$padj), ]
res_treat_sign <- res_treat_sign[res_treat_sign$padj < 0.05, ]

res_treat_log2FC <- res_treat$log2FoldChange
names(res_treat_log2FC) <- rownames(res_treat)
res_treat_log2FC[is.na(res_treat_log2FC)] <- 0

res_sign_log2FC <- res_sign$log2FoldChange
names(res_sign_log2FC) <- rownames(res_sign)

res_treat_log2FC <- res_treat_log2FC[match(names(res_sign_log2FC),
                                           names(res_treat_log2FC))]

is_na <- is.na(names(res_treat_log2FC))

res_sign_log2FC <- res_sign_log2FC[!is_na]
res_treat_log2FC <- res_treat_log2FC[!is_na]
# Check if what DE genes of old vs young follow the same direction than 
# rejuvenation after treatment. We need to multiply the old vs young log2FC
# because the comparison was towards aging, so positive log2FC are more
# expressed in old samples.

num_its <- 10000
propVec <- c()
for (n in 1:num_its){
        dir_vec <- c()
        res_treat_randLog2FC <- sample(res_treat$log2FoldChange,
                                       length(res_sign_log2FC))
        res_treat_randLog2FC[is.na(res_treat_randLog2FC)] <- 0
        names(res_treat_randLog2FC) <- names(res_sign_log2FC)
        for (i in seq_along(res_treat_log2FC)){
                gene <- names(res_sign_log2FC)[i]
                if ((res_treat_randLog2FC[gene] > 0 & - 1 * res_sign_log2FC[gene] > 0) ||
                    (res_treat_randLog2FC[gene] < 0 & - 1 * res_sign_log2FC[gene] < 0)){
                        dir_vec <- c(dir_vec, 1)
                }else{
                        dir_vec <- c(dir_vec, 0)
                }
        }
        prop <- sum(dir_vec)/length(dir_vec)
        propVec <- c(propVec, prop)
}

hist(propVec)

table(metDatComb$age)
dir_vec <- c()
for (i in seq_along(res_treat_log2FC)){
        #i <- 1
        gene <- names(res_treat_log2FC)[i]
        if ((res_treat_log2FC[gene] > 0 & -1 * res_sign_log2FC[gene] > 0) || (res_treat_log2FC[gene] < 0 & -1 * res_sign_log2FC[gene] < 0)){
                dir_vec <- c(dir_vec, 1)
        }else{
                dir_vec <- c(dir_vec, 0)
        }
}
names(dir_vec) <- names(res_treat_log2FC)

sum(propVec > sum(dir_vec)/length(dir_vec))/num_its

binom.test(x = sum(dir_vec),
           n = length(dir_vec),
           alternative = "greater",
           p = 0.3)

plotCounts(dds, gene = "ENSMUSG00000073590", intgroup = "group")


merged_log2_svAdj_signGenes <- merged_log2_svAdj[rownames(res_sign), ]

merged_log2_svAdj_signGenes_stand <- stand(merged_log2_svAdj_signGenes,
                                           axis = 1)

plot_pca(prcomp(t(merged_log2_svAdj_signGenes),
                scale. = T,
                center = T),
         metDatComb, col = "age", shape = "treatment",
         x = "PC1", y = "PC2",
         labs = NULL, topNFeats = NULL, biplot = F)
