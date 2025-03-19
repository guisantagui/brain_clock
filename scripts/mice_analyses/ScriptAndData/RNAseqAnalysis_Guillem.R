library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(ggrepel)

if(!require("WebGestaltR", quietly = T)) install.packages("WebGestaltR")
library(WebGestaltR)



setwd("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/scripts/ScriptAndData")

# files <- list.files(pattern = "featurecounts.txt") 
# 
# inp_list <- list()
# for(f in files){
#   n <- strsplit(f,"_")[[1]][1]
#   df <- read.table(f,header = T, stringsAsFactors = F)
#   colnames(df)[ncol(df)] <- n
#   rownames(df) <- df$Geneid
#   df <- df[,ncol(df),drop = F]
#   inp_list[[n]] <- df
# }
# 
# 
# inp_df <- do.call("cbind.data.frame",inp_list)
# 
# meta_df <- data.frame(Sample = colnames(inp_df),
#                       Group = c("2","4","4","2","4","4","4","1","3","3","3","3","1","3"),
#                       Condition = c("Cortex_Treated","Hippo_Treated","Hippo_Treated","Cortex_Treated",
#                                     "Hippo_Treated","Hippo_Treated","Hippo_Treated","Cortex_Control",
#                                     "Hippo_Control","Hippo_Control","Hippo_Control","Hippo_Control",
#                                     "Cortex_Control","Hippo_Control"),
#                       stringsAsFactors = F)
# rownames(meta_df) <- meta_df$Sample
# save(inp_df,meta_df,file = "Data.RData")

load("Data.RData")
table(meta_df$Condition)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
ensToSymb <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=gsub("\\.[0-9]+","",rownames(inp_df)),mart= mart)

### Hippocampus
meta_hippo_df <- meta_df[meta_df$Condition %in% c("Hippo_Control","Hippo_Treated"),]
inp_hippo_df <- inp_df[,meta_hippo_df$Sample]

dds_hippo <- DESeqDataSetFromMatrix(countData = as.matrix(inp_hippo_df),
                                    colData = meta_hippo_df,
                                    design = ~ Condition)

smallestGroupSize <- 5
keep <- rowSums(counts(dds_hippo) >= 5) >= smallestGroupSize
dds_hippo <- dds_hippo[keep,]

dds_hippo <- DESeq(dds_hippo,
                   test = "Wald",
                   fitType = "glmGamPoi",
                   sfType = "ratio")

diffExp_hippo_df <- as.data.frame(results(dds_hippo,contrast = c("Condition","Hippo_Treated","Hippo_Control")))
diffExp_hippo_df$Gene <- gsub("\\.[0-9]+","",rownames(diffExp_hippo_df))
diffExp_hippo_df <- merge(diffExp_hippo_df,ensToSymb,by.x = "Gene",by.y = "ensembl_gene_id", all.x = T, all.y = F)

fc_hippo_vec <- diffExp_hippo_df$log2FoldChange
names(fc_hippo_vec) <- diffExp_hippo_df$Gene

#CORTEX
meta_cortex_df <- meta_df[meta_df$Condition %in% c("Cortex_Control","Cortex_Treated"),]
inp_cortex_df <- inp_df[,meta_cortex_df$Sample]

dds_cortex <- DESeqDataSetFromMatrix(countData = as.matrix(inp_cortex_df),
                                    colData = meta_cortex_df,
                                    design = ~ Condition)

smallestGroupSize <- 2
keep <- rowSums(counts(dds_cortex) >= 5) >= smallestGroupSize
dds_cortex <- dds_cortex[keep,]

dds_cortex <- DESeq(dds_cortex,
                    test = "Wald",
                    fitType = "glmGamPoi",
                    sfType = "ratio")

diffExp_cortex_df <- as.data.frame(results(dds_cortex,contrast = c("Condition","Cortex_Treated","Cortex_Control")))
diffExp_cortex_df$Gene <- gsub("\\.[0-9]+","",rownames(diffExp_cortex_df))
diffExp_cortex_df <- merge(diffExp_cortex_df,ensToSymb,by.x = "Gene",by.y = "ensembl_gene_id", all.x = T, all.y = F)

fc_cortex_vec <- diffExp_cortex_df$log2FoldChange
names(fc_cortex_vec) <- diffExp_cortex_df$Gene

log10p_cortex_vec <- diffExp_cortex_df$pvalue
names(log10p_cortex_vec) <- diffExp_cortex_df$Gene
log10p_cortex_vec[is.na(log10p_cortex_vec)] <- 1
log10p_cortex_vec <- -log10(log10p_cortex_vec)

fc_cortex_df <- data.frame(V1 = names(fc_cortex_vec),
                           V2 = unname(fc_cortex_vec*log10p_cortex_vec),
                           stringsAsFactors = F)

res <- results(dds_cortex,contrast = c("Condition","Cortex_Treated","Cortex_Control"))
p <- EnhancedVolcano(diffExp_cortex_df,
                lab = diffExp_cortex_df$mgi_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c("Myoc","Thbs4","Xlr4b","Gabrb2","Aqp4","Mobp","Lrrk2","Tnnt2","Pla2g4e"),#diffExp_cortex_df$mgi_symbol[which(diffExp_cortex_df$padj < 1e-5)],
                xlim = c(min(diffExp_cortex_df$log2FoldChange[which(!is.na(diffExp_cortex_df$padj))])-1,max(diffExp_cortex_df$log2FoldChange[which(!is.na(diffExp_cortex_df$padj))])+1),
                ylim = c(0, max(-log10(diffExp_cortex_df$padj), na.rm = TRUE) + 2),
                ylab = bquote(~-Log[10] ~ italic('adj. P')),
                pCutoff = 0.05,
                FCcutoff = 0,
                title = '',
                subtitle = '',
                boxedLabels = T,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colConnectors = 'black',
                arrowheads = F,
                directionConnectors = 'both',
                #min.segment.length = 5,
                colGradient = c('red2','blue2'),
                max.overlaps = Inf) +
  theme(legend.position = 'none')
p

ggsave("Volcano_Cortex.png",p)

gsea_aging <- WebGestaltR(enrichMethod = "GSEA",
                          organism = "mmusculus",
                          enrichDatabaseFile = "m8.all.v2024.1.Mm.symbols.gmt",
                          enrichDatabaseType = "genesymbol",
                          interestGene = fc_cortex_df,#[which(diffExp_cortex_df$padj < 0.05),],
                          interestGeneType = "ensembl_gene_id",
                          referenceGene = gsub("\\.[0-9]+","",rownames(inp_df)),
                          referenceGeneType = "ensembl_gene_id",
                          sigMethod = "fdr",
                          fdrThr = 0.05,
                          maxNum = 10000,
                          minNum = 3,
                          isOutput = F)

saveRDS(gsea_aging,file = "GSEA_MSigDB_M8_Cortex.rds")

readLines("m8.all.v2024.1.Mm.symbols.gmt")


clockgenes <- read.table("ClockGeneswCoeffs.txt",header = T, sep = "\t", stringsAsFactors = F)
human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl',host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl',host = "https://dec2021.archive.ensembl.org/")
annot_table <- getLDS(
  mart = human,
  attributes = c('ensembl_gene_id'),
  martL = mouse,
  attributesL = c('ensembl_gene_id'),
  filters = 'ensembl_gene_id',
  values = clockgenes$ensembl_gene_id) #2 clock genes cannot be found

clockgenes <- merge(clockgenes,annot_table,by.x = 1, by.y = 1)


inp_tmp_cortex_df <- inp_cortex_df
inp_tmp_cortex_df <- t(t(inp_tmp_cortex_df*1e6)/colSums(inp_tmp_cortex_df))
cortex_treated <- rowMeans(inp_tmp_cortex_df[,meta_df$Sample[which(meta_df$Condition == "Cortex_Treated")]])
names(cortex_treated) <- gsub("\\.[0-9]+","",names(cortex_treated))
cortex_untreated <- rowMeans(inp_tmp_cortex_df[,meta_df$Sample[which(meta_df$Condition == "Cortex_Control")]])
names(cortex_untreated) <- gsub("\\.[0-9]+","",names(cortex_untreated))

cortex_log2fc <- log2((cortex_treated+1)/(cortex_untreated+1))

length(intersect(clockgenes$Gene.stable.ID.1,names(cortex_log2fc))) #421 genes are in the mouse data

clockcoeffs <- clockgenes$coefficients
names(clockcoeffs) <- clockgenes$Gene.stable.ID.1

commonGenes <- intersect(names(clockcoeffs),names(cortex_log2fc))

clockcoeffs <- clockcoeffs[commonGenes]
clockcoeffs <- gsub(",","\\.",clockcoeffs)
clockcoeffs <- as.numeric(clockcoeffs)
cortex_log2fc <- cortex_log2fc[commonGenes]

all(names(clockcoeffs) == names(cortex_log2fc))

library(ggrepel)
expChanges_inp <- data.frame("Untreated" = unname(cortex_untreated[commonGenes]),
                             "Treated" = unname(cortex_treated[commonGenes]),
                             "Genes" = commonGenes,
                             stringsAsFactors = F)
expChanges_inp$Diff <- expChanges_inp$Treated-expChanges_inp$Untreated
expChanges_inp <- merge(expChanges_inp,ensToSymb,by.x = "Genes", by.y = 1,all.x = T)
expChanges_inp$Clockcoeffs <- clockcoeffs
expChanges_inp$DiffWeighted <- expChanges_inp$Diff * expChanges_inp$Clockcoeffs
expChanges_inp$Label <- ifelse((nrow(expChanges_inp)-rank(abs(expChanges_inp$DiffWeighted))) <= 10,expChanges_inp$mgi_symbol,"")

p <- ggplot(expChanges_inp, aes(x=Untreated, y=Treated, label = Label)) + 
  geom_point()+
  geom_smooth(method = lm) +
  geom_text_repel(max.overlaps = Inf,box.padding = 0.5) +
  theme_minimal()
p

ggsave("ClockGeneExpression_Cortex_Label10.png",p)


cor(clockcoeffs, cortex_log2fc)

predAge_cortex <- c("Treated" = 0, "Untreated" = 0)
predAge_cortex["Treated"] <- sum(clockcoeffs * cortex_treated[commonGenes]) # 79.06771
predAge_cortex["Untreated"] <- sum(clockcoeffs * cortex_untreated[commonGenes]) # 113.4904

saveRDS(predAge_cortex,file = "PredAge_Cortex.rds")

table(sign(clockcoeffs),sign(cortex_treated[commonGenes]-cortex_untreated[commonGenes]))
#      -1   0   1
#  -1  97  18 127
#   1  77  13  89
# 31 genes do not change at all
# 186 genes change towards older phenotype
# 204 genes change towards younger phenotype

# Do predictions at the level of samples
inp_tmp_cortex_df_4Preds <- inp_tmp_cortex_df
rownames(inp_tmp_cortex_df_4Preds) <- gsub("\\.[0-9]+",
                                           "",
                                           rownames(inp_tmp_cortex_df_4Preds))

inp_tmp_cortex_df_4Preds <- inp_tmp_cortex_df_4Preds[commonGenes, ]

preds_df <- data.frame(sample = colnames(inp_tmp_cortex_df_4Preds),
                       treatment = meta_cortex_df$Condition[match(colnames(inp_tmp_cortex_df_4Preds), meta_cortex_df$Sample)],
                       pred_age = apply(inp_tmp_cortex_df_4Preds,
                                        2,
                                        function(x) sum(clockcoeffs * x)))
preds_df$treatment <- gsub(".*_", "", preds_df$treatment)

t.test(preds_df$pred_age[preds_df$treatment == "Treated"],
       preds_df$pred_age[preds_df$treatment == "Control"])

preds_bxplt <- ggplot(preds_df, mapping = aes(x = treatment, y = pred_age)) +
        geom_boxplot() +
        ylab("Predicted age") +
        theme(title = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              axis.title.x = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA,
                                          linewidth = 1),
              panel.grid.major = element_line(colour = "#d4d4d4"),
              legend.position = "right")
preds_bxplt
ggsave(plot = preds_bxplt, filename = "preds_bxplt.pdf", height = 5, width = 3)

idx_upandcorrect <- names(which(sign(clockcoeffs) == -1 & sign(cortex_treated[commonGenes]-cortex_untreated[commonGenes]) == 1))
idx_downandcorrect <- names(which(sign(clockcoeffs) == 1 & sign(cortex_treated[commonGenes]-cortex_untreated[commonGenes]) == -1))

inp <- expChanges_inp$DiffWeighted
names(inp) <- expChanges_inp$Genes
inp <- data.frame(V1 = expChanges_inp$Genes,
                  V2 = expChanges_inp$DiffWeighted)

out <- WebGestaltR(enrichMethod = "GSEA",
                  organism = "mmusculus",
                  enrichDatabase = "geneontology_Biological_Process",
                  interestGene = inp,#c(idx_downandcorrect,idx_upandcorrect),
                  interestGeneType = "ensembl_gene_id",
                  referenceGene = commonGenes,#gsub("\\.[0-9]+","",rownames(inp_df)),
                  referenceGeneType = "ensembl_gene_id",
                  sigMethod = "fdr",
                  fdrThr = 0.05,
                  maxNum = 10000,
                  minNum = 3,
                  isOutput = F)

saveRDS(out, file = "GSEA_GOBP_Clockgenes.rds")


ora_gobp <- WebGestaltR(enrichMethod = "ORA",
organism = "mmusculus",
enrichDatabase = "geneontology_Biological_Process",
interestGene = diffExp_cortex_df$Gene[which(diffExp_cortex_df$padj < 0.05)],
interestGeneType = "ensembl_gene_id",
referenceGene = gsub("\\.[0-9]+","",rownames(inp_df)),
referenceGeneType = "ensembl_gene_id",
sigMethod = "fdr",
fdrThr = 0.05,
minNum = 3,
maxNum = 10000,
isOutput = F)

saveRDS(ora_gobp, file = "ORA_GOBP_DiffExp_Cortex.rds")

# log10p_cortex_vec <- diffExp_cortex_df$pvalue
# names(log10p_cortex_vec) <- diffExp_cortex_df$Gene
# log10p_cortex_vec[is.na(log10p_cortex_vec)] <- 1
# log10p_cortex_vec <- -log10(log10p_cortex_vec)
# 
# inp <- data.frame(V1 = diffExp_cortex_df$Gene,
#                   V2 = diffExp_cortex_df$log2FoldChange*log10p_cortex_vec[diffExp_cortex_df$Gene])
# gsea_gobp <- WebGestaltR(enrichMethod = "GSEA",
#                         organism = "mmusculus",
#                         enrichDatabase = "geneontology_Biological_Process",
#                         interestGene = inp,#diffExp_cortex_df$Gene[which(diffExp_cortex_df$padj < 0.05)],
#                         interestGeneType = "ensembl_gene_id",
#                         referenceGene = gsub("\\.[0-9]+","",rownames(inp_df)),
#                         referenceGeneType = "ensembl_gene_id",
#                         sigMethod = "fdr",
#                         fdrThr = 0.05,
#                         minNum = 3,
#                         maxNum = 10000,
#                         isOutput = F)
# 
# saveRDS(ora_gobp, file = "GSEA_GOBP_DiffExp_Cortex.rds")

out <- WebGestaltR(enrichMethod = "ORA",
                  organism = "mmusculus",
                  enrichDatabaseFile = "msigdb.v2024.1.Mm.symbols.gmt",
                  enrichDatabaseType = "genesymbol",
                  interestGene = diffExp_cortex_df$Gene[which(diffExp_cortex_df$padj < 0.05)],
                  interestGeneType = "ensembl_gene_id",
                  referenceGene = gsub("\\.[0-9]+","",rownames(inp_df)),
                  referenceGeneType = "ensembl_gene_id",
                  sigMethod = "fdr",
                  fdrThr = 0.05,
                  minNum = 3,
                  maxNum = 100000,
                  isOutput = F)

saveRDS(out,file = "ORA_MSigDBComplete_Cortex.rds")

out <- WebGestaltR(enrichMethod = "GSEA",
                   organism = "mmusculus",
                   enrichDatabaseFile = "AnxietyGeneSets_Mouse.gmt",
                   enrichDatabaseType = "genesymbol",
                   interestGene = inp,#diffExp_cortex_df$Gene[which(diffExp_cortex_df$padj < 0.05)],
                   interestGeneType = "ensembl_gene_id",
                   referenceGene = gsub("\\.[0-9]+","",rownames(inp_df)),
                   referenceGeneType = "ensembl_gene_id",
                   sigMethod = "top",
                   fdrThr = 0.05,
                   minNum = 3,
                   maxNum = 100000,
                   isOutput = F)

