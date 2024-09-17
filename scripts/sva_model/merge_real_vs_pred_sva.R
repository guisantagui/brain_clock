
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

predSVA_file <- "/home/users/gsantamaria/projects/brain_clock/results/svaMod/allButLINCS_signChron/LINCS_svaPredCorr.csv"
realSVA_file <- "/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integ_LINCSSamps_wSC_all_sva_fast_excldLINCS_filtSignChron/combined_counts_wTBI_wPert111_wSC_log2_quantNorm_preproc_modFuncsAlpha1_coefs_noCerebell_onlyAge_svaAdj.csv"
outFile <- "/home/users/gsantamaria/projects/brain_clock/results/preprocessing/integRealSVA_LINCSpredSVA/integRealSVA_LINCSpredSVA.csv"
outDir <- dirname(outFile)

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}


predSVA <- readCsvFast(predSVA_file)
realSVA <- readCsvFast(realSVA_file)

predSVA <- predSVA[, match(colnames(realSVA), colnames(predSVA))]

sva_merged <- rbind.data.frame(realSVA, predSVA)

writeCsvFst(sva_merged, file = outFile)