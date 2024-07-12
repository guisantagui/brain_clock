# Script for debugging the rownames that appear unmapped to
# the metdat in the concat files

parsedDir <- "/mnt/lscratch/users/gsantamaria/test_large_files/NPC/parsed_mats/"

parsedFiles <- list.files(parsedDir, full.names = T)

metDatFiles <- parsedFiles[grep("metDat", parsedFiles)]
metDatFiles <- metDatFiles[!grepl("concat_metDat", metDatFiles)]

metDat_concat <- parsedFiles[grep("concat_metDat", parsedFiles)]

expMatFiles <- parsedFiles[grep("expMat", parsedFiles)]
expMatFiles <- expMatFiles[!grepl("concat_expMat", expMatFiles)]

expMat_concat <- parsedFiles[grep("concat_expMat", parsedFiles)]

# read.csv but faster
readCsvFst <- function(pth, header = T){
        df <- data.frame(data.table::fread(pth, header = header))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

readExpMetDatByIdx <- function(i){
        mDat <- readCsvFst(metDatFiles[i])
        eMat <- readCsvFst(expMatFiles[i])
        return(list(mDat = mDat, eMat = eMat))
}

mDatConc <- readCsvFst(metDat_concat)
eMatConc <- readCsvFst(expMat_concat)

for(i in 1:length(expMatFiles)){
        mdFile <- basename(metDatFiles[i])
        xmFile <- basename(expMatFiles[i])
        print(sprintf("%s vs %s", mdFile, xmFile))
        out <- readExpMetDatByIdx(i)
        print(sprintf("All samples in %s are in %s:", xmFile, mdFile))
        print(all(colnames(out$eMat) %in% make.names(out$mDat$specimenID)))
        
        print(sprintf("All samples in %s are in conc_expMat:", xmFile))
        print(all(colnames(out$eMat) %in% colnames(eMatConc)))
        
        print(sprintf("All samples in %s are in conc_metDat:", xmFile))
        print(all(colnames(out$eMat) %in% make.names(mDatConc$specimenID)))
        
        print(sprintf("All samples in %s are in conc_expMat:", mdFile))
        print(all(make.names(out$mDat$specimenID) %in% colnames(eMatConc)))

        print(sprintf("All samples in %s are in conc_metDat:", mdFile))
        print(all(make.names(out$mDat$specimenID) %in% make.names(mDatConc$specimenID)))
}




all(colnames(eMatConc) %in% make.names(mDatConc$specimenID))
sum(!colnames(eMatConc) %in% make.names(mDatConc$specimenID))

eMatGlob <- readCsvFst("/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_counts_wTBI_wPert111_log2_quantNorm_preproc_wLINCS.csv")
mDatGlob <- readCsvFst("/home/users/gsantamaria/projects/brain_clock/data/int_database_w111/combined_metDat_wTBI_wPert111_wLINCS.csv")

sum(!rownames(eMatGlob) %in% make.names(mDatGlob$specimenID))