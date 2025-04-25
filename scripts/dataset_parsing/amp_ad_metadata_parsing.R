################################################################################
# Brain clock: go through the individual metadata files of each one of the     #
# studies to determine which individuals are controls and which not and unify  #
# the information included in them.                                            #
################################################################################

setwd("../..")
outDir <- "./results/metadata_parsed/"

if (!dir.exists(outDir)){
    dir.create(outDir, recursive = T)
}

# LBP
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520

LBP <- "./data/metadata/LBP"
LBPInd <- read.csv(paste(LBP, "LBP_individual_metadata.csv", sep = "/"))

isCtrl <- LBPInd$diagnosis == "control"

LBPInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrl) + 1]

# Remove 89+
LBPInd <- LBPInd[LBPInd$ageDeath != "89+", ]
LBPInd$ageDeath <- as.numeric(LBPInd$ageDeath)
LBPInd$race <- gsub("Black or African American", "Black", LBPInd$race)

LBPIndCtrl <- LBPInd[LBPInd$diagn_4BrainClck == "Control", ]

write.csv(LBPInd, paste0(outDir, "LBP_ind_all.csv"))

# RNAseq_Harmonization
################################################################################

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study_Name=(RNAseq_Harmonization)

# This is a harmonized version of ROSMAP, MSBB and Mayo, so samples should be 
# the same, but since they have been harmonized it's better to use these
# metadata.

RNAseqHarm <- "./data/metadata/RNAseq_Harmonization"

MayoHarmInd <- read.csv(paste(RNAseqHarm,
                              "/RNAseq_Harmonization_Mayo_combined_metadata.csv",
                              sep = "/"))

# In the RNAseq Harmonization Mayo dataset there are several samples that
# don't have Braak score or Thal score, while in the original Mayo metadata
# they have. Replace them and filter based on that.

# Downloaded from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn9702085

rnaSeqReprocessing <- "./data/metadata/rnaSeqReprocessing"
rnaSeqReprocessing_files <- list.files(rnaSeqReprocessing, full.names = T)

mayoIndiv <- read.csv(paste(rnaSeqReprocessing,
                            "MayoRNAseq_individual_metadata.csv",
                            sep = "/"))

mayoIndivCtrls <- mayoIndiv[mayoIndiv$diagnosis == "control" & mayoIndiv$Thal == 0, ]

mayoIndivCtrls <- mayoIndivCtrls[!is.na(mayoIndivCtrls$diagnosis), ]

mayoIndivCtrls <- mayoIndivCtrls[!is.na(mayoIndivCtrls$individualID), ]

MayoHarmInd$Braak[is.na(MayoHarmInd$Braak)] <- mayoIndiv$Braak[match(MayoHarmInd$individualID[is.na(MayoHarmInd$Braak)],
                                                                     mayoIndiv$individualID)]

MayoHarmInd$thal[is.na(MayoHarmInd$thal)] <- mayoIndiv$Thal[match(MayoHarmInd$individualID[is.na(MayoHarmInd$thal)],
                                                                  mayoIndiv$individualID)]


isCtrlMayo <- MayoHarmInd$diagnosis == "control" & MayoHarmInd$Braak <= 3 & MayoHarmInd$thal == 0

isCtrlMayo[is.na(isCtrlMayo)] <- F


MayoHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMayo) + 1]

# Remove ambiguous ages and convert to numeric
MayoHarmInd <- MayoHarmInd[MayoHarmInd$ageDeath != "90+" & MayoHarmInd$ageDeath != "", ]

MayoHarmInd$ageDeath <- as.numeric(MayoHarmInd$ageDeath)

MayoHarmIndCtrl <- MayoHarmInd[MayoHarmInd$diagn_4BrainClck == "Control", ]

# The duplicated ones is because there is a cerebellum sample for that 
# individual. Probably we will need to remove it because it is a very different
# tissue, but we'll see.

MSBBHarmInd <- read.csv(paste(RNAseqHarm,
                              "/RNAseq_Harmonization_MSBB_combined_metadata.csv",
                              sep = "/"))

isCtrlMSBB <- MSBBHarmInd$Braak <= 3 & MSBBHarmInd$CERAD == 4
MSBBHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlMSBB) + 1]
MSBBHarmInd <- MSBBHarmInd[!is.na(MSBBHarmInd$diagn_4BrainClck), ]

# Remove ambiguous ages and convert to numeric
MSBBHarmInd <- MSBBHarmInd[MSBBHarmInd$ageDeath != "90+", ]
MSBBHarmInd$ageDeath <- as.numeric(MSBBHarmInd$ageDeath)

MSBBHarmIndCtrl <- MSBBHarmInd[MSBBHarmInd$diagn_4BrainClck == "Control", ]


ROSMAPHarmInd <- read.csv(paste(RNAseqHarm,
                              "/RNAseq_Harmonization_ROSMAP_combined_metadata.csv",
                              sep = "/"))

isCtrlROSMAP <- ROSMAPHarmInd$dcfdx_lv == 1 & ROSMAPHarmInd$cogdx == 1
isCtrlROSMAP <- isCtrlROSMAP & ROSMAPHarmInd$braaksc <= 3 & ROSMAPHarmInd$ceradsc == 4



ROSMAPHarmInd$diagn_4BrainClck <- c("Rest", "Control")[as.numeric(isCtrlROSMAP) + 1]
ROSMAPHarmInd <- ROSMAPHarmInd[!is.na(ROSMAPHarmInd$diagn_4BrainClck), ]

# Remove ambiguous ages and, convert to numeric and round
ROSMAPHarmInd <- ROSMAPHarmInd[ROSMAPHarmInd$age_death != "90+", ]
ROSMAPHarmInd$age_death <- round(as.numeric(ROSMAPHarmInd$age_death))


ROSMAPHarmIndCtrl <- ROSMAPHarmInd[ROSMAPHarmInd$diagn_4BrainClck == "Control", ]


# Save the datasets in separated files.
write.csv(MayoHarmInd, file = paste0(outDir, "Mayo_Harm_parsed.csv"))
write.csv(MSBBHarmInd, file = paste0(outDir, "MSBB_Harm_parsed.csv"))
write.csv(ROSMAPHarmInd, file = paste0(outDir, "ROSMAP_Harm_parsed.csv"))


# Combine the three of them in a single dataframe

# Find the columns with common names and change the colnames with discrepancies
commCols <- colnames(MayoHarmIndCtrl)[colnames(MayoHarmInd) %in% colnames(MSBBHarmInd) & colnames(MayoHarmInd) %in% colnames(ROSMAPHarmInd)]

colnames(ROSMAPHarmInd) <- gsub("apoe_genotype",
                                "apoeGenotype",
                                colnames(ROSMAPHarmInd))

ROSMAPHarmInd$msex <- c("female", "male")[ROSMAPHarmInd$msex + 1]
colnames(ROSMAPHarmInd) <- gsub("msex",
                                "sex",
                                colnames(ROSMAPHarmInd))

colnames(ROSMAPHarmInd) <- gsub("braaksc",
                                "Braak",
                                colnames(ROSMAPHarmInd))

colnames(ROSMAPHarmInd) <- gsub("age_death",
                                "ageDeath",
                                colnames(ROSMAPHarmInd))

ROSMAPHarmInd$race <- c("White",
                        "Black",
                        "Native_American",
                        "Native_Hawaiian",
                        "Asian",
                        "Other",
                        "Unknown")[ROSMAPHarmInd$race]

ROSMAPHarmInd$race[ROSMAPHarmInd$spanish == 1] <- "Hispanic"

MSBBHarmInd$race <- gsub("W", "White", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("B", "Black", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("H", "Hispanic", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("A", "Asian", MSBBHarmInd$race)
MSBBHarmInd$race <- gsub("U", "Unknown", MSBBHarmInd$race)

# Reassign commCols to account for the unifying changes and check that the
# format of the variables in the columns with common name is the same
commCols <- colnames(MayoHarmInd)[colnames(MayoHarmInd) %in% colnames(MSBBHarmInd) & colnames(MayoHarmInd) %in% colnames(ROSMAPHarmInd)]

# Add a column of substudy and merge datasets.
ROSMAPHarmInd$substudy <- rep("ROSMAP", nrow(ROSMAPHarmInd))
MayoHarmInd$substudy <- rep("Mayo", nrow(MayoHarmInd))
MSBBHarmInd$substudy <- rep("MSBB", nrow(MSBBHarmInd))



commCols <- colnames(MayoHarmInd)[colnames(MayoHarmInd) %in%
                                          colnames(MSBBHarmInd) &
                                          colnames(MayoHarmInd) %in%
                                          colnames(ROSMAPHarmInd)]

RNAseqHarmDF <- rbind.data.frame(ROSMAPHarmInd[, commCols],
                                 MayoHarmInd[, commCols],
                                 MSBBHarmInd[, commCols])

write.csv(RNAseqHarmDF, file = paste0(outDir, "RNAseq_Harmonization_ind_all.csv"))