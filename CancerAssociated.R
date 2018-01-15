#! /usr/bin/env Rscript

rm(list=ls())

library(splitstackshape)

setwd("/ufrc/renne/sunantha.s/research/CLASH_lncRNAs/Scripts")

###############################################################################################################################################
### Specify directories

INDIR <- "../ExcelFiles/"
REFDIR <- "../D11_Targets/"
OUTDIR <- "../CancerAssociated/"

file_5miR <- c("WT_BR1_5v3m.csv","WT_BR2_5v3m.csv","WT_BR3_5v3m.csv", "D11_BR1_5v3m.csv", "D11_BR2_5v3m.csv","D11_BR3_5v3m.csv")
file_3miR <- c("WT_BR1_5m3v.csv","WT_BR2_5m3v.csv","WT_BR3_5m3v.csv", "D11_BR1_5m3v.csv", "D11_BR2_5m3v.csv","D11_BR3_5m3v.csv")

###############################################################################################################################################
### Read all CLASH data files with 5'miRNA-3'lncRNA chimeras and assign column names

datafiles_5miR <- list()
i <- 1
for(file in file_5miR) {
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  temp <- temp[,c("GeneID_3","TranscriptID_3", "GeneName_3","GeneName_5")]
  temp <- cSplit(temp, sep="_", splitCols = "GeneName_5")
  temp <- subset.data.frame(temp, select=c(1:4))
  datafiles_5miR[[i]] <- temp
  i <- i+1
}
Targets_5miR <- do.call(rbind,datafiles_5miR)
colnames(Targets_5miR) <- c("GeneID", "TranscriptID", "GeneName", "miRNA")

###############################################################################################################################################
### Read all CLASH data files with 5'lncRNA-3'miRNA chimeras and assign column names

datafiles_3miR <- list()
i <- 1
for(file in file_3miR) {
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  temp <- temp[,c("GeneID_5","TranscriptID_5", "GeneName_5","GeneName_3")]
  temp <- cSplit(temp, sep="_", splitCols = "GeneName_3")
  temp <- subset.data.frame(temp, select=c(1:4))
  datafiles_3miR[[i]] <- temp
  i <- i+1
}
Targets_3miR <- do.call(rbind,datafiles_3miR)
colnames(Targets_3miR) <- c("GeneID", "TranscriptID", "GeneName", "miRNA")

###############################################################################################################################################
### Compare gene ids between targets and cancer databases

All_Targets <- unique(rbind(Targets_5miR, Targets_3miR))
Cancer <- read.csv(paste0(REFDIR,"Cancer_lncRNAs.csv"), stringsAsFactors = F)
Cancer_targets <- merge(All_Targets,Cancer, by.x="GeneID", by.y="GeneID")

write.csv(All_Targets, paste0(OUTDIR,"All_KSHV_Targets.csv"), row.names = F)
write.csv(Cancer_targets, paste0(OUTDIR, "KSHV_Cancer_lncRNAs.csv"), row.names = F)

###############################################################################################################################################
