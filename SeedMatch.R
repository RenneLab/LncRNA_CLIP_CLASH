#! usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/Scripts")

library(data.table)
library(splitstackshape)
library(stringr)

#######################################################################################################################################
### Specify directories and group files by the miRNAs expressed in those cell types

INDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/Annotated/"
OUTDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/SeedMatches/"

EBV_files_EBV <- c(list.files(INDIR, pattern = "\\EH.*.csv"),"KP_BC1_1.csv", "KP_BC1_2.csv")
KSHV_files_KSHV <- c("KH_BCBL1_1.csv","KH_BCBL1_2.csv","KH_BCBL1_3.csv","KH_BCBL1_4.csv","KP_BC1_1.csv", "KP_BC1_2.csv")
EBV_files_B95_8 <- list.files(INDIR, pattern = "\\EP.*.csv")
KSHV_files_BC3 <- c("KH_BC3_1.csv", "KH_BC3_2.csv", "KH_BC3_3.csv", "KP_BC3_1.csv", "KP_BC3_2.csv")

#######################################################################################################################################
### Specify function for miRNA seed match identification within cluster sequences

SeedMatchFinder <- function(RNA_df,miRNA_df, Seedlength) {
  temp <- data.frame()
  count <- 1
  if(Seedlength == 7) {column <- 5}
  if(Seedlength == 6) {column <- 6}
  for (i in 1:length(RNA_df[,2])) {
    for (j in 1:length(miRNA_df[,column])) {
      if(str_count(RNA_df[i,2], miRNA_df[j,column]) != 0) {
        temp[count,1] <- RNA_df[i,1]
        temp[count,2] <- RNA_df[i,2]
        temp[count,3] <- miRNA_df[j,2]
        temp[count,4] <- str_count(RNA_df[i,2], miRNA_df[j,column])
        count <- count+1
      }        
    }
  }
  return(unique(temp))
}

#######################################################################################################################################
### Loop over every group (based on miRNAs expressed in those cell types) and execute the SeedMatchFinder function.

for(file in EBV_files_EBV) {
  filename <- sub(".csv", "", file)
  RNA_df <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  miRNA_df <- read.csv("../miRNAs/EBV_miRNAs.csv", stringsAsFactors = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,7)
  write.csv(temp, paste0(OUTDIR,filename, "_EBV_7mer.csv"), row.names = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,6)
  write.csv(temp, paste0(OUTDIR,filename, "_EBV_6mer.csv"), row.names = F)
}

for(file in KSHV_files_KSHV) {
  filename <- sub(".csv", "", file)
  RNA_df <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  miRNA_df <- read.csv("../miRNAs/KSHV_miRNAs.csv", stringsAsFactors = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,7)
  write.csv(temp, paste0(OUTDIR,filename, "_KSHV_7mer.csv"), row.names = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,6)
  write.csv(temp, paste0(OUTDIR,filename, "_KSHV_6mer.csv"), row.names = F)
}

for(file in EBV_files_B95_8) {
  filename <- sub(".csv", "", file)
  RNA_df <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  miRNA_df <- read.csv("../miRNAs/B95_8.csv", stringsAsFactors = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,7)
  write.csv(temp, paste0(OUTDIR,filename, "_B95_8_7mer.csv"), row.names = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,6)
  write.csv(temp, paste0(OUTDIR,filename, "_B95_8_6mer.csv"), row.names = F)
}

for(file in KSHV_files_BC3) {
  filename <- sub(".csv", "", file)
  RNA_df <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  miRNA_df <- read.csv("../miRNAs/BC3.csv", stringsAsFactors = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,7)
  write.csv(temp, paste0(OUTDIR,filename, "_BC3_7mer.csv"), row.names = F)
  temp <- SeedMatchFinder(RNA_df,miRNA_df,6)
  write.csv(temp, paste0(OUTDIR,filename, "_BC3_6mer.csv"), row.names = F)
}
#######################################################################################################################################
