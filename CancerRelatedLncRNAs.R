#! usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/Scripts")

library(data.table)
library(splitstackshape)
library(stringr)

INDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/Pooled/"
REFDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/CancerAssociated/"

#######################################################################################################################################
### Use Ensembl data to add gene ID information to cancer associated lncRNAs from databases

Ensembl_data <- read.csv(paste0(REFDIR,"mart_export.txt"), stringsAsFactors = F)

Cancer <- read.csv(paste0(REFDIR,"lnc_Cancer.csv"), stringsAsFactors = F)
Cancer_lnc <- Cancer[,c(4,7)]
Cancer_lnc <- Cancer_lnc[!(Cancer_lnc$X.1 == "N/A"), c(2,1)]
colnames(Cancer_lnc) <- c("Cancer", "GeneID")


Disease <- read.csv(paste0(REFDIR,"lnc_Diseases.csv"), stringsAsFactors = F, header=F)
Disease_lnc <- Disease[,c(12,3)]
Disease_lnc <- cSplit(Disease_lnc, sep=".", splitCols = 1)
cancer <- grep("cancer|Cancer|tumor|Tumor|oma|malignan", Disease_lnc$V3)
temp <- Disease_lnc[cancer,]
glaucoma <- grep("glaucoma|Glaucoma", temp$V3)
Disease_lnc <- temp[!glaucoma,]
Disease_lnc$GeneID <- NA

for(i in 1:nrow(Disease_lnc)) {
  index <- c(grep(Disease_lnc$V12_1[i], Ensembl_data$RefSeq.ncRNA.ID),grep(Disease_lnc$V12_1[i], Ensembl_data$RefSeq.ncRNA.predicted.ID))
  if(length(index) > 0) {
    Disease_lnc$GeneID[i] <- Ensembl_data$Gene.stable.ID[index]
    }
}
Disease_lnc <- Disease_lnc[!is.na(Disease_lnc$GeneID),]
Disease_lnc <- subset.data.frame(Disease_lnc, select=c("V3", "GeneID"))
colnames(Disease_lnc) <- c("Cancer", "GeneID")

Cancer_db <- rbind(Cancer_lnc,Disease_lnc)
Cancer_db$Cancer <- tolower(Cancer_db$Cancer)
Cancer_db <- unique(Cancer_db)
Cancer_db <-aggregate(Cancer_db$Cancer ~ Cancer_db$GeneID, data = Cancer_db, function(x) (paste(x, collapse=",")))
colnames(Cancer_db) <- c("GeneID", "Cancer")

write.csv(Cancer_db, paste0(REFDIR,"Cancer_lncRNAs.csv"), row.names = F)

#######################################################################################################################################
### Compare cancer database information with targets identified from CLIP-Seqs

Cancer_db <- read.csv(paste0(REFDIR,"Cancer_lncRNAs.csv"), stringsAsFactors = F)

files <- list.files(INDIR, pattern="E.*\\.csv|K.*\\.csv")
datalist <- list()

counting <- data.frame(matrix(NA, nrow=length(files), ncol= 2))

colnames(counting) <- c("SampleName", "lncRNA")


count =1

for (file in files) {
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  temp$GeneID <- temp$gene_id
  temp <- cSplit(temp, sep=".", splitCols = "GeneID")
  temp <- subset.data.frame(temp,select=c(1:31))
  colnames(temp)[31] <- "GeneID"
  temp <- merge(temp,Cancer_db, by="GeneID")
  write.csv(temp,paste0(REFDIR,file), row.names = F)
  datalist[[count]] <- temp
  lncRNA <- unique(subset.data.frame(temp, temp$transcript_class=="lncRNA", select = "gene_id"))
  counting[count,1] <- file
  counting[count,2] <- nrow(lncRNA)
  count =count+1
}

counting[count,1] <- "Unique"
counting[count,2] <- length(unique(df$Lnc_Name))

df <- do.call(rbind,datalist)
write.csv(df, paste0(REFDIR,"AllTargets_CancerAssociated.csv"), row.names=F)

write.csv(counting, paste0(REFDIR,"Counts.csv"), row.names = F)
#######################################################################################################################################

