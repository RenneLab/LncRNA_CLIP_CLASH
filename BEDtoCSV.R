#! usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/Scripts")

library(data.table)
library(splitstackshape)

###############################################################################################################################
### Add annotations to the cluster BED files

INDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/"
OUTDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/Annotated/"
  
SeqFiles <- list.files(path = INDIR ,pattern = "\\_seqs.bed$")
AnnFiles <- list.files(path = INDIR ,pattern = "\\_ann.bed$")

sample_names <- sub("*_ann.bed", "", AnnFiles)

names <- c("ClusterID", "Sequence", "Chr", "Start", "Stop", "Read_Count", "Strand", "p_value", "FDR", "Lnc_Chr", "Lnc_Start", 
           "Lnc_End", "Lnc_Name", "Lnc_Strand", "Source", "Type", "transcript_class", "Overlap", "gene_id","transcript_id", "gene_type", 
           "gene_status", "gene_name", "transcript_type", "transcript_status", "transcript_name")

for (i in 1: length(SeqFiles)) {
  seq <- read.delim(paste0(INDIR,SeqFiles[i]), stringsAsFactors = F, header = FALSE)
  colnames(seq) <- c("ClusterID","Sequence")
  seq$Sequence <- toupper(seq$Sequence)
  ann <- read.delim(paste0(INDIR,AnnFiles[i]), stringsAsFactors = F, header = FALSE)
  names(ann)[names(ann) == "V4"] <- "ClusterID"
  ann_seq <- merge(seq,ann, by="ClusterID")
  ann_seq <- cSplit(ann_seq, 19, ";")
  ann_seq <- subset.data.frame(ann_seq, select = -c(14,18,29:ncol(ann_seq)), use.names =F, fill=F)
  ann_seq$V18_01 <- gsub("gene_id ", "", ann_seq$V18_01)
  ann_seq$V18_02 <- gsub("transcript_id ", "", ann_seq$V18_02)
  ann_seq$V18_03 <- gsub("gene_type ", "", ann_seq$V18_03)
  ann_seq$V18_04 <- gsub("gene_status ", "", ann_seq$V18_04)
  ann_seq$V18_05 <- gsub("gene_name ", "", ann_seq$V18_05)
  ann_seq$V18_06 <- gsub("transcript_type ", "", ann_seq$V18_06)
  ann_seq$V18_07 <- gsub("transcript_status ", "", ann_seq$V18_07)
  ann_seq$V18_08 <- gsub("transcript_name ", "", ann_seq$V18_08)
  colnames(ann_seq) <- names
  ann_seq$Read_Length <- nchar(ann_seq$Sequence)
  ann_seq <- ann_seq[ann_seq$Read_Count>4,]
  write.csv(ann_seq, paste0(OUTDIR,sample_names[i],".csv"), row.names = FALSE)
}

###############################################################################################################################
### Count the number of clusters in different RNA categories

files <- list.files(OUTDIR, pattern = "\\.csv$")

counting <- data.frame(matrix(NA, nrow=length(files), ncol= 5))

colnames(counting) <- c("SampleName", "Unannotated", "mRNA", "smallRNA", "lncRNA")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
  counting[count,1] <- file
  counting[count,2] <- sum(temp$transcript_class == ".")
  counting[count,3] <- sum(temp$transcript_class == "mRNA")
  counting[count,4] <- sum(temp$transcript_class == "smallRNA")
  counting[count,5] <- sum(temp$transcript_class == "lncRNA")
  count= count+1
}
write.csv(counting, paste0(OUTDIR,"Counts/","Classes.csv"), row.names=F)

###############################################################################################################################