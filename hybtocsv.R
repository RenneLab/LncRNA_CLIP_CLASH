#! /usr/bin/env Rscript

rm(list=ls())

#Script to convert .hyb files into .csv files for using later in CLASH analysis

setwd("/ufrc/renne/sunantha.s/research/CLASH_lncRNAs/Scripts")

library(splitstackshape)
library(data.table)

###########################################################################################################################################
### Read all input files and create a counts file. The counts file will be used to count the number of cellular/viral miRNA hybrids formed with
### lncRNAs in both possible orientations.

INDIR ="../Hyb/"
OUTDIR="../ExcelFiles/"

files <- list.files(path = INDIR ,pattern = "\\hybrids_ua.hyb$")
sample_names <- sub("*_comp_lncRNA_kshv_miR_hybrids_ua.hyb", "", files)

names <- c("SeqID","ReadSequences","Info_5","Read_start_5","Read_end_5","Transcript_start_5",
           "Transcript_end_5","MappingScore_5","Info_3","Read_start_3","Read_end_3","Transcript_start_3",
           "Transcript_end_3","MappingScore_3")

Counts <- data.frame(matrix(NA, nrow=length(files), ncol=6))
colnames(Counts) <- c("Sample","Total","miRNA_lncRNA", "lncRNA_miRNA", "KSHV_lncRNA", "lncRNA_KSHV")
write.csv(Counts, paste0(OUTDIR,"Counts/","HybridCounts.csv"), row.names=F)

###########################################################################################################################################
### Define a function that will take the hyb file name and sample name as input to create four different files: cellular/ viral miRNA and 
### 5' ligated/ 3' ligated orientations. 

formatConvert <- function (file,samp) {
temp <- read.delim(file, header=F)
temp <- temp[,-c(3,16)]
colnames(temp) <- names
ind_5mi3m <- intersect(grep("hsa-",temp$Info_5), grep("ENSG", temp$Info_3))
ind_5m3mi <- intersect(grep("ENSG", temp$Info_5),grep("hsa-",temp$Info_3))
ind_5v3m <- intersect(grep("kshv-",temp$Info_5), grep("ENSG", temp$Info_3))
ind_5m3v <- intersect(grep("ENSG", temp$Info_5),grep("kshv-",temp$Info_3))

Counts <- read.csv(paste0(OUTDIR,"Counts/","HybridCounts.csv"), stringsAsFactors = F)
Counts$Sample[i] <- samp
Counts$Total[i] <- length(temp[,1])
Counts$miRNA_lncRNA[i] <- length(ind_5mi3m)
Counts$lncRNA_miRNA[i] <- length(ind_5m3mi)
Counts$KSHV_lncRNA[i] <- length(ind_5v3m)
Counts$lncRNA_KSHV[i] <- length(ind_5m3v)

if(length(ind_5mi3m) != 0) {
file_5mi3m <- cSplit(temp[ind_5mi3m,], 9, "_")
file_5mi3m <- subset.data.frame(file_5mi3m, select=c(1:16))
colnames(file_5mi3m)[16] <- "GeneName_3"
colnames(file_5mi3m)[15] <- "TranscriptID_3"
colnames(file_5mi3m)[14] <- "GeneID_3"
names(file_5mi3m)[names(file_5mi3m) == 'Info_5'] <- "GeneName_5"
write.csv(file_5mi3m , paste0(OUTDIR,samp,"_5mi3m.csv"), row.names=F)
}

if(length(ind_5v3m) != 0) {
file_5v3m <- cSplit(temp[ind_5v3m,], 9, "_")
file_5v3m <- subset.data.frame(file_5v3m, select=c(1:16))
colnames(file_5v3m)[16] <- "GeneName_3"
colnames(file_5v3m)[15] <- "TranscriptID_3"
colnames(file_5v3m)[14] <- "GeneID_3"
names(file_5v3m)[names(file_5v3m) == 'Info_5'] <- "GeneName_5"
write.csv(file_5v3m, paste0(OUTDIR,samp,"_5v3m.csv"), row.names=F)
}

if(length(ind_5m3mi) != 0) {
file_5m3mi <- cSplit(temp[ind_5m3mi,], 3, "_")
file_5m3mi <- subset.data.frame(file_5m3mi, select=c(1:16))
colnames(file_5m3mi)[16] <- "GeneName_5"
colnames(file_5m3mi)[15] <- "TranscriptID_5"
colnames(file_5m3mi)[14] <- "GeneID_5"
names(file_5m3mi)[names(file_5m3mi) == 'Info_3'] <- "GeneName_3"
write.csv(file_5m3mi , paste0(OUTDIR,samp,"_5m3mi.csv"), row.names=F)
}


if(length(ind_5m3v) != 0) {
file_5m3v <- cSplit(temp[ind_5m3v,], 3, "_")
file_5m3v <- subset.data.frame(file_5m3v, select=c(1:16))
colnames(file_5m3v)[16] <- "GeneName_5"
colnames(file_5m3v)[15] <- "TranscriptID_5"
colnames(file_5m3v)[14] <- "GeneID_5"
names(file_5m3v)[names(file_5m3v) == 'Info_3'] <- "GeneName_3"
write.csv(file_5m3v , paste0(OUTDIR,samp,"_5m3v.csv"), row.names=F)
}

write.csv(Counts, paste0(OUTDIR,"Counts/","HybridCounts.csv"), row.names=F)

}

############################################################################################################################################
### Loop function over all files and samples

for (i in 1:length(files)) {
  formatConvert(paste0(INDIR,files[i]), sample_names[i])
}

###########################################################################################################################################
