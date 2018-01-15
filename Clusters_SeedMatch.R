#! usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/Scripts")

library(data.table)
library(splitstackshape)
library(stringr)

#######################################################################################################################################
### Specify directories

INDIRC <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/Annotated/"
INDIRM <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/SeedMatches/"
OUTDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/Targets/"

cfiles <- list.files(INDIRC, pattern = "\\.csv")
mfiles7 <- list.files(INDIRM, pattern = "_7mer\\.csv")
mfiles6 <- list.files(INDIRM, pattern = "_6mer\\.csv")

cfiles[23] <- cfiles[21] #Ensuring file names match between cfiles and mfiles, needs clean up
cfiles[24] <- cfiles[22]
cfiles[22] <- cfiles[20]
cfiles[21] <- cfiles[20]
cfiles[20] <- cfiles[19]

sample_names <- sub("*.csv", "", cfiles)

#######################################################################################################################################
### Define a function to merge

Merger <- function (mfile, cfile, filename) {

temp <- read.csv(paste0(INDIRM, mfile), stringsAsFactors = F)
temp$V5 <- paste0(temp$V3,"(",temp$V4,")")
temp1 <- aggregate(temp$V5 ~ temp$V1, data = temp, function(x) (paste(x, collapse=",")))
colnames(temp1) <- c("ClusterID", "miRNA")
temp2 <- aggregate(temp$V4 ~ temp$V1, data = temp, function(x) (sum(x)))
colnames(temp2) <- c("ClusterID", "SeedMatchesTotal")
temp3 <- merge(temp1,temp2, by= "ClusterID")

info <- read.csv(paste0(INDIRC, cfile),stringsAsFactor = F)
info <- merge(info, temp3, by= "ClusterID")
write.csv(info, paste0(OUTDIR,filename), row.names= F)

}

#######################################################################################################################################
### Run function over all 7mer and 6mer clusters

for (i in 1:length(cfiles)) {
  Merger(mfiles6[i], cfiles[i], mfiles6[i])
  Merger(mfiles7[i], cfiles[i], mfiles7[i])
}

#######################################################################################################################################
### Classify the RNA in four broad types: mRNA, smallRNA, lncRNA or unannotated and count each class

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
write.csv(counting, paste0(OUTDIR,"Counts/","Classes_Targets.csv"), row.names=F)

#######################################################################################################################################
### Ignore unannotated ones

# counting <- data.frame(matrix(NA, nrow=length(files), ncol= 4))
# 
# colnames(counting) <- c("SampleName", "mRNA", "smallRNA", "lncRNA")
# 
# count =1
# for(file in files) {
#   temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
#   temp <- na.omit(unique(temp[,c("transcript_id","transcript_class")]))
#   counting[count,1] <- file
#   counting[count,2] <- sum(temp$transcript_class == "mRNA")
#   counting[count,3] <- sum(temp$transcript_class == "smallRNA")
#   counting[count,4] <- sum(temp$transcript_class == "lncRNA")
#   count= count+1
# }
# write.csv(counting, paste0(OUTDIR,"Counts/","Classes_Targets_genes.csv"), row.names=F)

#######################################################################################################################################
### Obtain RNA class counts for genes (instead of clusters as done above)

files <- list.files(path = OUTDIR ,pattern = "\\.csv$")

counting <- data.frame(matrix(NA, nrow=length(files), ncol= 4))

colnames(counting) <- c("SampleName", "mRNA", "smallRNA", "lncRNA")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
  temp <- na.omit(unique(temp[,c("transcript_id","transcript_class")]))
  counting[count,1] <- file
  counting[count,2] <- sum(temp$transcript_class == "mRNA")
  counting[count,3] <- sum(temp$transcript_class == "smallRNA")
  counting[count,4] <- sum(temp$transcript_class == "lncRNA")
  count= count+1
}
write.csv(counting, paste0(OUTDIR,"Counts/","Classes_Targets_genes.csv"), row.names=F)

#######################################################################################################################################
