#! /usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLASH_lncRNAs/Scripts")

###########################################################################################################################################
#Read all CLASH data files with 5'miRNA-3'mRNA chimeras and assign column names

INDIR <- "../ExcelFiles/"
OUTDIR <- "../Binding/"

lncRNA_length <- read.table("/ufrc/renne/sunantha.s/research/Ensembl/Human_ncRNA/Human_ncRNA_lengths.fa", sep="|", stringsAsFactors = F)
lncRNA_length <- lncRNA_length[,c(2,4)]
colnames(lncRNA_length) <- c("Transcript_ID", "Transcript_length")

###########################################################################################################################################
#Change the comparison column from factor to character (**check why**) and if loop to classify each mRNA fragment into 6 classes: 5'UTR,
#5'UTR-CDS, "CDS", "CDS-3'UTR" and "unknown". Write the output files to a different folder.

files <- list.files(path = INDIR ,pattern = "\\.csv$")

files_5miR <- c(grep("5mi3m", files, value=T), grep("5v3m", files, value=T))
files_3miR <- c(grep("5m3mi", files, value=T), grep("5m3v", files, value=T))

for (file in files_5miR) {
  
filename <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
filename$TranscriptID_3 <- as.character(filename$TranscriptID_3)

filename <- merge(filename, lncRNA_length, by.x= "TranscriptID_3", by.y= "Transcript_ID", all.x=T)
  for (i in 1:nrow(filename)) {
    if(is.na(filename$Transcript_start_3[i]) || is.na(filename$Transcript_length[i])) {
      filename$Binding_Region[i] <- "unknown"
    } else if (filename$Transcript_start_3[i] < (1/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "A"
    } else if (filename$Transcript_start_3[i] > (1/5)*filename$Transcript_length[i] && filename$Transcript_start_3[i] < (2/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "B"
    } else if (filename$Transcript_start_3[i] > (2/5)*filename$Transcript_length[i] && filename$Transcript_start_3[i] < (3/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "C"
    } else if (filename$Transcript_start_3[i] > (3/5)*filename$Transcript_length[i] && filename$Transcript_start_3[i] < (4/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "D"
    } else if (filename$Transcript_start_3[i] > (4/5)*filename$Transcript_length[i] && filename$Transcript_start_3[i] < (5/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "E"
    }
  }
  write.csv(filename, paste0(OUTDIR, file), row.names=F)
}

for (file in files_3miR)  {
  
  filename <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  filename$TranscriptID_5 <- as.character(filename$TranscriptID_5)
  
  filename <- merge(filename, lncRNA_length, by.x= "TranscriptID_5", by.y= "Transcript_ID", all.x=T)
  for (i in 1:nrow(filename)) {
    if(is.na(filename$Transcript_start_5[i]) || is.na(filename$Transcript_length[i])) {
      filename$Binding_Region[i] <- "unknown"
    } else if (filename$Transcript_start_5[i] < (1/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "A"
    } else if (filename$Transcript_start_5[i] > (1/5)*filename$Transcript_length[i] && filename$Transcript_start_5[i] < (2/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "B"
    } else if (filename$Transcript_start_5[i] > (2/5)*filename$Transcript_length[i] && filename$Transcript_start_5[i] < (3/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "C"
    } else if (filename$Transcript_start_5[i] > (3/5)*filename$Transcript_length[i] && filename$Transcript_start_5[i] < (4/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "D"
    } else if (filename$Transcript_start_5[i] > (4/5)*filename$Transcript_length[i] && filename$Transcript_start_5[i] < (5/5)*filename$Transcript_length[i]) {
      filename$Binding_Region[i] <- "E"
    }
  }
  write.csv(filename, paste0(OUTDIR, file), row.names=F)
}

###########################################################################################################################################
### Count the number of binding events in each bin for each sample analyzed

files <- list.files(OUTDIR)
files <- files[-1]

counting <- data.frame(matrix(NA, nrow=length(files), ncol=6))
colnames(counting) <- c("filename","20%","40%", "60%", "80%", "100%")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file))
  counting[count,1] <- file
  counting[count,2] <- sum(temp$Binding_Region == "A")
  counting[count,3] <- sum(temp$Binding_Region == "B")
  counting[count,4] <- sum(temp$Binding_Region == "C")
  counting[count,5] <- sum(temp$Binding_Region == "D")
  counting[count,6] <- sum(temp$Binding_Region == "E")
  count= count+1
}
write.csv(counting, paste0(OUTDIR,"Counts/","BindingRegion.csv"), row.names=F)

###########################################################################################################################################