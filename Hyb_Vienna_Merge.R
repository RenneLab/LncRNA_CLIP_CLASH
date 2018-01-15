#! /usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLASH_lncRNAs/Scripts")

library(splitstackshape)
library(dplyr)
library(data.table)

###########################################################################################################################################
### Specify the input hyb and vienna files

INDIRV <- "../Hyb/"
INDIRH <- "../ExcelFiles/"
OUTDIR <- "../HybVienna/"

vfiles <- list.files(path = INDIRV ,pattern = "\\hybrids_ua.viennad.line$")
sample_names <- sub("*_comp_lncRNA_kshv_miR_hybrids_ua.viennad.line", "", vfiles)
classes <- c("5mi3m", "5v3m", "5m3mi", "5m3v")

###########################################################################################################################################
### Function to clean up and merge hyb and vienna files 

Merger <- function(file,samp){

vienna <- read.delim(file,header=F, stringsAsFactors = F)
kshv <- c(grep("kshv-", vienna$V1))
vienna <- cSplit(vienna,1,"_")
vienna <- select(vienna, c(12,1,2,6,8,9,10,11))
vienna <- rbindlist(list(subset.data.frame(vienna[-kshv,], select=c(1,2,3,4,5,6)),
              subset.data.frame(vienna[kshv,], select=c(1,2,3,4,5,6))), use.names=F, fill=F)
names <- c("SeqID", "Seq_diag", "Seq_diag_5", "Seq_diag_3", "Diagram", "Folding_Energy")
colnames(vienna) <- names

for (class in classes) {

if(file.exists(paste0(INDIRH,samp,"_",class,".csv"))) {
  
hyb <- read.csv(paste0(INDIRH,samp,"_",class,".csv"), stringsAsFactors = F)
hyb <- cSplit(hyb, 1, "_")
names(hyb)[names(hyb) == 'SeqID_1'] <- "SeqID"
names(hyb)[names(hyb) == 'SeqID_2'] <- "Num_PCR_Dups"

hyb_vienna <- merge(as.data.frame(vienna), as.data.frame(hyb), by="SeqID")

if (class == "5mi3m" | class == "5v3m") {
  hyb_vienna$miRNA_length <- hyb_vienna$Read_end_5-hyb_vienna$Read_start_5+1
  hyb_vienna$miR_Diagram <- substr(hyb_vienna$Diagram, start=1,stop = hyb_vienna$miRNA_length)
  }
if (class == "5m3mi" | class == "5m3v") {
  hyb_vienna$miRNA_length <- hyb_vienna$Read_end_3-hyb_vienna$Read_start_3+1
  hyb_vienna$miR_Diagram <- substr(hyb_vienna$Diagram, start= hyb_vienna$Read_start_3-hyb_vienna$Read_start_5+1, stop = nchar(hyb_vienna$Diagram))
  }

write.csv(hyb_vienna, paste0(OUTDIR,samp,"_hyb_vienna_",class,".csv"))

}
}
}

###########################################################################################################################################
### Run the function over all hyb and vienne file pairs

for (i in 1:length(vfiles)) {
  Merger(paste0(INDIRV,vfiles[i]), sample_names[i])
}

###########################################################################################################################################