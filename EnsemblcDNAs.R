#! /usr/bin/env Rscript

rm(list=ls())

#R script to identify the part of mRNA (5'UTR, CDS or 3'UTR) that the miRNA binds to, using CLASH data.
#Author: Sunantha Sethuraman
#Date: 2016-06-01

setwd("/lfs/scratch/sunantha.s/research/Ensembl/HumanTrancriptome")

###############################################################################################################################################
# Reading the data downloaded from Ensembl into a dataframe

Human_mRNAs <- read.table("mRNAsList.txt", sep = "|", stringsAsFactors =F, header = F)

###############################################################################################################################################
### Loading necessary packages

library(splitstackshape)
library(dplyr)

###############################################################################################################################################
### Pre-processing Ensembl data: splitting exon boundary information to columns

Human_mRNAs <- cSplit(Human_mRNAs, 4, ";")
Exon_starts <- select(Human_mRNAs, 16:377)
Human_mRNAs <- select(Human_mRNAs,1:15)

Human_mRNAs <- cSplit(Human_mRNAs, 4, ";")
Exon_ends <- select(Human_mRNAs, 15:376)
Human_mRNAs <- select(Human_mRNAs,1:14)

###############################################################################################################################################
### Looping through columns to identify CDS start and CDS end (lowest exon bounday was considered CDS start and highest exon boundary was 
### considered CDS end)

for (i in 1:nrow(Exon_starts)){
  row <- as.numeric(Exon_starts[i,])
  row <- row[!is.na(row)]
  Human_mRNAs$CDS_start[i] <- min(row)
}

for (i in 1:nrow(Exon_ends)){
  row <- as.numeric(Exon_ends[i,])
  row <- row[!is.na(row)]
  Human_mRNAs$CDS_end[i] <- max(row)
}

#Name the Ensembl dataframe columns and write to file "Human_mRNAs.csv"; Only columns necessary for comparison were subset to Human_CDS dataframe

names <- c("Gene_ID","Transcript_ID", "Gene_Name", "5UTR_Start", "5UTR_End", 
           "3UTR_Start", "3UTR_End", "Transcript_Type", "Peptide_ID", "Strand",
           "Transcript_Start", "Transcript_End", "Transcription_StartSite", 
           "cDNA_Length", "CDS_Start", "CDS_End")
colnames(Human_mRNAs) <- names

Human_mRNAs$Gene_ID <- sub(">(ENSG00000[0-9]+)", "\\1", Human_mRNAs[,1])
Human_CDS <- Human_mRNAs[,c("Transcript_ID","cDNA_Length", "CDS_Start", "CDS_End")]

write.csv(Human_mRNAs, "Human_mRNAs.csv", row.names=F)
###############################################################################################################################################