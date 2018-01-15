#! usr/bin/env RScript

rm(list=ls())

library(data.table)
library(splitstackshape)

#######################################################################################################################################
### Extracting information from GENCODE V19 for annotation
### Two files are written, one for exon annotation, the other for transcript annotation

reffile <- read.table("../PIPECLIP/ReferenceFiles/gencode.v19.annotation.bed", sep = "\t",stringsAsFactors = F, fill =T)
ref_transcript <- reffile[which(reffile$V8 == "transcript"),]
write.table(ref_transcript, "../PIPECLIP/ReferenceFiles/AllTranscripts.bed", quote = FALSE, row.names = FALSE, col.name =FALSE, sep = "\t")
ref_exon <- reffile[which(reffile$V8 == "exon"),]
write.table(ref_exon, "../PIPECLIP/ReferenceFiles/AllExons.bed", quote = FALSE, row.names = FALSE, col.name =FALSE, sep = "\t")

#######################################################################################################################################
### Extracting classes: exons

exons <- cSplit(ref_exon, 10, sep=";")
exons$V10_06 <- gsub("transcript_type ", "", exons$V10_06)
TranscriptTypes <- levels(factor(exons$V10_06))
codingRNA <- "protein_coding"
smallRNA <- c("Mt_rRNA", "miRNA", "Mt_tRNA", "rRNA", "snoRNA", "snRNA")
lncRNA <- TranscriptTypes[!TranscriptTypes %in% c(codingRNA,smallRNA)]
ref_exon$transcript_class[which(exons$V10_06 %in% lncRNA)] = "lncRNA"
ref_exon$transcript_class[which(exons$V10_06 %in% codingRNA)] = "mRNA"
ref_exon$transcript_class[which(exons$V10_06 %in% smallRNA)] = "smallRNA"
write.table(ref_exon, "../PIPECLIP/ReferenceFiles/AllExons.bed", quote = FALSE, row.names = FALSE, col.name =FALSE, sep = "\t")

#######################################################################################################################################
### Extracting classes: transcripts

transcripts <- cSplit(ref_transcript, 10, sep=";")
transcripts$V10_06 <- gsub("transcript_type ", "", transcripts$V10_06)
TranscriptTypes <- levels(factor(transcripts$V10_06))
codingRNA <- "protein_coding"
smallRNA <- c("Mt_rRNA", "miRNA", "Mt_tRNA", "rRNA", "snoRNA", "snRNA")
lncRNA <- TranscriptTypes[!TranscriptTypes %in% c(codingRNA,smallRNA)]
ref_transcript$transcript_class[which(transcripts$V10_06 %in% lncRNA)] = "lncRNA"
ref_transcript$transcript_class[which(transcripts$V10_06 %in% codingRNA)] = "mRNA"
ref_transcript$transcript_class[which(transcripts$V10_06 %in% smallRNA)] = "smallRNA"
write.table(ref_transcript, "../PIPECLIP/ReferenceFiles/AllTranscripts.bed", quote = FALSE, row.names = FALSE, col.name =FALSE, sep = "\t")

#######################################################################################################################################
### Pre-processing miRNA information obtained from miRBase for both KSHV and EBV miRNAs
#######################################################################################################################################
### KSHV miRNAs

KSHV_miR <- read.table("../miRNAs/KSHV_miR.fa", sep= "\t",stringsAsFactors = F)
KSHV_miR <- KSHV_miR[-47,]
KSHV_miR_seqs <- KSHV_miR[seq(2,50,2)]
KSHV_miR_names <- KSHV_miR[seq(1,49,2)]
KSHV_miR_seqs <- gsub("U","T", KSHV_miR_seqs)
KSHV_miR_names <- gsub(">","", KSHV_miR_names)
KSHV_miR_names <- gsub("Kaposi.*","", KSHV_miR_names)
KSHV_miR <- data.frame(cbind(KSHV_miR_names, KSHV_miR_seqs))
KSHV_miR <- cSplit(KSHV_miR,1," ")
colnames(KSHV_miR) <- c("miRNA_mature","miRNA","miRBase_ID")
KSHV_miR$RevComp <- KSHV_miR$miRNA_mature
KSHV_miR$RevComp <- gsub("A","t", KSHV_miR$RevComp)
KSHV_miR$RevComp <- gsub("T","a", KSHV_miR$RevComp)
KSHV_miR$RevComp <- gsub("G","c", KSHV_miR$RevComp)
KSHV_miR$RevComp <- gsub("C","g", KSHV_miR$RevComp)
for (i in 1:length(KSHV_miR$Revcomp)) {
KSHV_miR$RevComp[i] <- paste(rev(unlist(strsplit(KSHV_miR$RevComp[i], split=""))), collapse="")
}
KSHV_miR$RevComp <- toupper(KSHV_miR$RevComp)
KSHV_miR$SeedMatch_7 <- substr(KSHV_miR$RevComp, start=2, stop=8)
KSHV_miR$SeedMatch_6 <- substr(KSHV_miR$RevComp, start=2, stop=7)
write.csv(KSHV_miR, "../miRNAs/KSHV_miRNAs.csv", row.names = F)

#######################################################################################################################################
### EBV miRNAs

EBV_miR <- read.table("../miRNAs/EBV_miR.fa", sep= "\t",stringsAsFactors = F)
EBV_miR <- EBV_miR[-c(17,64,85),]
EBV_miR_seqs <- EBV_miR[seq(2,88,2)]
EBV_miR_names <- EBV_miR[seq(1,87,2)]
EBV_miR_seqs <- gsub("U","T", EBV_miR_seqs)
EBV_miR_names <- gsub(">","", EBV_miR_names)
EBV_miR_names <- gsub("Epstein.*","", EBV_miR_names)
EBV_miR <- data.frame(cbind(EBV_miR_names, EBV_miR_seqs))
EBV_miR <- cSplit(EBV_miR,1," ")
colnames(EBV_miR) <- c("miRNA_mature","miRNA","miRBase_ID")
EBV_miR$RevComp <- EBV_miR$miRNA_mature
EBV_miR$RevComp <- gsub("A","t", EBV_miR$RevComp)
EBV_miR$RevComp <- gsub("T","a", EBV_miR$RevComp)
EBV_miR$RevComp <- gsub("G","c", EBV_miR$RevComp)
EBV_miR$RevComp <- gsub("C","g", EBV_miR$RevComp)
for (i in 1:length(EBV_miR$Revcomp)) {
  EBV_miR$RevComp[i] <- paste(rev(unlist(strsplit(EBV_miR$RevComp[i], split=""))), collapse="")
}
EBV_miR$RevComp <- toupper(EBV_miR$RevComp)
EBV_miR$SeedMatch_7 <- substr(EBV_miR$RevComp, start=2, stop=8)
EBV_miR$SeedMatch_6 <- substr(EBV_miR$RevComp, start=2, stop=7)
write.csv(EBV_miR, "../miRNAs/EBV_miRNAs.csv", row.names = F)

#######################################################################################################################################
