#! usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/Scripts")

library(data.table)
library(splitstackshape)
library(stringr)
library(gplots)

#######################################################################################################################################
### Specify input filenames

INDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/Targets/"
OUTDIR <- "/ufrc/renne/sunantha.s/research/CLIP_lncRNAs/MetaAnalysis/PIPECLIP/Pooled/"

EH_files_6 <- c("EH_Jijoye_1_EBV_6mer.csv", "EH_Jijoye_2_EBV_6mer.csv", "EH_Jijoye_3_EBV_6mer.csv","EH_Jijoye_4_EBV_6mer.csv",
                "EH_Jijoye_5_EBV_6mer.csv", "EH_Jijoye_6_EBV_6mer.csv")
EP_files_6 <- c("EP_EF3DAGO2_1_B95_8_6mer.csv", "EP_LCL35_1_B95_8_6mer.csv", "EP_LCLBAC_1_B95_8_6mer.csv", "EP_LCLBACD1_1_B95_8_6mer.csv",
                "EP_LCLBACD3_1_B95_8_6mer.csv", "KP_BC1_1_EBV_6mer.csv", "KP_BC1_2_EBV_6mer.csv")
KH_files_6 <- c("KH_BC3_1_BC3_6mer.csv", "KH_BC3_2_BC3_6mer.csv", "KH_BC3_3_BC3_6mer.csv", "KH_BCBL1_1_KSHV_6mer.csv", "KH_BCBL1_2_KSHV_6mer.csv",
                "KH_BCBL1_3_KSHV_6mer.csv", "KH_BCBL1_4_KSHV_6mer.csv")
KP_files_6 <- c("KP_BC1_1_KSHV_6mer.csv","KP_BC1_2_KSHV_6mer.csv", "KP_BC3_1_BC3_6mer.csv", "KP_BC3_2_BC3_6mer.csv")

EH_files_7 <- c("EH_Jijoye_1_EBV_7mer.csv", "EH_Jijoye_2_EBV_7mer.csv", "EH_Jijoye_3_EBV_7mer.csv","EH_Jijoye_4_EBV_7mer.csv",
                "EH_Jijoye_5_EBV_7mer.csv", "EH_Jijoye_6_EBV_7mer.csv")
EP_files_7 <- c("EP_EF3DAGO2_1_B95_8_7mer.csv", "EP_LCL35_1_B95_8_7mer.csv", "EP_LCLBAC_1_B95_8_7mer.csv", "EP_LCLBACD1_1_B95_8_7mer.csv",
                "EP_LCLBACD3_1_B95_8_7mer.csv", "KP_BC1_1_EBV_7mer.csv", "KP_BC1_2_EBV_7mer.csv")
KH_files_7 <- c("KH_BC3_1_BC3_7mer.csv", "KH_BC3_2_BC3_7mer.csv", "KH_BC3_3_BC3_7mer.csv", "KH_BCBL1_1_KSHV_7mer.csv", "KH_BCBL1_2_KSHV_7mer.csv",
                "KH_BCBL1_3_KSHV_7mer.csv", "KH_BCBL1_4_KSHV_7mer.csv")
KP_files_7 <- c("KP_BC1_1_KSHV_7mer.csv","KP_BC1_2_KSHV_7mer.csv", "KP_BC3_1_BC3_7mer.csv", "KP_BC3_2_BC3_7mer.csv")

groups <- list(EH_files_6, EP_files_6, KH_files_6, KP_files_6, EH_files_7, EP_files_7, KH_files_7, KP_files_7)
groupnames <- list("EH_files_6", "EP_files_6", "KH_files_6", "KP_files_6", "EH_files_7", "EP_files_7", "KH_files_7", "KP_files_7")


#######################################################################################################################################
### For every group, pool data from corresponding input files

for(g in 1:length(groups)) {
  group <- groups[g]
  groupname <- groupnames[g]
  datalist <- list()
  i <- 1
for(file in unlist(group)) {
  temp <- read.csv(paste0(INDIR, file), stringsAsFactors = F)
  temp$Sample <- file
  datalist[[i]] <- temp
  i <- i+1
}
  df <- do.call(rbind,datalist)
  write.csv(df, paste0(OUTDIR,groupname,".csv"), row.names=F)
}

#######################################################################################################################################
### For 6mer files, filter unique genes and create a Venn diagram of overlap between the four categories

EH_6mer <- read.csv(paste0(OUTDIR,"EH_files_6.csv"), stringsAsFactors = F)
EP_6mer <- read.csv(paste0(OUTDIR,"EP_files_6.csv"), stringsAsFactors = F)
KH_6mer <- read.csv(paste0(OUTDIR,"KH_files_6.csv"), stringsAsFactors = F)
KP_6mer <- read.csv(paste0(OUTDIR,"KP_files_6.csv"), stringsAsFactors = F)

EH_6mer_m <- unique(EH_6mer[which(EH_6mer$transcript_class=="mRNA"),"gene_name"])
EP_6mer_m <- unique(EP_6mer[which(EP_6mer$transcript_class=="mRNA"),"gene_name"])
KH_6mer_m <- unique(KH_6mer[which(KH_6mer$transcript_class=="mRNA"),"gene_name"])
KP_6mer_m <- unique(KP_6mer[which(KP_6mer$transcript_class=="mRNA"),"gene_name"])

l <- list(EH_6mer_m, EP_6mer_m, KH_6mer_m, KP_6mer_m)
names(l) <- c("EH_6mer", "EP_6mer", "KH_6mer", "KP_6mer")
venn(list(EH_6mer_m, EP_6mer_m, KH_6mer_m, KP_6mer_m))

mRNA_table <- crossprod(table(stack(l)))
write.csv(mRNA_table, paste0(OUTDIR,"mRNA_6mer.csv"))

EH_6mer_l <- unique(EH_6mer[which(EH_6mer$transcript_class=="lncRNA"),"gene_name"])
EP_6mer_l <- unique(EP_6mer[which(EP_6mer$transcript_class=="lncRNA"),"gene_name"])
KH_6mer_l <- unique(KH_6mer[which(KH_6mer$transcript_class=="lncRNA"),"gene_name"])
KP_6mer_l <- unique(KP_6mer[which(KP_6mer$transcript_class=="lncRNA"),"gene_name"])

l <- list(EH_6mer_l, EP_6mer_l, KH_6mer_l, KP_6mer_l)
names(l) <- c("EH_6mer", "EP_6mer", "KH_6mer", "KP_6mer")
venn(list(EH_6mer_l, EP_6mer_l, KH_6mer_l, KP_6mer_l))

lncRNA_table <- crossprod(table(stack(l)))
write.csv(lncRNA_table, paste0(OUTDIR,"lncRNA_6mer.csv"))

#######################################################################################################################################
### For 7mer files, filter unique genes and create a Venn diagram of overlap between the four categories

EH_7mer <- read.csv(paste0(OUTDIR,"EH_files_7.csv"), stringsAsFactors = F)
EP_7mer <- read.csv(paste0(OUTDIR,"EP_files_7.csv"), stringsAsFactors = F)
KH_7mer <- read.csv(paste0(OUTDIR,"KH_files_7.csv"), stringsAsFactors = F)
KP_7mer <- read.csv(paste0(OUTDIR,"KP_files_7.csv"), stringsAsFactors = F)

EH_7mer_m <- unique(EH_7mer[which(EH_7mer$transcript_class=="mRNA"),"gene_name"])
EP_7mer_m <- unique(EP_7mer[which(EP_7mer$transcript_class=="mRNA"),"gene_name"])
KH_7mer_m <- unique(KH_7mer[which(KH_7mer$transcript_class=="mRNA"),"gene_name"])
KP_7mer_m <- unique(KP_7mer[which(KP_7mer$transcript_class=="mRNA"),"gene_name"])

l <- list(EH_7mer_m, EP_7mer_m, KH_7mer_m, KP_7mer_m)
names(l) <- c("EH_7mer", "EP_7mer", "KH_7mer", "KP_7mer")
venn(list(EH_7mer_m, EP_7mer_m, KH_7mer_m, KP_7mer_m))

mRNA_table <- crossprod(table(stack(l)))
write.csv(mRNA_table, paste0(OUTDIR,"mRNA_7mer.csv"))

EH_7mer_l <- unique(EH_7mer[which(EH_7mer$transcript_class=="lncRNA"),"gene_name"])
EP_7mer_l <- unique(EP_7mer[which(EP_7mer$transcript_class=="lncRNA"),"gene_name"])
KH_7mer_l <- unique(KH_7mer[which(KH_7mer$transcript_class=="lncRNA"),"gene_name"])
KP_7mer_l <- unique(KP_7mer[which(KP_7mer$transcript_class=="lncRNA"),"gene_name"])

l <- list(EH_7mer_l, EP_7mer_l, KH_7mer_l, KP_7mer_l)
names(l) <- c("EH_7mer", "EP_7mer", "KH_7mer", "KP_7mer")
venn(list(EH_7mer_l, EP_7mer_l, KH_7mer_l, KP_7mer_l))

lncRNA_table <- crossprod(table(stack(l)))
write.csv(lncRNA_table, paste0(OUTDIR,"lncRNA_7mer.csv"))

#######################################################################################################################################