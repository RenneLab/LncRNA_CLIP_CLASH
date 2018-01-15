#! /usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLASH_lncRNAs/Scripts")

library(splitstackshape)
library(dplyr)
library(data.table)

###########################################################################################################################################
### Change some fundamental settings to avoid R converting long integers into exponential notation
### This is an important step for this script to work as intended

digitsum <- function(x) sum(floor(x / 10^(0:(nchar(x) - 1))) %% 10)
options(scipen = 999)

###########################################################################################################################################
### Loop over each file to categorize them into specific seed types

INDIR = "../HybVienna/"
OUTDIR= "../SeedMatch/"

files <- list.files(path = INDIR ,pattern = "\\.csv$")

for (file in files) {
  Type <- list()
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  miR_Diagram <- temp$miR_Diagram
  miR_Diagram <- gsub("\\(", "1", miR_Diagram)
  miR_Diagram <- gsub("\\)", "1", miR_Diagram)
  miR_Diagram <- gsub("\\.", "0", miR_Diagram)
  for (i in 1:length(miR_Diagram)){
       if (nchar(miR_Diagram[i]) > 7) {
         p1 <- as.numeric(substr(miR_Diagram[i],2,2))
         p2 <- as.numeric(substr(miR_Diagram[i],3,3))
         p3 <- as.numeric(substr(miR_Diagram[i],4,4))
         p4 <- as.numeric(substr(miR_Diagram[i],5,5))
         p5 <- as.numeric(substr(miR_Diagram[i],6,6))
         p6 <- as.numeric(substr(miR_Diagram[i],7,7))
         p7 <- as.numeric(substr(miR_Diagram[i],8,8))
         if (sum(p1,p2,p3,p4,p5,p6,p7) == 7) {Type[i] <- "seed_7"
         } else if (sum(p1,p2,p3,p4,p5,p6) == 6) {Type[i] <- "seed_6"
         } else if (sum(p1,p2,p3,p4,p5,p6,p7) == 6) {Type[i] <- "seed_7m1"
         } else if (sum(p1,p2,p3,p4,p5,p6,p7) == 5) {Type[i] <- "seed_7m2"
         } else {Type[i] <- "other"}
       }
    else {Type[i] <- "seed_truncated"}
  }
  temp$Seed_Type <- unlist(Type)
  write.csv(temp,paste0(OUTDIR,file), row.names = F)
}

###########################################################################################################################################
### Loop over each file to categorize them into specific compensatory base-pairing types

files <- list.files(path = OUTDIR ,pattern = "\\.csv$")

for (file in files) {
  Type_3p <- list()
  temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
  miR_Diagram <- temp$miR_Diagram
  miR_Diagram <- gsub("\\(", "1", miR_Diagram)
  miR_Diagram <- gsub("\\)", "1", miR_Diagram)
  miR_Diagram <- gsub("\\.", "0", miR_Diagram)
  for (i in 1:length(miR_Diagram)){
    if (nchar(miR_Diagram[i]) > 11 && ! is.na(nchar(miR_Diagram[i]))) {
      Bind_3p  <- digitsum(as.numeric(substr(miR_Diagram[i],11,nchar(miR_Diagram[i]))))
      if (Bind_3p == 0) {Type_3p[i] <- "Absent"
      } else if (Bind_3p < 5) {Type_3p[i] <- "Weak"
      } else if (Bind_3p < 8) {Type_3p[i] <- "Moderate"
      } else if (Bind_3p > 7) {Type_3p[i] <- "Strong"
      }
    }
    else {Type_3p[i] <- "miR_truncated"}
  }
  temp$Bind_3p_Type <- unlist(Type_3p)
  write.csv(temp,paste0(OUTDIR,file), row.names = F)
}


###########################################################################################################################################
### Count seed types for all samples

files <- list.files(OUTDIR, pattern = "\\.csv$")

counting <- data.frame(matrix(NA, nrow=length(files), ncol=7))
colnames(counting) <- c("Sample","seed_7", "seed_6", "seed_7m1", "seed_7m2", "other", "seed_truncated")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
  counting[count,1] <- file
  counting[count,2] <- sum(temp$Seed_Type == "seed_7")
  counting[count,3] <- sum(temp$Seed_Type == "seed_6")
  counting[count,4] <- sum(temp$Seed_Type == "seed_7m1")
  counting[count,5] <- sum(temp$Seed_Type == "seed_7m2")
  counting[count,6] <- sum(temp$Seed_Type == "other")
  counting[count,7] <- sum(temp$Seed_Type == "seed_truncated")
  count= count+1
}
write.csv(counting, paste0(OUTDIR,"Counts/","SeedMatches.csv"), row.names=F)


# Running these counts only on complete diagrams

counting <- data.frame(matrix(NA, nrow=length(files), ncol=7))
colnames(counting) <- c("Sample","seed_7", "seed_6", "seed_7m1", "seed_7m2", "other", "seed_truncated")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
  temp <- temp[which(nchar(temp$Diagram) == temp$Read_end_3 - temp$Read_start_5 +1),]
  counting[count,1] <- file
  counting[count,2] <- sum(temp$Seed_Type == "seed_7")
  counting[count,3] <- sum(temp$Seed_Type == "seed_6")
  counting[count,4] <- sum(temp$Seed_Type == "seed_7m1")
  counting[count,5] <- sum(temp$Seed_Type == "seed_7m2")
  counting[count,6] <- sum(temp$Seed_Type == "other")
  counting[count,7] <- sum(temp$Seed_Type == "seed_truncated")
  count= count+1
  test <- temp
}
write.csv(counting, paste0(OUTDIR,"Counts/","SeedMatches_CompleteDiags.csv"), row.names=F)


###########################################################################################################################################
### Count seed pairing type for all samples

files <- list.files(OUTDIR, pattern = "\\.csv$")

counting <- data.frame(matrix(NA, nrow=length(files), ncol=6))
colnames(counting) <- c("Sample","Absent", "Weak", "Moderate", "Strong", "miR_truncated")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
  counting[count,1] <- file
  counting[count,2] <- sum(temp$Bind_3p_Type == "Absent")
  counting[count,3] <- sum(temp$Bind_3p_Type == "Weak")
  counting[count,4] <- sum(temp$Bind_3p_Type == "Moderate")
  counting[count,5] <- sum(temp$Bind_3p_Type == "Strong")
  counting[count,6] <- sum(temp$Bind_3p_Type == "miR_truncated")
  count= count+1
}
write.csv(counting, paste0(OUTDIR,"Counts/","Binding_3p.csv"), row.names=F)

# Running these counts only on complete diagrams

counting <- data.frame(matrix(NA, nrow=length(files), ncol=6))
colnames(counting) <- c("Sample","Absent", "Weak", "Moderate", "Strong", "miR_truncated")

count =1
for(file in files) {
  temp <- read.csv(paste0(OUTDIR,file), stringsAsFactors = F)
  temp <- temp[which(nchar(temp$Diagram) == temp$Read_end_3 - temp$Read_start_5 +1),]
  counting[count,1] <- file
  counting[count,2] <- sum(temp$Bind_3p_Type == "Absent")
  counting[count,3] <- sum(temp$Bind_3p_Type == "Weak")
  counting[count,4] <- sum(temp$Bind_3p_Type == "Moderate")
  counting[count,5] <- sum(temp$Bind_3p_Type == "Strong")
  counting[count,6] <- sum(temp$Bind_3p_Type == "miR_truncated")
  count= count+1
}
write.csv(counting, paste0(OUTDIR,"Counts/","Binding_3p_CompleteDiags.csv"), row.names=F)

###########################################################################################################################################