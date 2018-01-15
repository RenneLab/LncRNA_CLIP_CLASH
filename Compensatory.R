#! /usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLASH_lncRNAs/Scripts")

###############################################################################################################################################
# Obtain pair-wise counts of seed type and compensatory base pairing from all files

INDIR <- "../SeedMatch/"
OUTDIR <- "../SeedMatch/Counts/"

files_5mi3m <- list.files(path = INDIR ,pattern = "5mi3m\\.csv$")
files_5m3mi <- list.files(path = INDIR ,pattern = "5m3mi\\.csv$")
files_5v3m <- list.files(path = INDIR ,pattern = "5v3m\\.csv$")
files_5m3v <- list.files(path = INDIR ,pattern = "5m3v\\.csv$")

list_5mi3m <- list()
count <- 1
for(file in files_5mi3m) {
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  temp <- temp[,c(26,27)]
  list_5mi3m[[count]] <- temp
  count <- count+1
}
Comp_5mi3m <- do.call(rbind, list_5mi3m)
Comp_5mi3m <- as.data.frame(table(Comp_5mi3m))
write.csv(Comp_5mi3m, paste0(OUTDIR,"Comp_5mi3m.csv"), row.names = F)

list_5m3mi <- list()
count <- 1
for(file in files_5m3mi) {
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  temp <- temp[,c(26,27)]
  list_5m3mi[[count]] <- temp
  count <- count+1
}
Comp_5m3mi <- do.call(rbind, list_5m3mi)
Comp_5m3mi <- as.data.frame(table(Comp_5m3mi))
write.csv(Comp_5m3mi, paste0(OUTDIR,"Comp_5m3mi.csv"), row.names = F)

list_5v3m <- list()
count <- 1
for(file in files_5v3m) {
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  temp <- temp[,c(26,27)]
  list_5v3m[[count]] <- temp
  count <- count+1
}
Comp_5v3m <- do.call(rbind, list_5v3m)
Comp_5v3m <- as.data.frame(table(Comp_5v3m))
write.csv(Comp_5v3m, paste0(OUTDIR,"Comp_5v3m.csv"), row.names = F)

list_5m3v <- list()
count <- 1
for(file in files_5m3v) {
  temp <- read.csv(paste0(INDIR,file), stringsAsFactors = F)
  temp <- temp[,c(26,27)]
  list_5m3v[[count]] <- temp
  count <- count+1
}
Comp_5m3v <- do.call(rbind, list_5m3v)
Comp_5m3v <- as.data.frame(table(Comp_5m3v))
write.csv(Comp_5m3v, paste0(OUTDIR,"Comp_5m3v.csv"), row.names = F)
###############################################################################################################################################