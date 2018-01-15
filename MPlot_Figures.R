#! /usr/bin/env Rscript

rm(list=ls())

setwd("/ufrc/renne/sunantha.s/research/CLASH_lncRNAs/Scripts")

library(splitstackshape)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(devEMF)

###########################################################################################################################################
### For every sample create a directory with subdirectories for every biorep and four files for every 
### class (5'miRNA/ 3'miRNA and cellular/viral miRNA)

INDIR ="../HybVienna/"
OUTDIR = "../Figures/"

files <- list.files(path = INDIR ,pattern = "\\.csv$")
sample_names <- unique(sub("_hyb_vienna_.*.csv", "", files))
classes <- c("5mi3m", "5v3m", "5m3mi", "5m3v")

for (samp in sample_names) {
  dir.create(paste0(OUTDIR,samp), showWarnings = F)
  for(class in classes) {
    dir.create(paste0(OUTDIR,samp,"/",class), showWarnings = F)
  }
}

###########################################################################################################################################
### Define a function to calculate percentage binding at each position of the miRNA and plot the profile
### Plotting is done for individual miRNAs and also a cumulative plot is generated averaging over all miRNAs

FoldingPattern <- function(samp,class) {
  
  if (file.exists(paste0(INDIR,samp,"_hyb_vienna_",class,".csv"))) {
    
    hyb_vienna <- read.csv(paste0(INDIR,samp,"_hyb_vienna_",class,".csv"))
    hyb_vienna$miR_Diagram <- gsub("\\(", "1", hyb_vienna$miR_Diagram)
    hyb_vienna$miR_Diagram <- gsub("\\)", "1", hyb_vienna$miR_Diagram)
    hyb_vienna$miR_Diagram <- gsub("\\.", "0", hyb_vienna$miR_Diagram)
    
    if (class == "5mi3m" | class == "5v3m") {
      unique_miR <- as.vector(unique(hyb_vienna$GeneName_5))}
    if (class == "5m3mi" | class == "5m3v") {
      unique_miR <- as.vector(unique(hyb_vienna$GeneName_3))}
    
    #Looping over each miRNA and creating plots by miRNA
    
    Tpos <- rep.int(0,24)
    Tcount <- rep.int(0,24)
    Tpercent <- rep.int(0,24)
    
    
    for (miR in unique_miR) {
      
      if (class == "5mi3m" | class == "5v3m") {
        miR_hyb <- hyb_vienna[hyb_vienna$GeneName_5 == miR,]}
      
      if (class == "5m3mi" | class == "5m3v") {
        miR_hyb <- hyb_vienna[hyb_vienna$GeneName_3 == miR,]}
      
      max_miR_length <- max(nchar(as.character(miR_hyb$miR_Diagram)))
      
      if (max_miR_length != 0) {
        
        pos <- rep.int(0,max_miR_length)
        count <- rep.int(0,max_miR_length)
        percent <- rep.int(0,max_miR_length)
        
        for(j in 1:length(miR_hyb$miR_Diagram)) {
          for(i in 1:max_miR_length) {
            if(substr(miR_hyb$miR_Diagram[j],i,i) == 1) {
              pos[i] <- pos[i]+1
              count[i] <- count[i]+1
              Tpos[i] <- Tpos[i]+1
              Tcount[i] <- Tcount[i] +1
            } else if(substr(miR_hyb$miR_Diagram[j],i,i) == 0) {
              count[i] <- count[i]+1
              Tcount[i] <- Tcount[i] +1
            }
          }
        }
        
        for(i in 1:max_miR_length) {
          percent[i] <- 100*(pos[i]/count[i])
        }
        names <- rep("",max_miR_length)
        for(i in 1:max_miR_length) {
          names[i] <- paste0("pos",i)
        }
        names(percent) <- names
        
        #Plotting
        
        svg(paste0(OUTDIR,samp,"/",class,"/",miR,"_",class,".svg"))
        par(las=2)
        plot(percent,xaxt="n", type="p", pch=19, ylab = "% targets", xlab= "miRNA nucleotide position", col="orangered", ylim=c(0,100))
        lines(percent,type="c", lty=1, lwd=2, col="blue")
        axis(1,at=1:max_miR_length,labels=1:max_miR_length)
        dev.off()
        
      }
    }
    
    # Averaging over all miRNAs and plotting
    
    for(i in 1:24) {
      if(Tcount[i] !=0) {
        Tpercent[i] <- 100*(Tpos[i]/Tcount[i])
      }
      Tnames <- rep("",length(Tpercent))
      for(i in 1:length(Tpercent)) {
        Tnames[i] <- paste0("pos",i)
      }
      names(Tpercent) <- Tnames
      
      #Plotting
      
      svg(paste0(OUTDIR,samp,"/",class,"/",samp,"_",class,".svg"))
      par(las=2)
      plot(Tpercent,xaxt="n", type="p", pch=19, ylab = "% targets", xlab= "miRNA nucleotide position",col="orangered", ylim=c(0,100))
      lines(Tpercent,type="c", lty=1, lwd=2, col="blue")
      axis(1,at=1:length(Tpercent),labels=1:length(Tpercent))
      dev.off()
    }
  }
}

###########################################################################################################################################
### Run the function over all samples and all classes 

for (samp in sample_names) {
  for(class in classes) {
    FoldingPattern(samp,class)
  }
}
###########################################################################################################################################
