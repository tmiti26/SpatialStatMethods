rm(list=ls())
library(sp)
#library(rmapshaper)
library(ggplot2)
library(magrittr)
#library(sf)
#library(spatstat)
library(goftest)
library(readr)
library(kSamples)
library(plotrix)
library(tidyverse)
#library(pBrackets)

#setup and import data
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))

numbTissue <- 99
N <- 12
numberList <- seq(from = 1, to = numbTissue, by = 1)

## ---- data-driven x-max: the distance where the nearest-stroma distance frequency
## distributions of all conditions have merged and dropped to ~0 (< thr of the peak).
## Reads one variation (i=1) of myNNP per condition purely to set the axis.
condListX <- c("3D", "R", "Exp30", "Exp60", "Exp90", "Exp120", "Exp180", "Exp320")
myMergeXmax <- function(condList = condListX, thr = 0.02, pad = 3, gmax = 200){
  gx <- seq(0, gmax, length.out = 400); maxDens <- rep(0, length(gx))
  for (cc in condList){
    fn <- paste0(1, "_myNNP", cc, "_", N, ".csv")
    if (!file.exists(fn)) next
    d <- read.csv(fn)[, 1]; d <- d[is.finite(d)]
    if (length(d) < 3) next
    maxDens <- pmax(maxDens, density(d, from = 0, to = gmax, n = length(gx))$y)
  }
  peak <- max(maxDens); if (!is.finite(peak) || peak <= 0) return(gmax)
  idx <- which(maxDens > thr * peak)
  if (length(idx) == 0) return(gmax)
  min(gx[max(idx)] + pad, gmax)
}


for(i in numberList){
  number <- i
  
  #numb <- as.numeric(paste(readLines("polNumb.txt"), collapse=" "))
  linn <-readLines(file("Rparams.txt",open="r"))
  close(file("Rparams.txt",open="r"))
  quadratSize <- as.numeric(linn[1]) #1000
  PCF_r <- as.numeric(linn[2]) # 450
  J_r <- as.numeric(linn[3]) #250
  Z_r_length <- as.numeric(linn[4]) #500
  dens_r_length <- as.numeric(linn[5]) #500
  PCF_r_length <- as.integer(PCF_r*2)
  J_r_length <- as.integer(J_r*2)
  
  ############################################## Sample Data Analysis Quadrats Start ######################################
  #########################################################################################################################
  #reading from files the average values
  
  #nam <- paste("averObsP", number, sep = "")
  assign(paste("RDFP", number, sep = ""),  (read.csv(file = paste(i,"_myenvRDFP3D_", N, ".csv",sep=""))))
  #assign(paste("RDFN", number, sep = ""),  (read.csv(file = paste(i,"_myenvRDFN3D_", N, ".csv",sep=""))))
  
  assign(paste("RDFPRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPR_", N, ".csv",sep=""))))
  #assign(paste("RDFNRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNR_", N, ".csv",sep=""))))
  
  assign(paste("RDFPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp30_", N, ".csv",sep=""))))
  #assign(paste("RDFNExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp30_", N, ".csv",sep=""))))
  
  assign(paste("RDFPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp60_", N, ".csv",sep=""))))
  #assign(paste("RDFNExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp60_", N, ".csv",sep=""))))
  
  assign(paste("RDFPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp90_", N, ".csv",sep=""))))
  #assign(paste("RDFNExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp90_", N, ".csv",sep=""))))
  
  #assign(paste("RDFPExp120_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp120_", N, ".csv",sep=""))))
  #assign(paste("RDFNExp120_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp120_", N, ".csv",sep=""))))
  
  # assign(paste("RDFPExp150_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp150_", N, ".csv",sep=""))))
  # assign(paste("RDFNExp150_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp150_", N, ".csv",sep=""))))
  
  assign(paste("RDFPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp180_", N, ".csv",sep=""))))
  #assign(paste("RDFNExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp180_", N, ".csv",sep=""))))
  
  # assign(paste("RDFPExp240_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp240_", N, ".csv",sep=""))))
  # assign(paste("RDFNExp240_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp240_", N, ".csv",sep=""))))
  
  assign(paste("RDFPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp320_", N, ".csv",sep=""))))
  #assign(paste("RDFNExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFNExp320_", N, ".csv",sep=""))))
  
}

#********************************* RDF ***********************************
x <- RDFP1$r

for (i in 1:numbTissue){
  # Create new variable names
  pop.name = paste0("RDFP",i)
  #pop.nameN =  paste0("RDFN",i)
  pop.nameR = paste0("RDFPRand",i)
  #pop.nameRN = paste0("RDFNRand",i)
  pop.namePExp30 = paste0("RDFPExp30_",i)
  #pop.nameNExp30 =  paste0("RDFNExp30_",i)
  pop.namePExp60 = paste0("RDFPExp60_",i)
  #pop.nameNExp60 =  paste0("RDFNExp60_",i)
  pop.namePExp90 = paste0("RDFPExp90_",i)
  #pop.nameNExp90 =  paste0("RDFNExp90_",i)
  #pop.namePExp120 = paste0("RDFPExp120_",i)
  #pop.nameNExp120 =  paste0("RDFNExp120_",i)
  # pop.namePExp150 = paste0("RDFPExp150_",i)
  # pop.nameNExp150 =  paste0("RDFNExp150_",i)
  pop.namePExp180 = paste0("RDFPExp180_",i)
  #pop.nameNExp180 =  paste0("RDFNExp180_",i)
  # pop.namePExp240 = paste0("RDFPExp240_",i)
  # pop.nameNExp240 =  paste0("RDFNExp240_",i)
  pop.namePExp320 = paste0("RDFPExp320_",i)
  #pop.nameNExp320 =  paste0("RDFNExp320_",i)
  
  assign(paste0("diffRand", i),((eval(parse(text = pop.nameR))$obs)/(eval(parse(text = pop.nameR))$mmean)))
  #assign(paste0("diffNRand", i),((eval(parse(text = pop.nameRN))$obs)/(eval(parse(text = pop.nameRN))$mmean)))
  assign(paste0("diffRL", i),((eval(parse(text = pop.nameR))$lo)/(eval(parse(text = pop.nameR))$mmean)))
  assign(paste0("diffRH", i),((eval(parse(text = pop.nameR))$hi)/(eval(parse(text = pop.nameR))$mmean)))
  #assign(paste0("diffNRL", i),((eval(parse(text = pop.nameRN))$lo)/(eval(parse(text = pop.nameRN))$mmean)))
  #assign(paste0("diffNRH", i),((eval(parse(text = pop.nameRN))$hi)/(eval(parse(text = pop.nameRN))$mmean)))
  
  assign(paste0("diff", i),((eval(parse(text = pop.name))$obs)/(eval(parse(text = pop.name))$mmean)))
  #assign(paste0("diffN", i),((eval(parse(text = pop.nameN))$obs)/(eval(parse(text = pop.nameN))$mmean)))
  assign(paste0("diffL", i),((eval(parse(text = pop.name))$lo)/(eval(parse(text = pop.name))$mmean)))
  assign(paste0("diffH", i),((eval(parse(text = pop.name))$hi)/(eval(parse(text = pop.name))$mmean)))
  #assign(paste0("diffNL", i),((eval(parse(text = pop.nameN))$lo)/(eval(parse(text = pop.nameN))$mmean)))
  #assign(paste0("diffNH", i),((eval(parse(text = pop.nameN))$hi)/(eval(parse(text = pop.nameN))$mmean)))
  
  assign(paste0("diffExp30_", i),((eval(parse(text = pop.namePExp30))$obs)/(eval(parse(text = pop.namePExp30))$mmean)))
  #assign(paste0("diffNExp30_", i),((eval(parse(text = pop.nameNExp30))$obs)/(eval(parse(text = pop.nameNExp30))$mmean)))
  assign(paste0("diffLExp30_", i),((eval(parse(text = pop.namePExp30))$lo)/(eval(parse(text = pop.namePExp30))$mmean)))
  assign(paste0("diffHExp30_", i),((eval(parse(text = pop.namePExp30))$hi)/(eval(parse(text = pop.namePExp30))$mmean)))
  #assign(paste0("diffNLExp30_", i),((eval(parse(text = pop.nameNExp30))$lo)/(eval(parse(text = pop.nameNExp30))$mmean)))
  #assign(paste0("diffNHExp30_", i),((eval(parse(text = pop.nameNExp30))$hi)/(eval(parse(text = pop.nameNExp30))$mmean)))
  
  assign(paste0("diffExp60_", i),((eval(parse(text = pop.namePExp60))$obs)/(eval(parse(text = pop.namePExp60))$mmean)))
  #assign(paste0("diffNExp60_", i),((eval(parse(text = pop.nameNExp60))$obs)/(eval(parse(text = pop.nameNExp60))$mmean)))
  assign(paste0("diffLExp60_", i),((eval(parse(text = pop.namePExp60))$lo)/(eval(parse(text = pop.namePExp60))$mmean)))
  assign(paste0("diffHExp60_", i),((eval(parse(text = pop.namePExp60))$hi)/(eval(parse(text = pop.namePExp60))$mmean)))
  #assign(paste0("diffNLExp60_", i),((eval(parse(text = pop.nameNExp60))$lo)/(eval(parse(text = pop.nameNExp60))$mmean)))
  #assign(paste0("diffNHExp60_", i),((eval(parse(text = pop.nameNExp60))$hi)/(eval(parse(text = pop.nameNExp60))$mmean)))
  
  assign(paste0("diffExp90_", i),((eval(parse(text = pop.namePExp90))$obs)/(eval(parse(text = pop.namePExp90))$mmean)))
  #assign(paste0("diffNExp90_", i),((eval(parse(text = pop.nameNExp90))$obs)/(eval(parse(text = pop.nameNExp90))$mmean)))
  assign(paste0("diffLExp90_", i),((eval(parse(text = pop.namePExp90))$lo)/(eval(parse(text = pop.namePExp90))$mmean)))
  assign(paste0("diffHExp90_", i),((eval(parse(text = pop.namePExp90))$hi)/(eval(parse(text = pop.namePExp90))$mmean)))
  #assign(paste0("diffNLExp90_", i),((eval(parse(text = pop.nameNExp90))$lo)/(eval(parse(text = pop.nameNExp90))$mmean)))
  #assign(paste0("diffNHExp90_", i),((eval(parse(text = pop.nameNExp90))$hi)/(eval(parse(text = pop.nameNExp90))$mmean)))
  
  # assign(paste0("diffExp150_", i),((eval(parse(text = pop.namePExp150))$obs)/(eval(parse(text = pop.namePExp150))$mmean)))
  # assign(paste0("diffNExp150_", i),((eval(parse(text = pop.nameNExp150))$obs)/(eval(parse(text = pop.nameNExp150))$mmean)))
  # assign(paste0("diffLExp150_", i),((eval(parse(text = pop.namePExp150))$lo)/(eval(parse(text = pop.namePExp150))$mmean)))
  # assign(paste0("diffHExp150_", i),((eval(parse(text = pop.namePExp150))$hi)/(eval(parse(text = pop.namePExp150))$mmean)))
  # assign(paste0("diffNLExp150_", i),((eval(parse(text = pop.nameNExp150))$lo)/(eval(parse(text = pop.nameNExp150))$mmean)))
  # assign(paste0("diffNHExp150_", i),((eval(parse(text = pop.nameNExp150))$hi)/(eval(parse(text = pop.nameNExp150))$mmean)))
  
  assign(paste0("diffExp180_", i),((eval(parse(text = pop.namePExp180))$obs)/(eval(parse(text = pop.namePExp180))$mmean)))
  #assign(paste0("diffNExp180_", i),((eval(parse(text = pop.nameNExp180))$obs)/(eval(parse(text = pop.nameNExp180))$mmean)))
  assign(paste0("diffLExp180_", i),((eval(parse(text = pop.namePExp180))$lo)/(eval(parse(text = pop.namePExp180))$mmean)))
  assign(paste0("diffHExp180_", i),((eval(parse(text = pop.namePExp180))$hi)/(eval(parse(text = pop.namePExp180))$mmean)))
  #assign(paste0("diffNLExp180_", i),((eval(parse(text = pop.nameNExp180))$lo)/(eval(parse(text = pop.nameNExp180))$mmean)))
  #assign(paste0("diffNHExp180_", i),((eval(parse(text = pop.nameNExp180))$hi)/(eval(parse(text = pop.nameNExp180))$mmean)))
  
  # assign(paste0("diffExp240_", i),((eval(parse(text = pop.namePExp240))$obs)/(eval(parse(text = pop.namePExp240))$mmean)))
  # assign(paste0("diffNExp240_", i),((eval(parse(text = pop.nameNExp240))$obs)/(eval(parse(text = pop.nameNExp240))$mmean)))
  # assign(paste0("diffLExp240_", i),((eval(parse(text = pop.namePExp240))$lo)/(eval(parse(text = pop.namePExp240))$mmean)))
  # assign(paste0("diffHExp240_", i),((eval(parse(text = pop.namePExp240))$hi)/(eval(parse(text = pop.namePExp240))$mmean)))
  # assign(paste0("diffNLExp240_", i),((eval(parse(text = pop.nameNExp240))$lo)/(eval(parse(text = pop.nameNExp240))$mmean)))
  # assign(paste0("diffNHExp240_", i),((eval(parse(text = pop.nameNExp240))$hi)/(eval(parse(text = pop.nameNExp240))$mmean)))
  
  assign(paste0("diffExp320_", i),((eval(parse(text = pop.namePExp320))$obs)/(eval(parse(text = pop.namePExp320))$mmean)))
  #assign(paste0("diffNExp320_", i),((eval(parse(text = pop.nameNExp320))$obs)/(eval(parse(text = pop.nameNExp320))$mmean)))
  assign(paste0("diffLExp320_", i),((eval(parse(text = pop.namePExp320))$lo)/(eval(parse(text = pop.namePExp320))$mmean)))
  assign(paste0("diffHExp320_", i),((eval(parse(text = pop.namePExp320))$hi)/(eval(parse(text = pop.namePExp320))$mmean)))
  #assign(paste0("diffNLExp320_", i),((eval(parse(text = pop.nameNExp320))$lo)/(eval(parse(text = pop.nameNExp320))$mmean)))
  #assign(paste0("diffNHExp320_", i),((eval(parse(text = pop.nameNExp320))$hi)/(eval(parse(text = pop.nameNExp320))$mmean)))
  
}

RDFPRandtot <- rep(0, length(diff1))
RDFNRandtot <- rep(0, length(diff1))
RDFPRandtotL <- rep(0, length(diff1))
RDFPRandtotH <- rep(0, length(diff1))
RDFNRandtotH <- rep(0, length(diff1))
RDFNRandtotL <- rep(0, length(diff1))
RDFPtot <- rep(0, length(diff1))
RDFNtot <- rep(0, length(diff1))
RDFPtotL <- rep(0, length(diff1))
RDFPtotH <- rep(0, length(diff1))
RDFNtotH <- rep(0, length(diff1))
RDFNtotL <- rep(0, length(diff1))
RDFPExp30 <- rep(0, length(diff1))
RDFNExp30 <- rep(0, length(diff1))
RDFPExp30L <- rep(0, length(diff1))
RDFPExp30H <- rep(0, length(diff1))
RDFNExp30H <- rep(0, length(diff1))
RDFNExp30L <- rep(0, length(diff1))
RDFPExp60 <- rep(0, length(diff1))
RDFNExp60 <- rep(0, length(diff1))
RDFPExp60L <- rep(0, length(diff1))
RDFPExp60H <- rep(0, length(diff1))
RDFNExp60H <- rep(0, length(diff1))
RDFNExp60L <- rep(0, length(diff1))
RDFPExp90 <- rep(0, length(diff1))
RDFNExp90 <- rep(0, length(diff1))
RDFPExp90L <- rep(0, length(diff1))
RDFPExp90H <- rep(0, length(diff1))
RDFNExp90H <- rep(0, length(diff1))
RDFNExp90L <- rep(0, length(diff1))
RDFPExp120 <- rep(0, length(diff1))
RDFNExp120 <- rep(0, length(diff1))
RDFPExp120L <- rep(0, length(diff1))
RDFPExp120H <- rep(0, length(diff1))
RDFNExp120H <- rep(0, length(diff1))
RDFNExp120L <- rep(0, length(diff1))
# RDFPExp150 <- rep(0, length(diff1))
# RDFNExp150 <- rep(0, length(diff1))
# RDFPExp150L <- rep(0, length(diff1))
# RDFPExp150H <- rep(0, length(diff1))
# RDFNExp150H <- rep(0, length(diff1))
# RDFNExp150L <- rep(0, length(diff1))
RDFPExp180 <- rep(0, length(diff1))
RDFNExp180 <- rep(0, length(diff1))
RDFPExp180L <- rep(0, length(diff1))
RDFPExp180H <- rep(0, length(diff1))
RDFNExp180H <- rep(0, length(diff1))
RDFNExp180L <- rep(0, length(diff1))
# RDFPExp240 <- rep(0, length(diff1))
# RDFNExp240 <- rep(0, length(diff1))
# RDFPExp240L <- rep(0, length(diff1))
# RDFPExp240H <- rep(0, length(diff1))
# RDFNExp240H <- rep(0, length(diff1))
# RDFNExp240L <- rep(0, length(diff1))
RDFPExp320 <- rep(0, length(diff1))
RDFNExp320 <- rep(0, length(diff1))
RDFPExp320L <- rep(0, length(diff1))
RDFPExp320H <- rep(0, length(diff1))
RDFNExp320H <- rep(0, length(diff1))
RDFNExp320L <- rep(0, length(diff1))

for (i in 1:numbTissue){
  # Create new variable names
  pop.name = paste0("diff",i)
  #pop.nameN =  paste0("diffN",i)
  pop.nameR = paste0("diffRand",i)
  #pop.nameRN = paste0("diffNRand",i)
  
  pop.namePExp30 = paste0("diffExp30_",i)
  #pop.nameNExp30 =  paste0("diffNExp30_",i)
  pop.namePExp60 = paste0("diffExp60_",i)
  #pop.nameNExp60 =  paste0("diffNExp60_",i)
  pop.namePExp90 = paste0("diffExp90_",i)
  #pop.nameNExp90 =  paste0("diffNExp90_",i)
  #pop.namePExp120 = paste0("diffExp120_",i)
  #pop.nameNExp120 =  paste0("diffNExp120_",i)
  # pop.namePExp150 = paste0("diffExp150_",i)
  # pop.nameNExp150 =  paste0("diffNExp150_",i)
  pop.namePExp180 = paste0("diffExp180_",i)
  #pop.nameNExp180 =  paste0("diffNExp180_",i)
  #pop.namePExp240 = paste0("diffExp240_",i)
  # pop.nameNExp240 =  paste0("diffNExp240_",i)
  pop.namePExp320 = paste0("diffExp320_",i)
  #pop.nameNExp320 =  paste0("diffNExp320_",i)
  
  RDFPRandtot <- RDFPRandtot + (eval(parse(text = pop.nameR)))/numbTissue
  #RDFNRandtot <- RDFNRandtot + (eval(parse(text = pop.nameRN)))/numbTissue
  
  RDFPtot <- RDFPtot + (eval(parse(text = pop.name)))/numbTissue
  #RDFNtot <- RDFNtot + (eval(parse(text = pop.nameN)))/numbTissue
  RDFPExp30 <- RDFPExp30 + (eval(parse(text = pop.namePExp30)))/numbTissue
  #RDFNExp30 <- RDFNExp30 + (eval(parse(text = pop.nameNExp30)))/numbTissue
  RDFPExp60 <- RDFPExp60 + (eval(parse(text = pop.namePExp60)))/numbTissue
  #RDFNExp60 <- RDFNExp60 + (eval(parse(text = pop.nameNExp60)))/numbTissue
  RDFPExp90 <- RDFPExp90 + (eval(parse(text = pop.namePExp90)))/numbTissue
  #RDFNExp90 <- RDFNExp90 + (eval(parse(text = pop.nameNExp90)))/numbTissue
  #RDFPExp120 <- RDFPExp120 + (eval(parse(text = pop.namePExp120)))/numbTissue
  #RDFNExp120 <- RDFNExp120 + (eval(parse(text = pop.nameNExp120)))/numbTissue
  # RDFPExp150 <- RDFPExp150 + (eval(parse(text = pop.namePExp150)))/numbTissue
  # RDFNExp150 <- RDFNExp150 + (eval(parse(text = pop.nameNExp150)))/numbTissue
  RDFPExp180 <- RDFPExp180 + (eval(parse(text = pop.namePExp180)))/numbTissue
  #RDFNExp180 <- RDFNExp180 + (eval(parse(text = pop.nameNExp180)))/numbTissue
  # RDFPExp240 <- RDFPExp240 + (eval(parse(text = pop.namePExp240)))/numbTissue
  # RDFNExp240 <- RDFNExp240 + (eval(parse(text = pop.nameNExp240)))/numbTissue
  RDFPExp320 <- RDFPExp320 + (eval(parse(text = pop.namePExp320)))/numbTissue
  #RDFNExp320 <- RDFNExp320 + (eval(parse(text = pop.nameNExp320)))/numbTissue
  
}

for (i in 1:numbTissue){
  pop.nameL = paste0("diffL",i)
  #pop.nameNL =  paste0("diffNL",i)
  pop.nameH = paste0("diffH",i)
  pop.nameNH =  paste0("diffNH",i)
  pop.nameRL = paste0("diffRL",i)
  #pop.nameRNL = paste0("diffNRL",i)
  pop.nameRH = paste0("diffRH",i)
  #pop.nameRNH = paste0("diffNRH",i)
  
  RDFPtotL <- RDFPtotL + (eval(parse(text = pop.nameL)))/numbTissue
  #RDFNtotL <- RDFNtotL + (eval(parse(text = pop.nameNL)))/numbTissue
  RDFPRandtotL <- RDFPRandtotL + (eval(parse(text = pop.nameRL)))/numbTissue
  #RDFNRandtotL <- RDFNRandtotL + (eval(parse(text = pop.nameRNL)))/numbTissue
  
  RDFPtotH <- RDFPtotH + (eval(parse(text = pop.nameH)))/numbTissue
  #RDFNtotH <- RDFNtotH + (eval(parse(text = pop.nameNH)))/numbTissue
  RDFPRandtotH <- RDFPRandtotH + (eval(parse(text = pop.nameRH)))/numbTissue
  #RDFNRandtotH <- RDFNRandtotH + (eval(parse(text = pop.nameRNH)))/numbTissue
  
}


###################Plot#####################
xmaxPCF  <- c()
xmaxPCFRand <- c()
ymaxPCF <- c()
ymaxPCFRand <- c()
xmaxPCFN  <- c()
xmaxPCFRandN <- c()
ymaxPCFN <- c()
ymaxPCFRandN <- c()
xmaxPCFExp30  <- c()
ymaxPCFExp30 <- c()
xmaxPCFNExp30  <- c()
ymaxPCFNExp30 <- c()
xmaxPCFExp60  <- c()
ymaxPCFExp60 <- c()
xmaxPCFNExp60  <- c()
ymaxPCFNExp60 <- c()
xmaxPCFExp90  <- c()
ymaxPCFExp90 <- c()
xmaxPCFNExp90  <- c()
ymaxPCFNExp90 <- c()
xmaxPCFExp120  <- c()
ymaxPCFExp120 <- c()
xmaxPCFNExp120  <- c()
ymaxPCFNExp120 <- c()
xmaxPCFExp180  <- c()
ymaxPCFExp180 <- c()
xmaxPCFNExp180  <- c()
ymaxPCFNExp180 <- c()
xmaxPCFExp320  <- c()
ymaxPCFExp320 <- c()
xmaxPCFNExp320  <- c()
ymaxPCFNExp320 <- c()

for (i in 1:numbTissue){
  # Create new variable names
  pop.name = paste0("diff",i)
  pop.nameR = paste0("diffRand",i)
  # pop.nameN =  paste0("diffN",i)
  # pop.nameRN = paste0("diffNRand",i)
  
  pop.namePExp30 = paste0("diffExp30_",i)
  #pop.nameNExp30 =  paste0("diffNExp30_",i)
  pop.namePExp60 = paste0("diffExp60_",i)
  #pop.nameNExp60 =  paste0("diffNExp60_",i)
  pop.namePExp90 = paste0("diffExp90_",i)
  #pop.nameNExp90 =  paste0("diffNExp90_",i)
  # pop.namePExp120 = paste0("diffExp120_",i)
  # pop.nameNExp120 =  paste0("diffNExp120_",i)
  # pop.namePExp150 = paste0("diffExp150_",i)
  # pop.nameNExp150 =  paste0("diffNExp150_",i)
  pop.namePExp180 = paste0("diffExp180_",i)
  #pop.nameNExp180 =  paste0("diffNExp180_",i)
  # pop.namePExp240 = paste0("diffExp240_",i)
  # pop.nameNExp240 =  paste0("diffNExp240_",i)
  pop.namePExp320 = paste0("diffExp320_",i)
  #pop.nameNExp320 =  paste0("diffNExp320_",i)
  
  xmaxPCF  <- c(xmaxPCF , max(eval(parse(text = pop.name))[1:200],na.rm = TRUE))
  xmaxPCFRand  <- c(xmaxPCFRand , max(eval(parse(text = pop.nameR))[1:200],na.rm = TRUE))
  ymaxPCF  <- c(ymaxPCF , RDFPRand1$r[which.max(eval(parse(text = pop.name))[1:200])])
  ymaxPCFRand  <- c(ymaxPCFRand , RDFPRand1$r[which.max(eval(parse(text = pop.nameR))[1:200])])
  #xmaxPCFN  <- c(xmaxPCFN , max(eval(parse(text = pop.nameN))[1:200],na.rm = TRUE))
  #xmaxPCFRandN  <- c(xmaxPCFRandN , max(eval(parse(text = pop.nameRN))[1:200],na.rm = TRUE))
  #ymaxPCFN  <- c(ymaxPCFN , RDFPRand1$r[which.max(eval(parse(text = pop.nameN))[1:200])])
  #ymaxPCFRandN  <- c(ymaxPCFRandN , RDFPRand1$r[which.max(eval(parse(text = pop.nameRN))[1:200])])
  
  xmaxPCFExp30  <- c(xmaxPCFExp30 , max(eval(parse(text = pop.namePExp30))[1:200],na.rm = TRUE))
  ymaxPCFExp30  <- c(ymaxPCFExp30 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp30))[1:200])])
  # xmaxPCFNExp30  <- c(xmaxPCFNExp30 , max(eval(parse(text = pop.nameNExp30))[1:200],na.rm = TRUE))
  # ymaxPCFNExp30  <- c(ymaxPCFNExp30 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp30))[1:200])])
  
  xmaxPCFExp60  <- c(xmaxPCFExp60 , max(eval(parse(text = pop.namePExp60))[1:200],na.rm = TRUE))
  ymaxPCFExp60  <- c(ymaxPCFExp60 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp60))[1:200])])
  # xmaxPCFNExp60  <- c(xmaxPCFNExp60 , max(eval(parse(text = pop.nameNExp60))[1:200],na.rm = TRUE))
  # ymaxPCFNExp60  <- c(ymaxPCFNExp60 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp60))[1:200])])
  
  xmaxPCFExp90  <- c(xmaxPCFExp90 , max(eval(parse(text = pop.namePExp90))[1:200],na.rm = TRUE))
  ymaxPCFExp90  <- c(ymaxPCFExp90 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp90))[1:200])])
  # xmaxPCFNExp90  <- c(xmaxPCFNExp90 , max(eval(parse(text = pop.nameNExp90))[1:200],na.rm = TRUE))
  # ymaxPCFNExp90  <- c(ymaxPCFNExp90 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp90))[1:200])])
  
  
  xmaxPCFExp180  <- c(xmaxPCFExp180 , max(eval(parse(text = pop.namePExp180))[1:200],na.rm = TRUE))
  ymaxPCFExp180  <- c(ymaxPCFExp180 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp180))[1:200])])
  # xmaxPCFNExp180  <- c(xmaxPCFNExp180 , max(eval(parse(text = pop.nameNExp180))[1:200],na.rm = TRUE))
  # ymaxPCFNExp180  <- c(ymaxPCFNExp180 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp180))[1:200])])
  
  xmaxPCFExp320  <- c(xmaxPCFExp320 , max(eval(parse(text = pop.namePExp320))[1:200],na.rm = TRUE))
  ymaxPCFExp320  <- c(ymaxPCFExp320 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp320))[1:200])])
  # xmaxPCFNExp320  <- c(xmaxPCFNExp320 , max(eval(parse(text = pop.nameNExp320))[1:200],na.rm = TRUE))
  # ymaxPCFNExp320  <- c(ymaxPCFNExp320 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp320))[1:200])])
  # 
}

# xmaxPCF  <- c(xmaxPCF , max(eval(parse(text = pop.name))[1:200],na.rm = TRUE))
# xmaxPCFRand  <- c(xmaxPCFRand , max(eval(parse(text = pop.nameR))[1:200],na.rm = TRUE))
# ymaxPCF  <- c(ymaxPCF , RDFPRand1$r[which.max(eval(parse(text = pop.name))[1:200])])
# ymaxPCFRand  <- c(ymaxPCFRand , RDFPRand1$r[which.max(eval(parse(text = pop.nameR))[1:200])])
# xmaxPCFN  <- c(xmaxPCFN , max(eval(parse(text = pop.nameN))[1:200],na.rm = TRUE))
# xmaxPCFRandN  <- c(xmaxPCFRandN , max(eval(parse(text = pop.nameRN))[1:200],na.rm = TRUE))
# ymaxPCFN  <- c(ymaxPCFN , RDFPRand1$r[which.max(eval(parse(text = pop.nameN))[1:200])])
# ymaxPCFRandN  <- c(ymaxPCFRandN , RDFPRand1$r[which.max(eval(parse(text = pop.nameRN))[1:200])])
#
# xmaxPCFExp30  <- c(xmaxPCFExp30 , max(eval(parse(text = pop.namePExp30))[1:200],na.rm = TRUE))
# ymaxPCFExp30  <- c(ymaxPCFExp30 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp30))[1:200])])
# xmaxPCFNExp30  <- c(xmaxPCFNExp30 , max(eval(parse(text = pop.nameNExp30))[1:200],na.rm = TRUE))
# ymaxPCFNExp30  <- c(ymaxPCFNExp30 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp30))[1:200])])
#
# xmaxPCFExp60  <- c(xmaxPCFExp60 , max(eval(parse(text = pop.namePExp60))[1:200],na.rm = TRUE))
# ymaxPCFExp60  <- c(ymaxPCFExp60 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp60))[1:200])])
# xmaxPCFNExp60  <- c(xmaxPCFNExp60 , max(eval(parse(text = pop.nameNExp60))[1:200],na.rm = TRUE))
# ymaxPCFNExp60  <- c(ymaxPCFNExp60 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp60))[1:200])])
#
# xmaxPCFExp90  <- c(xmaxPCFExp90 , max(eval(parse(text = pop.namePExp90))[1:200],na.rm = TRUE))
# ymaxPCFExp90  <- c(ymaxPCFExp90 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp90))[1:200])])
# xmaxPCFNExp90  <- c(xmaxPCFNExp90 , max(eval(parse(text = pop.nameNExp90))[1:200],na.rm = TRUE))
# ymaxPCFNExp90  <- c(ymaxPCFNExp90 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp90))[1:200])])
#
# xmaxPCFExp120  <- c(xmaxPCFExp120 , max(eval(parse(text = pop.namePExp120))[1:200],na.rm = TRUE))
# ymaxPCFExp120  <- c(ymaxPCFExp120 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp120))[1:200])])
# xmaxPCFNExp120  <- c(xmaxPCFNExp120 , max(eval(parse(text = pop.nameNExp120))[1:200],na.rm = TRUE))
# ymaxPCFNExp120  <- c(ymaxPCFNExp120 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp120))[1:200])])
#
# xmaxPCFExp180  <- c(xmaxPCFExp180 , max(eval(parse(text = pop.namePExp180))[1:200],na.rm = TRUE))
# ymaxPCFExp180  <- c(ymaxPCFExp180 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp180))[1:200])])
# xmaxPCFNExp180  <- c(xmaxPCFNExp180 , max(eval(parse(text = pop.nameNExp180))[1:200],na.rm = TRUE))
# ymaxPCFNExp180  <- c(ymaxPCFNExp180 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp180))[1:200])])
#
# xmaxPCFExp320  <- c(xmaxPCFExp320 , max(eval(parse(text = pop.namePExp320))[1:200],na.rm = TRUE))
# ymaxPCFExp320  <- c(ymaxPCFExp320 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp320))[1:200])])
# xmaxPCFNExp320  <- c(xmaxPCFNExp320 , max(eval(parse(text = pop.nameNExp320))[1:200],na.rm = TRUE))
# ymaxPCFNExp320  <- c(ymaxPCFNExp320 , RDFPRand1$r[which.max(eval(parse(text = pop.nameNExp320))[1:200])])
#


##################################################RDF####################
pdf(file= "RFDtot.pdf")
#plot(RDFP1$r[10:400], RDFPtot[10:400],  col = "black",
#     ylab = ("g(r)"), font.lab = 2, cex.lab = 1.3,
#     xlab = expression(plain(paste("r"," ", "(", mu,"m",")"))), main = "PCF Difference")
#dev.off()

#################individual ######
#test7 <- shapiro.test(xmaxPCFExp150)
#test9 <- shapiro.test(xmaxPCFExp240)

WCPG <- t.test(xmaxPCF, xmaxPCFRand, paired = TRUE, alternative = "two.sided" )
#WCPN <- t.test(xmaxPCF, xmaxPCFN, paired = TRUE, alternative = "two.sided" )
#WCPNR <- t.test(xmaxPCFN, xmaxPCFRandN, paired = TRUE, alternative = "two.sided" )
#WCPNRR <- t.test(xmaxPCFRand, xmaxPCFN, paired = TRUE, alternative = "two.sided" )
WCPGExp30 <- t.test(xmaxPCF, xmaxPCFExp30, paired = TRUE, alternative = "two.sided" )
#WCPGExp60 <- t.test(xmaxPCF, xmaxPCFExp60, paired = TRUE, alternative = "two.sided" )
WCPGExp90 <- t.test(xmaxPCF, xmaxPCFExp90, paired = TRUE, alternative = "two.sided" )
#WCPGExp120 <- t.test(xmaxPCF, xmaxPCFExp120, paired = TRUE, alternative = "two.sided" )
#WCPGExp150 <- t.test(xmaxPCF, xmaxPCFExp150, paired = TRUE, alternative = "two.sided" )
WCPGExp180 <- t.test(xmaxPCF, xmaxPCFExp180, paired = TRUE, alternative = "two.sided" )
#WCPGExp240 <- t.test(xmaxPCF, xmaxPCFExp240, paired = TRUE, alternative = "two.sided" )
WCPGExp320 <- t.test(xmaxPCF, xmaxPCFExp320, paired = TRUE, alternative = "two.sided" )

myAsterics <- function(myPval){
  if (myPval >= 0.05){return ("ns")}
  else if((myPval <= 0.05) & (myPval > 0.005)){return("*")}
  else if((myPval <= 0.005) & (myPval > 0.0005)){return("**")}
  else {return("***")}
}

# Places asterisks above each significant box (p < 0.05). The vertical gap is a
# fixed fraction of the y-axis range (frac), so the spacing is visually identical
# across all figures regardless of their y-scale. No bracket lines.
add_sig_stars <- function(p_list, data_list, tested_x, ylim, frac = 0.05) {
  offset <- frac * (ylim[2] - ylim[1])
  layers <- list()
  for (nm in names(p_list)) {
    pval <- p_list[[nm]]$p.value
    if (pval < 0.05) {
      y_star <- max(data_list[[nm]], na.rm = TRUE) + offset
      layers <- c(layers, list(
        annotate("text", x = tested_x[[nm]], y = y_star,
                 label = myAsterics(pval), size = 8, fontface = "bold")
      ))
    }
  }
  layers
}

data <- data.frame(Category = factor(rep(c("BrdU+ Max Bias",  "BrdU+ Exp30", "BrdU+ Exp90",
                                           "BrdU+ Exp180", "BrdU+ Exp320", "BrdU+ Random"), each = numbTissue),
                                     levels = c("BrdU+ Max Bias",  "BrdU+ Exp30", "BrdU+ Exp90",
                                                "BrdU+ Exp180", "BrdU+ Exp320","BrdU+ Random")),
                   Value = c("BrdU+ Max Bias" = xmaxPCF,
                             "BrdU+ Exp30" = xmaxPCFExp30, "BrdU+ Exp90" = xmaxPCFExp90,
                             "BrdU+ Exp180" = xmaxPCFExp180, "BrdU+ Exp320" = xmaxPCFExp320,
                             "BrdU+ Random" = xmaxPCFRand))  # Trim C to match length of A

WCPG <- t.test(xmaxPCF, xmaxPCFRand, paired = TRUE, alternative = "two.sided" )
WCPGR30 <- t.test(xmaxPCFRand, xmaxPCFExp30, paired = TRUE, alternative = "two.sided" )
WCPGR90 <- t.test(xmaxPCFRand, xmaxPCFExp90, paired = TRUE, alternative = "two.sided" )
WCPGR180 <- t.test(xmaxPCFRand, xmaxPCFExp180, paired = TRUE, alternative = "two.sided" )
WCPGR320 <- t.test(xmaxPCFRand, xmaxPCFExp320, paired = TRUE, alternative = "two.sided" )

# x positions: 1=Max Bias, 2=Exp30, 3=Exp90, 4=Exp180, 5=Exp320, 6=Random
pP_list <- list(
  "BrdU+ Max Bias" = list(p.value = WCPG$p.value),
  "BrdU+ Exp30"    = list(p.value = WCPGR30$p.value),
  "BrdU+ Exp90"    = list(p.value = WCPGR90$p.value),
  "BrdU+ Exp180"   = list(p.value = WCPGR180$p.value),
  "BrdU+ Exp320"   = list(p.value = WCPGR320$p.value))

tested_xP <- c("BrdU+ Max Bias" = 1, "BrdU+ Exp30" = 2, "BrdU+ Exp90" = 3,
               "BrdU+ Exp180" = 4, "BrdU+ Exp320" = 5)

# Create the ggplot
bxV <- data$Value[is.finite(data$Value)]
bxR <- max(bxV) - min(bxV); if (bxR == 0) bxR <- 1
yLoBox <- min(bxV) - 0.05 * bxR    # pad below the lowest point
yHiBox <- max(bxV) + 0.20 * bxR    # headroom above the tallest box for the asterisk + text
ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.5, alpha = 0.9, shape = 1) +
  scale_x_discrete(limits = c("BrdU+ Max Bias",  "BrdU+ Exp30", "BrdU+ Exp90",
                              "BrdU+ Exp180", "BrdU+ Exp320", "BrdU+ Random"), labels = function(x) sub("^BrdU[+/-]+ ", "", x)) +
  scale_fill_manual(values = c("grey", "grey", "grey", "grey", "grey", "grey")) +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "black")) +
  labs(y = expression(plain(paste("Maximum RDF")))) +
  coord_cartesian(ylim = c(yLoBox, yHiBox)) +  # Increase Y-axis limit for annotations
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = margin(14, 16, 22, 40, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 32, face = "plain"),  # Bold X-axis labels
    axis.text.y = element_text(size = 32, face = "plain"),  # Bold Y-axis labels
    axis.title.x = element_blank(),#element_text(size = 14, face = "bold"),  # Bold X-axis title
    axis.title.y = element_text(size = 32, face = "plain"),  # Bold Y-axis title
    legend.position = "none"  # Remove legend if not needed
  ) +
  add_sig_stars(pP_list,
                data_list = list("BrdU+ Max Bias" = xmaxPCF,
                                 "BrdU+ Exp30"    = xmaxPCFExp30,
                                 "BrdU+ Exp90"    = xmaxPCFExp90,
                                 "BrdU+ Exp180"   = xmaxPCFExp180,
                                 "BrdU+ Exp320"   = xmaxPCFExp320),
                tested_x = tested_xP,
                ylim = c(yLoBox, yHiBox))
#+
# # Add horizontal comparison lines
# geom_segment(aes(x = 1, xend = 2, y = 4.5, yend = 4.5), color = "black", size = 0.5) +
# geom_segment(aes(x = 2, xend = 3, y = 4.9, yend = 4.9), color = "black", size = 0.5) +
# geom_segment(aes(x = 2, xend = 4, y = 5.3, yend = 5.3), color = "black", size = 0.5) +
# geom_segment(aes(x = 2, xend = 5, y = 5.7, yend = 5.7), color = "black", size = 0.5) +
# geom_segment(aes(x = 2, xend = 6, y = 6.1, yend = 6.1), color = "black", size = 0.5) +
# # Add p-value asterisks above the lines
# geom_text(aes(x = (1+2)/2, y = 4.7, label = myAsterics(round(WCPG$p.value, 3))), size =4, fontface = "bold") +
# geom_text(aes(x = (2+3)/2, y = 5.0, label = myAsterics(round(WCPGR30$p.value, 3))), size = 4, fontface = "bold") +
# geom_text(aes(x = (2+4)/2, y = 5.4, label = myAsterics(round(WCPGR90$p.value, 3))), size = 4, fontface = "bold") +
# geom_text(aes(x = (2+5)/2, y = 5.8, label = myAsterics(round(WCPGR180$p.value, 3))), size = 4, fontface = "bold") +
# geom_text(aes(x = (2+6)/2, y = 6.3, label = myAsterics(round(WCPGR320$p.value, 3))), size = 4, fontface = "bold")

# Add p-value asterisks above the lines
#geom_text(aes(x = 2, y = 4.7, label = myAsterics(round(WCPG$p.value, 3))), size =5, fontface = "bold") +
#  geom_text(aes(x = 3, y = 4.7, label = myAsterics(round(WCPGR30$p.value, 3))), size = 5, fontface = "bold") +
#  geom_text(aes(x = 4, y = 4.7, label = myAsterics(round(WCPGR90$p.value, 3))), size = 5, fontface = "bold") +
#  geom_text(aes(x = 5, y = 4.7, label = myAsterics(round(WCPGR180$p.value, 3))), size = 5, fontface = "bold") +
#  geom_text(aes(x = 6, y = 4.7, label = myAsterics(round(WCPGR320$p.value, 3))), size = 5, fontface = "bold")

ggsave("gXmaxggP.eps", device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("gXmaxggP.pdf", device = "pdf",     dpi = 1200, width = 10, height = 10)

######### curves
# Compute density for each dataset
df <- data.frame(
  x = c(RDFP1$r,
        RDFP1$r, #density(NNPExp60)$x,
        RDFP1$r, #density(NNPExp120)$x, #density(NNPExp150)$x,
        RDFP1$r, #density(NNPExp240)$x,
        RDFP1$r,
        RDFP1$r),#, density(NNP0$x)$x),
  y = c(RDFPtot,
        RDFPExp30, #density(NNPExp60)$y,
        RDFPExp90, #density(NNPExp120)$y, #density(NNPExp150)$y,
        RDFPExp180, #density(NNPExp240)$y,
        RDFPExp320,
        RDFPRandtot
  ),#, density(NNP0$x)$y)
  group = rep(c("BrdU+ Max Bias",
                "BrdU+ Exp30",
                #"NNPExp60",
                "BrdU+ Exp90",
                # "NNPExp120",
                # "NNPExp150",
                "BrdU+ Exp180",
                # "NNPExp240",
                "BrdU+ Exp320",
                "BrdU+ Random"
  ),#, "NNP0"),
  each = length(RDFPtot)))

df$group <- factor(df$group, levels = c(
  "BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
  "BrdU+ Exp180",  "BrdU+ Exp320", "BrdU+ Random"))

pink_gradient <- c(
  "#7B2D8B",  # Medium purple
  "#B5006A",  # Dark magenta-rose
  "#E8003A",  # Bright crimson
  "#FF7800",  # Vivid orange
  "#FF99CC"   # Pale pink
)

names(pink_gradient) <- c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                          "BrdU+ Exp180",  "BrdU+ Exp320")
color_mapping <- c(
  pink_gradient,
  "BrdU+ Random" = "black"
)

line_mapping <- c(
  "BrdU+ Max Bias" = "solid",
  "BrdU+ Exp30"    = "solid",
  "BrdU+ Exp90"    = "solid",
  "BrdU+ Exp180"   = "solid",
  "BrdU+ Exp320"   = "solid",
  "BrdU+ Random"   = "solid"
)

# Horizontal reference lines
# Create the ggplot
## dynamic y-range: fit all condition curves + the null envelope within the x-window (+5% pad)
xLo <- 0; xHi <- myMergeXmax()
inWin <- df$x >= xLo & df$x <= xHi
inR   <- RDFP1$r >= xLo & RDFP1$r <= xHi
yAll  <- c(df$y[inWin], RDFPRandtotL[inR], RDFPRandtotH[inR]); yAll <- yAll[is.finite(yAll)]
yPad  <- 0.05 * (max(yAll) - min(yAll)); yLo <- min(yAll) - yPad; yHi <- max(yAll) + yPad

ggplot(df, aes(x = x, y = y, color = group, linetype = group)) +
  # random-labelling null envelope (identical for every condition -> shown once, in grey)
  geom_ribbon(data = data.frame(x = RDFP1$r, ymin = RDFPRandtotL, ymax = RDFPRandtotH),
              aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
              fill = "grey60", alpha = 0.25) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = line_mapping) +
  labs(
    #title = "Radial Distribution Function",
    x = expression(plain(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
    y  = expression(plain(paste("Radial Distribution Function")))
  ) +
  coord_cartesian(ylim = c(yLo, yHi), xlim = c(xLo, xHi)) +  # Increase Y-axis limit for annotations
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    axis.ticks.length = unit(0.25, "cm"),
    text = element_text(size = 32),
    plot.title = element_text(face = "italic", hjust = 0.5),
    axis.title = element_text(face = "plain", size = 32),
    plot.margin = margin(14, 16, 22, 40, "pt"),
    legend.position = "none"
  ) +
  guides(
    color = guide_legend(title = ""),
    linetype = guide_legend(title = "")
  )

ggsave(filename = "RDFCombAll.eps",
       device = cairo_ps,
       dpi = 1200,
       width = 10,
       height = 10,
       family = "Arial")
ggsave(filename = "RDFCombAll.pdf",
       device = "pdf",
       dpi = 1200,
       width = 10,
       height = 10)
