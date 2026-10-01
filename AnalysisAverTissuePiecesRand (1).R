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
numberList <- seq(from = 1, to = numbTissue, by = 1)

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
assign(paste("RDFP", number, sep = ""),  (read.csv(file = paste(i,"_myenvRDFP3D30_5",".csv",sep=""))))
assign(paste("NNP", number, sep = ""),  (read.csv(file = paste(i,"_myNNP3D30_5",".csv",sep=""))))
assign(paste("RhoP", number, sep = ""),  (read.csv(file = paste(i,"_myRhoP3D30_5",".csv",sep=""))))
assign(paste("GcrossP", number, sep = ""),  (read.csv(file = paste(i,"_myenvGcrossP3D30_5",".csv",sep=""))))
assign(paste("JcrossP", number, sep = ""),  (read.csv(file = paste(i,"_myenvJcrossP3D30_5",".csv",sep=""))))
assign(paste("LcrossP", number, sep = ""),  (read.csv(file = paste(i,"_myenvLcrossP3D30_5",".csv",sep=""))))

assign(paste("RDFPRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPR30_5",".csv",sep=""))))
assign(paste("NNPRand", number, sep = ""),  (read.csv(file = paste(number,"_myNNPR30_5",".csv",sep=""))))
assign(paste("RhoPRand", number, sep = ""),  (read.csv(file = paste(number,"_myRhoPR30_5",".csv",sep=""))))
assign(paste("GcrossPRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvGcrossPR30_5",".csv",sep=""))))
assign(paste("JcrossPRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPR30_5",".csv",sep=""))))
assign(paste("LcrossPRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvLcrossPR30_5",".csv",sep=""))))

assign(paste("RDFPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp30_5",".csv",sep=""))))
assign(paste("NNPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp30_5",".csv",sep=""))))
assign(paste("RhoPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myRhoPExp30_5",".csv",sep=""))))
assign(paste("GcrossPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvGcrossPExp30_5",".csv",sep=""))))
assign(paste("JcrossPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp30_5",".csv",sep=""))))
assign(paste("LcrossPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvLcrossPExp30_5",".csv",sep=""))))

assign(paste("RDFPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp60_5",".csv",sep=""))))
assign(paste("NNPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp60_5",".csv",sep=""))))
assign(paste("RhoPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myRhoPExp60_5",".csv",sep=""))))
assign(paste("GcrossPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvGcrossPExp60_5",".csv",sep=""))))
assign(paste("JcrossPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp60_5",".csv",sep=""))))
assign(paste("LcrossPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvLcrossPExp60_5",".csv",sep=""))))

assign(paste("RDFPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp90_5",".csv",sep=""))))
assign(paste("NNPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp90_5",".csv",sep=""))))
assign(paste("RhoPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myRhoPExp90_5",".csv",sep=""))))
assign(paste("GcrossPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvGcrossPExp90_5",".csv",sep=""))))
assign(paste("JcrossPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp90_5",".csv",sep=""))))
assign(paste("LcrossPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvLcrossPExp90_5",".csv",sep=""))))



assign(paste("RDFPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp180_5",".csv",sep=""))))
assign(paste("NNPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp180_5",".csv",sep=""))))
assign(paste("RhoPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myRhoPExp180_5",".csv",sep=""))))
assign(paste("GcrossPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvGcrossPExp180_5",".csv",sep=""))))
assign(paste("JcrossPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp180_5",".csv",sep=""))))
assign(paste("LcrossPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvLcrossPExp180_5",".csv",sep=""))))


assign(paste("RDFPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvRDFPExp320_5",".csv",sep=""))))
assign(paste("NNPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp320_5",".csv",sep=""))))
assign(paste("RhoPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myRhoPExp320_5",".csv",sep=""))))
assign(paste("GcrossPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvGcrossPExp320_5",".csv",sep=""))))
assign(paste("JcrossPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp320_5",".csv",sep=""))))
assign(paste("LcrossPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvLcrossPExp320_5",".csv",sep=""))))

}

imyNumb <- 5
assign(paste("RDFP", 0, sep = ""),  (read.csv(file = paste(imyNumb,"_myenvRDFP",".csv",sep=""))))
assign(paste("NNP", 0, sep = ""),  (read.csv(file = paste(imyNumb,"_myNNP",".csv",sep=""))))
assign(paste("RhoP", 0, sep = ""),  (read.csv(file = paste(imyNumb,"_myRhoP",".csv",sep=""))))
assign(paste("GcrossP", 0, sep = ""),  (read.csv(file = paste(imyNumb,"_myenvGcrossP",".csv",sep=""))))
assign(paste("JcrossP", 0, sep = ""),  (read.csv(file = paste(imyNumb,"_myenvJcrossP",".csv",sep=""))))
assign(paste("LcrossP", 0, sep = ""),  (read.csv(file = paste(imyNumb,"_myenvLcrossP",".csv",sep=""))))

JcrossP0
########## getting data for distribution to the nearest neighbour
NNPtotal <- c()
NNPRandtotal <- c()
NNPExp60 <- c()
NNPExp30 <- c()
NNPExp90 <- c()
NNPExp180 <- c()
NNPExp320 <- c()


for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("NNP",i)
  pop.nameR = paste0("NNPRand",i)
  pop.namePExp60 = paste0("NNPExp60_",i)
  pop.namePExp30 = paste0("NNPExp30_",i)
  pop.namePExp90 = paste0("NNPExp90_",i)
  pop.namePExp180 = paste0("NNPExp180_",i)
  pop.namePExp320 = paste0("NNPExp320_",i)
  
  NNPtotal <-  c(NNPtotal , eval(parse(text = pop.name))$x)
  NNPRandtotal <-  c(NNPRandtotal , eval(parse(text = pop.nameR))$x)
  NNPExp30 <-  c(NNPExp30 , eval(parse(text = pop.namePExp30))$x)
  NNPExp60 <-  c(NNPExp60 , eval(parse(text = pop.namePExp60))$x)
  NNPExp90 <-  c(NNPExp90 , eval(parse(text = pop.namePExp90))$x)
  NNPExp180 <-  c(NNPExp180 , eval(parse(text = pop.namePExp180))$x)
  NNPExp320 <-  c(NNPExp320 , eval(parse(text = pop.namePExp320))$x)
}

write.table(NNPtotal, file = paste(0,"_NNPtotal_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(NNPRandtotal, file = paste(0,"_NNPRandtotal_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
##################
######## medians of the distribution
mediansP  <- c()
mediansPRand <- c()

mediansPExp30  <- c()
mediansPExp60  <- c()
mediansPExp90  <- c()
mediansPExp180  <- c()
mediansPExp320  <- c()


for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("NNP",i)
  pop.nameR = paste0("NNPRand",i)
  pop.namePExp60 = paste0("NNPExp60_",i)
  pop.namePExp30 = paste0("NNPExp30_",i)
  pop.namePExp90 = paste0("NNPExp90_",i)
  pop.namePExp180 = paste0("NNPExp180_",i)
  pop.namePExp320 = paste0("NNPExp320_",i)
  
  mediansP  <- c(mediansP , median(eval(parse(text = pop.name))$x))
  mediansPRand  <- c(mediansPRand , median(eval(parse(text = pop.nameR))$x))
  mediansPExp30  <- c(mediansPExp30 , median(eval(parse(text = pop.namePExp30))$x))
  mediansPExp60  <- c(mediansPExp60 , median(eval(parse(text = pop.namePExp60))$x))
  mediansPExp90  <- c(mediansPExp90 , median(eval(parse(text = pop.namePExp90))$x))
  mediansPExp180  <- c(mediansPExp180 , median(eval(parse(text = pop.namePExp180))$x))
  mediansPExp320  <- c(mediansPExp320 , median(eval(parse(text = pop.namePExp320))$x))
  
}

write.table(mediansP, file = paste(0,"_mediansP_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(mediansPRand, file = paste(0,"_mediansPRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
#################
######### peaks of the distribution
xmaxP <- c()
xmaxPRand <- c()

xmaxPExp30 <- c()
xmaxPExp60 <- c()
xmaxPExp90 <- c()
xmaxPExp180 <- c()
xmaxPExp320 <- c()

for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("NNP",i)
  pop.nameR = paste0("NNPRand",i)
  
  pop.namePExp60 = paste0("NNPExp60_",i)
  pop.namePExp30 = paste0("NNPExp30_",i)
  pop.namePExp90 = paste0("NNPExp90_",i)
  pop.namePExp180 = paste0("NNPExp180_",i)
  pop.namePExp320 = paste0("NNPExp320_",i)
  
  xmaxP <- c(xmaxP,density(eval(parse(text = pop.name))$x)$x[which.max(density(eval(parse(text = pop.name))$x)$y)])
  xmaxPRand <- c(xmaxPRand,density(eval(parse(text = pop.nameR))$x)$x[which.max(density(eval(parse(text = pop.nameR))$x)$y)])
  xmaxPExp30 <- c(xmaxPExp30,density(eval(parse(text = pop.namePExp30))$x)$x[which.max(density(eval(parse(text = pop.namePExp30))$x)$y)])
  xmaxPExp60 <- c(xmaxPExp60,density(eval(parse(text = pop.namePExp60))$x)$x[which.max(density(eval(parse(text = pop.namePExp60))$x)$y)])
  xmaxPExp90 <- c(xmaxPExp90,density(eval(parse(text = pop.namePExp90))$x)$x[which.max(density(eval(parse(text = pop.namePExp90))$x)$y)])
  xmaxPExp180 <- c(xmaxPExp180,density(eval(parse(text = pop.namePExp180))$x)$x[which.max(density(eval(parse(text = pop.namePExp180))$x)$y)])
  xmaxPExp320 <- c(xmaxPExp320,density(eval(parse(text = pop.namePExp320))$x)$x[which.max(density(eval(parse(text = pop.namePExp320))$x)$y)])
  
}

write.table(xmaxP, file = paste(0,"_xmaxP_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(xmaxPRand, file = paste(0,"_xmaxPRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
#################
################### standard deviations of the distributions


for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("NNP",i)
  pop.nameR = paste0("NNPRand",i)
  
  pop.namePExp60 = paste0("NNPExp60_",i)
  pop.namePExp30 = paste0("NNPExp30_",i)
  pop.namePExp90 = paste0("NNPExp90_",i)
  pop.namePExp180 = paste0("NNPExp180_",i)
  pop.namePExp320 = paste0("NNPExp320_",i)
  
  
  
  
}

####################
############## kolgomorov-smirnoff stats of the distributions
ksStat  <- c()
ksPValue <- c()
ksStatRand  <- c()
ksPValueRand <- c()

ksStatExp30  <- c()
ksPValueExp30  <- c()
ksStatExp60  <- c()
ksPValueExp60 <- c()
ksStatExp90  <- c()
ksPValueExp90 <- c()
ksStatExp180  <- c()
ksPValueExp180 <- c()
ksStatExp320  <- c()
ksPValueExp320 <- c()

ksStatRExp30 <- c()
ksStatRExp60 <- c()
ksStatRExp90 <- c()
ksStatRExp180 <- c()
ksStatRExp320 <- c()

ksStatBias   <- c()
ksStatRandom <- c()

for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("NNP",i)
  pop.nameR = paste0("NNPRand",i)
  pop.namePExp60 = paste0("NNPExp60_",i)
  pop.namePExp30 = paste0("NNPExp30_",i)
  pop.namePExp90 = paste0("NNPExp90_",i)
  pop.namePExp180 = paste0("NNPExp180_",i)
  pop.namePExp320 = paste0("NNPExp320_",i)
 
  
  ksStatBias   <- c(ksStatBias,   round(ks.test(eval(parse(text = pop.name))$x,  NNPRandtotal)$statistic,  digits = 5))  # BrdU+ complete-bias vs random
  ksStatRandom <- c(ksStatRandom, round(ks.test(eval(parse(text = pop.nameR))$x, NNPRandtotal)$statistic, digits = 5))  # BrdU+ random vs random (floor)
  ksStatRand <- c(ksStatRand, round(ks.test(NNP0$x, eval(parse(text = pop.nameR))$x)$statistic, digits = 5))
  #ksPValueRand <- c(ksPValueRand, round(ks.test(eval(parse(text = NNP0))$x, eval(parse(text = pop.nameR))$x)$p.value, digits = 5))
  
  #ksPValueN <- c(ksPValue, round(ks.test(eval(parse(text = NNP0))$x, eval(parse(text = pop.nameN))$x)$p.value, digits = 5))
  
  #ksPValueRandN <- c(ksPValueRand, round(ks.test(eval(parse(text = NNP0))$x, eval(parse(text = pop.nameRN))$x)$p.value, digits = 5))
  
  ksStat <- c(ksStat, round(ks.test(NNP0$x, eval(parse(text = pop.name))$x)$statistic, digits = 5))
  #ksPValue <- c(ksPValue, round(ks.test(eval(parse(text = NNP0))$x, eval(parse(text = pop.name))$x)$p.value, digits = 5))
  
  #ksPValueRandNN <- c(ksPValueRand, round(ks.test(eval(parse(text = NNN0))$x, eval(parse(text = pop.nameR))$x)$p.value, digits = 5))
  
  #ksPValueNN <- c(ksPValue, round(ks.test(eval(parse(text = NNN0))$x, eval(parse(text = pop.nameN))$x)$p.value, digits = 5))
  
  #ksPValueRandNN <- c(ksPValueRand, round(ks.test(eval(parse(text = NNN0))$x, eval(parse(text = pop.nameRN))$x)$p.value, digits = 5))
  
  #ksPValueNN <- c(ksPValue, round(ks.test(eval(parse(text = NNN0))$x, eval(parse(text = pop.name))$x)$p.value, digits = 5))

  
  ksStatRExp30 <- c(ksStatRExp30, round(ks.test(eval(parse(text = pop.namePExp30))$x, eval(parse(text = pop.nameR))$x)$statistic, digits = 5))
  
  ksStatRExp60 <- c(ksStatRExp60, round(ks.test(eval(parse(text = pop.namePExp60))$x, eval(parse(text = pop.nameR))$x)$statistic, digits = 5))
  
  ksStatRExp90 <- c(ksStatRExp90, round(ks.test(eval(parse(text = pop.namePExp90))$x, eval(parse(text = pop.nameR))$x)$statistic, digits = 5))
  
  
  
  ksStatRExp180 <- c(ksStatRExp180, round(ks.test(eval(parse(text = pop.namePExp180))$x, eval(parse(text = pop.nameR))$x)$statistic, digits = 5))
 
  
  ksStatRExp320 <- c(ksStatRExp320, round(ks.test(eval(parse(text = pop.namePExp320))$x, eval(parse(text = pop.nameR))$x)$statistic, digits = 5))
  
}

write.table(ksStat, file = paste(0,"_KS_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(ksStatRand, file = paste(0,"_KSRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(ksPValue, file = paste(0,"_KSPV_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(ksPValueRand, file = paste(0,"_KSPVRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

#********************************* RDF ***********************************
x <- RDFP1$r
for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("RDFP",i)
  pop.nameR = paste0("RDFPRand",i)
  pop.namePExp30 = paste0("RDFPExp30_",i)
  pop.namePExp60 = paste0("RDFPExp60_",i)
  pop.namePExp90 = paste0("RDFPExp90_",i)
  pop.namePExp180 = paste0("RDFPExp180_",i)
  pop.namePExp320 = paste0("RDFPExp320_",i)
  
  
  assign(paste0("diffRand", i),((eval(parse(text = pop.nameR))$obs)/(eval(parse(text = pop.nameR))$mmean)))
  assign(paste0("diffRL", i),((eval(parse(text = pop.nameR))$lo)/(eval(parse(text = pop.nameR))$mmean)))
  assign(paste0("diffRH", i),((eval(parse(text = pop.nameR))$hi)/(eval(parse(text = pop.nameR))$mmean)))
  
  
  assign(paste0("diff", i),((eval(parse(text = pop.name))$obs)/(eval(parse(text = pop.name))$mmean)))
  assign(paste0("diffL", i),((eval(parse(text = pop.name))$lo)/(eval(parse(text = pop.name))$mmean)))
  assign(paste0("diffH", i),((eval(parse(text = pop.name))$hi)/(eval(parse(text = pop.name))$mmean)))
  
  assign(paste0("diffExp30_", i),((eval(parse(text = pop.namePExp30))$obs)/(eval(parse(text = pop.namePExp30))$mmean)))
  assign(paste0("diffLExp30_", i),((eval(parse(text = pop.namePExp30))$lo)/(eval(parse(text = pop.namePExp30))$mmean)))
  assign(paste0("diffHExp30_", i),((eval(parse(text = pop.namePExp30))$hi)/(eval(parse(text = pop.namePExp30))$mmean)))
  
  assign(paste0("diffExp60_", i),((eval(parse(text = pop.namePExp60))$obs)/(eval(parse(text = pop.namePExp60))$mmean)))
  assign(paste0("diffLExp60_", i),((eval(parse(text = pop.namePExp60))$lo)/(eval(parse(text = pop.namePExp60))$mmean)))
  assign(paste0("diffHExp60_", i),((eval(parse(text = pop.namePExp60))$hi)/(eval(parse(text = pop.namePExp60))$mmean)))
  
  assign(paste0("diffExp90_", i),((eval(parse(text = pop.namePExp90))$obs)/(eval(parse(text = pop.namePExp90))$mmean)))
  assign(paste0("diffLExp90_", i),((eval(parse(text = pop.namePExp90))$lo)/(eval(parse(text = pop.namePExp90))$mmean)))
  assign(paste0("diffHExp90_", i),((eval(parse(text = pop.namePExp90))$hi)/(eval(parse(text = pop.namePExp90))$mmean)))
  
  
  
  assign(paste0("diffExp180_", i),((eval(parse(text = pop.namePExp180))$obs)/(eval(parse(text = pop.namePExp180))$mmean)))
  assign(paste0("diffLExp180_", i),((eval(parse(text = pop.namePExp180))$lo)/(eval(parse(text = pop.namePExp180))$mmean)))
  assign(paste0("diffHExp180_", i),((eval(parse(text = pop.namePExp180))$hi)/(eval(parse(text = pop.namePExp180))$mmean)))
  
  
  assign(paste0("diffExp320_", i),((eval(parse(text = pop.namePExp320))$obs)/(eval(parse(text = pop.namePExp320))$mmean)))
  assign(paste0("diffLExp320_", i),((eval(parse(text = pop.namePExp320))$lo)/(eval(parse(text = pop.namePExp320))$mmean)))
  assign(paste0("diffHExp320_", i),((eval(parse(text = pop.namePExp320))$hi)/(eval(parse(text = pop.namePExp320))$mmean)))
  
}

RDFPRandtot <- rep(0, length(diff1))
RDFPRandtotL <- rep(0, length(diff1))
RDFPRandtotH <- rep(0, length(diff1))

RDFPtot <- rep(0, length(diff1))
RDFPtotL <- rep(0, length(diff1))
RDFPtotH <- rep(0, length(diff1))

RDFPExp30 <- rep(0, length(diff1))
RDFPExp30L <- rep(0, length(diff1))
RDFPExp30H <- rep(0, length(diff1))

RDFPExp60 <- rep(0, length(diff1))
RDFPExp60L <- rep(0, length(diff1))
RDFPExp60H <- rep(0, length(diff1))

RDFPExp90 <- rep(0, length(diff1))
RDFPExp90L <- rep(0, length(diff1))
RDFPExp90H <- rep(0, length(diff1))



RDFPExp180 <- rep(0, length(diff1))
RDFPExp180L <- rep(0, length(diff1))
RDFPExp180H <- rep(0, length(diff1))


RDFPExp320 <- rep(0, length(diff1))
RDFPExp320L <- rep(0, length(diff1))
RDFPExp320H <- rep(0, length(diff1))

for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("diff",i)
  pop.nameR = paste0("diffRand",i)
  
  pop.namePExp30 = paste0("diffExp30_",i)
  pop.namePExp60 = paste0("diffExp60_",i)
  pop.namePExp90 = paste0("diffExp90_",i)
  pop.namePExp180 = paste0("diffExp180_",i)
  pop.namePExp320 = paste0("diffExp320_",i)
  
  RDFPRandtot <- RDFPRandtot + (eval(parse(text = pop.nameR)))/numbTissue
  
  RDFPtot <- RDFPtot + (eval(parse(text = pop.name)))/numbTissue
  RDFPExp30 <- RDFPExp30 + (eval(parse(text = pop.namePExp30)))/numbTissue
  RDFPExp60 <- RDFPExp60 + (eval(parse(text = pop.namePExp60)))/numbTissue
  RDFPExp90 <- RDFPExp90 + (eval(parse(text = pop.namePExp90)))/numbTissue
  RDFPExp180 <- RDFPExp180 + (eval(parse(text = pop.namePExp180)))/numbTissue
  RDFPExp320 <- RDFPExp320 + (eval(parse(text = pop.namePExp320)))/numbTissue
  
}
  
for (i in 1:numbTissue){ 
  pop.nameL = paste0("diffL",i)
  pop.nameH = paste0("diffH",i)
  pop.nameRL = paste0("diffRL",i)
  pop.nameRH = paste0("diffRH",i)
  
  RDFPtotL <- RDFPtotL + (eval(parse(text = pop.nameL)))/numbTissue
  RDFPRandtotL <- RDFPRandtotL + (eval(parse(text = pop.nameRL)))/numbTissue
  
  RDFPtotH <- RDFPtotH + (eval(parse(text = pop.nameH)))/numbTissue
  RDFPRandtotH <- RDFPRandtotH + (eval(parse(text = pop.nameRH)))/numbTissue
  
}

write.table(RDFPtot, file = paste(0,"_RDFPtot_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(RDFPRandtot, file = paste(0,"_RDFPRandtot_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")


###################Plot#####################

xmaxPCF  <- c()
xmaxPCFRand <- c()
ymaxPCF <- c()
ymaxPCFRand <- c()


xmaxPCFExp30  <- c()
ymaxPCFExp30 <- c()

xmaxPCFExp60  <- c()
ymaxPCFExp60 <- c()

xmaxPCFExp90  <- c()
ymaxPCFExp90 <- c()



xmaxPCFExp180  <- c()
ymaxPCFExp180 <- c()


xmaxPCFExp320  <- c()
ymaxPCFExp320 <- c()

for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("diff",i)
  pop.nameR = paste0("diffRand",i)
  
  pop.namePExp30 = paste0("diffExp30_",i)
  pop.namePExp60 = paste0("diffExp60_",i)
  pop.namePExp90 = paste0("diffExp90_",i)
  pop.namePExp180 = paste0("diffExp180_",i)
  pop.namePExp320 = paste0("diffExp320_",i)
  
  xmaxPCF  <- c(xmaxPCF , max(eval(parse(text = pop.name))[1:200],na.rm = TRUE))
  xmaxPCFRand  <- c(xmaxPCFRand , max(eval(parse(text = pop.nameR))[1:200],na.rm = TRUE))
  ymaxPCF  <- c(ymaxPCF , RDFPRand1$r[which.max(eval(parse(text = pop.name))[1:200])])
  ymaxPCFRand  <- c(ymaxPCFRand , RDFPRand1$r[which.max(eval(parse(text = pop.nameR))[1:200])])
  
  xmaxPCFExp30  <- c(xmaxPCFExp30 , max(eval(parse(text = pop.namePExp30))[1:200],na.rm = TRUE))
  ymaxPCFExp30  <- c(ymaxPCFExp30 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp30))[1:200])])
  
  xmaxPCFExp60  <- c(xmaxPCFExp60 , max(eval(parse(text = pop.namePExp60))[1:200],na.rm = TRUE))
  ymaxPCFExp60  <- c(ymaxPCFExp60 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp60))[1:200])])
  
  xmaxPCFExp90  <- c(xmaxPCFExp90 , max(eval(parse(text = pop.namePExp90))[1:200],na.rm = TRUE))
  ymaxPCFExp90  <- c(ymaxPCFExp90 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp90))[1:200])])
  
  
  
  xmaxPCFExp180  <- c(xmaxPCFExp180 , max(eval(parse(text = pop.namePExp180))[1:200],na.rm = TRUE))
  ymaxPCFExp180  <- c(ymaxPCFExp180 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp180))[1:200])])
  
  
  xmaxPCFExp320  <- c(xmaxPCFExp320 , max(eval(parse(text = pop.namePExp320))[1:200],na.rm = TRUE))
  ymaxPCFExp320  <- c(ymaxPCFExp320 , RDFPRand1$r[which.max(eval(parse(text = pop.namePExp320))[1:200])])
  
}

write.table(xmaxPCF, file = paste(0,"_xmaxPCF_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(xmaxPCFRand, file = paste(0,"_xmaxPCFRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(ymaxPCF, file = paste(0,"_ymaxPCF_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(ymaxPCFRand, file = paste(0,"_ymaxPCFRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")


#############################
#*********************************** J function ***********************
x <- JcrossP1$r
for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("JcrossP",i)
  pop.nameR = paste0("JcrossPRand",i)
  
  pop.namePExp30 = paste0("JcrossPExp30_",i)
  pop.namePExp60 = paste0("JcrossPExp60_",i)
  pop.namePExp90 = paste0("JcrossPExp90_",i)
  pop.namePExp180 = paste0("JcrossPExp180_",i)
  pop.namePExp320 = paste0("JcrossPExp320_",i)
  
  assign(paste0("Jdiff", i),((eval(parse(text = pop.name))$obs)/(eval(parse(text = pop.name))$mmean)))
  assign(paste0("JdiffRand", i),((eval(parse(text = pop.nameR))$obs)/(eval(parse(text = pop.nameR))$mmean)))
  
  assign(paste0("JdiffPExp30_", i),((eval(parse(text = pop.namePExp30))$obs)/(eval(parse(text = pop.namePExp30))$mmean)))
  assign(paste0("JdiffPExp60_", i),((eval(parse(text = pop.namePExp60))$obs)/(eval(parse(text = pop.namePExp60))$mmean)))
  assign(paste0("JdiffPExp90_", i),((eval(parse(text = pop.namePExp90))$obs)/(eval(parse(text = pop.namePExp90))$mmean)))
  assign(paste0("JdiffPExp180_", i),((eval(parse(text = pop.namePExp180))$obs)/(eval(parse(text = pop.namePExp180))$mmean)))
  assign(paste0("JdiffPExp320_", i),((eval(parse(text = pop.namePExp320))$obs)/(eval(parse(text = pop.namePExp320))$mmean)))
  
}

JPtot <- rep(0, length(diff1))
JPRandtot <- rep(0, length(diff1))

JPExp30 <- rep(0, length(diff1))
JPExp60 <- rep(0, length(diff1))
JPExp90 <- rep(0, length(diff1))
JPExp180 <- rep(0, length(diff1))
JPExp320 <- rep(0, length(diff1))


for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("Jdiff",i)
  pop.nameR = paste0("JdiffRand",i)
  
  pop.namePExp30 = paste0("JdiffPExp30_",i)
  pop.namePExp60 = paste0("JdiffPExp60_",i)
  pop.namePExp90 = paste0("JdiffPExp90_",i)
  pop.namePExp180 = paste0("JdiffPExp180_",i)
  pop.namePExp320 = paste0("JdiffPExp320_",i)
  
  JPtot <- JPtot + (eval(parse(text = pop.name)))/numbTissue
  JPRandtot <- JPRandtot + (eval(parse(text = pop.nameR)))/numbTissue
  
  JPExp30 <- JPExp30 + (eval(parse(text = pop.namePExp30)))/numbTissue
  JPExp60 <- JPExp60 + (eval(parse(text = pop.namePExp60)))/numbTissue
  JPExp90 <- JPExp90 + (eval(parse(text = pop.namePExp90)))/numbTissue
  JPExp180 <- JPExp180 + (eval(parse(text = pop.namePExp180)))/numbTissue
  JPExp320 <- JPExp320 + (eval(parse(text = pop.namePExp320)))/numbTissue
}

write.table(JPtot, file = paste(0,"_JPtot_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(JPRandtot, file = paste(0,"_JPRandtot_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(JcrossP1$r, file = paste(0,"_xJ_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

xmaxJ  <- c()
xmaxJRand <- c()
xmaxJN  <- c()
xmaxJNRand <- c()

xmaxJPExp30 <- c()
xmaxJNExp30 <- c()
xmaxJPExp60 <- c()
xmaxJNExp60 <- c()
xmaxJPExp90 <- c()
xmaxJNExp90 <- c()
xmaxJPExp180 <- c()
xmaxJNExp180 <- c()
xmaxJPExp320 <- c()
xmaxJNExp320 <- c()

for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("Jdiff",i)
  pop.nameR = paste0("JdiffRand",i)
  
  pop.namePExp30 = paste0("JdiffPExp30_",i)
  pop.namePExp60 = paste0("JdiffPExp60_",i)
  pop.namePExp90 = paste0("JdiffPExp90_",i)
  pop.namePExp180 = paste0("JdiffPExp180_",i)
  pop.namePExp320 = paste0("JdiffPExp320_",i)
  
  xmaxJ  <- c(xmaxJ , eval(parse(text = pop.name))[120])
  xmaxJRand  <- c(xmaxJRand , eval(parse(text = pop.nameR))[120])
  
  xmaxJPExp30  <- c(xmaxJPExp30 , eval(parse(text = pop.namePExp30))[120])
  xmaxJPExp60  <- c(xmaxJPExp60 , eval(parse(text = pop.namePExp60))[120])
  xmaxJPExp90  <- c(xmaxJPExp90 , eval(parse(text = pop.namePExp90))[120])
  xmaxJPExp180  <- c(xmaxJPExp180 , eval(parse(text = pop.namePExp180))[120])
  xmaxJPExp320  <- c(xmaxJPExp320 , eval(parse(text = pop.namePExp320))[120])
  
}

write.table(xmaxJ, file = paste(0,"_xmaxJ_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(xmaxJRand, file = paste(0,"_xmaxJRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(xmaxJN, file = paste(0,"_xmaxJN_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(xmaxJNRand, file = paste(0,"_xmaxJNRand_5",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
#*********************************************
##################################################plotting ########
#########distributions

xmaxp <- which.max(density(NNPtotal)$y)
xmaxpRand <- which.max(density(NNPRandtotal)$y)
xmaxpExp30 <- which.max(density(NNPExp30)$y)
xmaxpExp60 <- which.max(density(NNPExp60)$y)
xmaxpExp90 <- which.max(density(NNPExp90)$y)
xmaxpExp180 <- which.max(density(NNPExp180)$y)
xmaxpExp320 <- which.max(density(NNPExp320)$y)

medianP0 <- median(NNP0$x)
ksP0 <-  round(ks.test(NNP0$x, NNPRandtotal)$statistic, digits = 5)
xmaxP0 <- density(NNP0$x)$x[which.max(density(NNP0$x)$y)]
gmaxP0 <- max((RDFP0$obs/RDFP0$mmean),na.rm = TRUE)
jmaxP0 <- max(JcrossP0$obs,na.rm = TRUE)


# Compute density for each dataset (Exp120, Exp150, Exp240 dropped)
df <- data.frame(
  x = c(density(NNPtotal)$x, density(NNPExp30)$x, density(NNPExp60)$x,
        density(NNPExp90)$x, density(NNPExp180)$x, density(NNPExp320)$x,
        density(NNPRandtotal)$x, density(NNP0$x)$x),
  y = c(density(NNPtotal)$y, density(NNPExp30)$y, density(NNPExp60)$y,
        density(NNPExp90)$y, density(NNPExp180)$y, density(NNPExp320)$y,
        density(NNPRandtotal)$y, density(NNP0$x)$y),
  group = rep(c("NNPtotal", "NNPExp30", "NNPExp60", "NNPExp90",
                "NNPExp180", "NNPExp320", "NNPRandtotal", "NNP0"),
              each = length(density(NNPtotal)$x))
)
df$group <- factor(df$group, levels = c(
  "NNPtotal", "NNPExp30", "NNPExp60", "NNPExp90",
  "NNPExp180", "NNPExp320", "NNPRandtotal", "NNP0"))

# palette D1 (Bias = deep magenta, gradient to light; Sample/Random black)
color_mapping <- c(
  "NNPtotal"     = "#8A0046",   # Bias
  "NNPExp30"     = "#C2006B",
  "NNPExp60"     = "#E80074",
  "NNPExp90"     = "#FF4D6D",
  "NNPExp180"    = "#FF92C7",
  "NNPExp320"    = "#FFC65C",
  "NNPRandtotal" = "black",     # Random (dashed)
  "NNP0"         = "black"      # Sample  (solid)
)
line_mapping <- c(
  "NNPtotal" = "solid", "NNPExp30" = "solid", "NNPExp60" = "solid",
  "NNPExp90" = "solid", "NNPExp180" = "solid", "NNPExp320" = "solid",
  "NNPRandtotal" = "dashed", "NNP0" = "solid"
)

pDistr <- ggplot(df, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = line_mapping) +
  labs(x = expression(plain(paste("Distance from Stroma"," ","(",mu,"m",")"))),
       y = "Probability Density") +
  coord_cartesian(ylim = c(0, 0.035), xlim = c(0, 200)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_text(size = 40, face = "plain"),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  )
print(pDistr)
ggsave("DistrStr.eps", plot = pDistr, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("DistrStr.pdf", plot = pDistr, device = "pdf",      dpi = 1200, width = 10, height = 10)


######################################################################################

WCP <- t.test(mediansP, mediansPRand, paired = TRUE, alternative = "two.sided") #1 3


maxMediansPPRand <- max(round(max(mediansP)), round(max(mediansPRand))) # 1,3
maxAllMedians <- max(round(max(mediansP)),round(max(mediansPRand)))

maxMediansPExp30 <- max(round(max(mediansP)), round(max(mediansPExp30))) # 1,3
maxMediansPExp60 <- max(round(max(mediansP)), round(max(mediansPExp60))) # 1,3
maxMediansPExp90 <- max(round(max(mediansP)), round(max(mediansPExp90))) # 1,3
maxMediansPExp180 <- max(round(max(mediansP)), round(max(mediansPExp180))) # 1,3
maxMediansPExp320 <- max(round(max(mediansP)), round(max(mediansPExp320))) # 1,3


myAsterics <- function(myPval){
  if (myPval >= 0.05){return ("ns")}
  else if((myPval < 0.05) & (myPval >= 0.005)){return("*")}
  else if((myPval < 0.005) & (myPval >= 0.0005)){return("**")}
  else {return("***")}
}


library(ggplot2)

# Data with explicit factor ordering

data <- data.frame(Category = factor(rep(c("Max Bias", "Random",  
                                            "Sample" ), each = 10),
                                     levels = c("Max Bias", "Random",  
                                                 "Sample" )),
                   Value = c("Max Bias" = mediansP,  "Random" = mediansPRand, 
                             
                             "Sample" = c(medianP0, rep(NA, length(mediansP) - length(medianP0))))  # Extend B with NA
                              )  # Trim C to match length of A


# Create the boxplot with asterisks and horizontal lines
ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.8, alpha = 0.9, shape = 1, stroke = 0.7) +
  scale_x_discrete(limits = c("Max Bias", "Random", 
                              "Sample" )) +
  scale_fill_manual(values = c("grey", "grey", "grey")) +  
  scale_color_manual(values = c("black", "black", "black")) +
  labs(y = expression(plain(paste("Medians", " ", "(", mu, "m", ")")))) +
  coord_cartesian(ylim = c(0, 135)) +  # Increase Y-axis limit for annotations
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  ) +
  # Add horizontal comparison lines
 # geom_segment(aes(x = 1, xend = 2, y = 115, yend = 115), color = "black", size = 0.5) +
  #geom_segment(aes(x = 3, xend = 4, y = 115, yend = 115), color = "black", size = 0.5) +
 # geom_segment(aes(x = 1, xend = 3, y = 128, yend = 128), color = "black", size = 0.5) +
 # geom_segment(aes(x = 2, xend = 3, y = 121, yend = 121), color = "black", size = 0.5) +
  # Add p-value asterisks above the lines
  geom_text(aes(x = (1+2)/2, y = 116, label = myAsterics(round(WCP$p.value, 3))), size = 10, fontface = "bold") #+
 # geom_text(aes(x = (3+4)/2, y = 116, label = myAsterics(round(WCNNR$p.value, 3))), size = 4, fontface = "bold") +
  #geom_text(aes(x = (2+3)/2, y = 122, label = myAsterics(round(WCPN$p.value, 3))), size = 4, fontface = "bold") +
  #geom_text(aes(x = (1+3)/2, y = 129, label = myAsterics(round(WCPRN$p.value, 3))), size = 4, fontface = "bold") 

ggsave("MedianPairedgg.eps", device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("MedianPairedgg.pdf", device = "pdf",      dpi = 1200, width = 10, height = 10)
################
################


###############Peaks ##########


WCPPX <- t.test(xmaxP, xmaxPRand, paired = TRUE,alternative = "two.sided") # 1 3

maxXmaxPPRand <- max(round(max(xmaxP)), round(max(xmaxPRand))) # 1,3
maxAllXmax <- max(round(max(xmaxP)),round(max(xmaxPRand)))


maxXmaxPExp30 <- max(round(max(xmaxP)), round(max(xmaxPExp30))) # 1,3
maxXmaxPExp60 <- max(round(max(xmaxP)), round(max(xmaxPExp60))) # 1,3
maxXmaxPExp90 <- max(round(max(xmaxP)), round(max(xmaxPExp90))) # 1,3
maxXmaxPExp180 <- max(round(max(xmaxP)), round(max(xmaxPExp180))) # 1,3
maxXmaxPExp320 <- max(round(max(xmaxP)), round(max(xmaxPExp320))) # 1,3


myAsterics <- function(myPval){
  if (myPval >= 0.05){return ("ns")}
  else if((myPval < 0.05) & (myPval >= 0.005)){return("*")}
  else if((myPval < 0.005) & (myPval >= 0.0005)){return("**")}
  else {return("***")}
}


library(ggplot2)

# Data with explicit factor ordering

data <- data.frame(Category = factor(rep(c("Max Bias", "Random", 
                                            "Sample"), each = 10),
                                     levels = c("Max Bias", "Random",  
                                                 "Sample")),
                   Value = c("Max Bias" = xmaxP,  "Random" = xmaxPRand, 
                              
                             "Sample" = c(xmaxP0, rep(NA, length(xmaxP) - length(xmaxP0)))) ) # Trim C to match length of A


                   

# Create the boxplot with asterisks and horizontal lines
ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.8, alpha = 0.9, shape = 1, stroke = 0.7) +
  scale_x_discrete(limits = c("Max Bias", "Random", 
                              "Sample")) +
  scale_fill_manual(values = c("grey", "grey", "grey")) +  
  scale_color_manual(values = c("black", "black", "black")) +
  labs(y = expression(plain(paste("Peaks", " ", "(", mu, "m", ")")))) +
  coord_cartesian(ylim = c(0, 135)) +  # Increase Y-axis limit for annotations
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  ) +
  # Add horizontal comparison lines
  #geom_segment(aes(x = 1, xend = 2, y = 115, yend = 115), color = "black", size = 0.5) +
  #geom_segment(aes(x = 3, xend = 4, y = 115, yend = 115), color = "black", size = 0.5) +
  #geom_segment(aes(x = 1, xend = 3, y = 128, yend = 128), color = "black", size = 0.5) +
  #geom_segment(aes(x = 2, xend = 3, y = 121, yend = 121), color = "black", size = 0.5) +
  # Add p-value asterisks above the lines
  geom_text(aes(x = (1+2)/2, y = 116, label = myAsterics(round(WCPPX$p.value, 3))), size = 10, fontface = "bold") #+
  #geom_text(aes(x = (3+4)/2, y = 116, label = myAsterics(round(WCPNPX$p.value, 3))), size = 4, fontface = "bold") +
  #geom_text(aes(x = (2+3)/2, y = 122, label = myAsterics(round(WCPNPX$p.value, 3))), size = 4, fontface = "bold") +
  #geom_text(aes(x = (1+3)/2, y = 129, label = myAsterics(round(WCPNPX$p.value, 3))), size = 4, fontface = "bold") 

ggsave("PeakggP.eps", device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("PeakggP.pdf", device = "pdf",      dpi = 1200, width = 10, height = 10)


############################################################################
######KS statistics


# ---- KS between BrdU+ and BrdU- per condition: Bias / Random (10 pieces), Sample (observed) ----
# (ksStatBias, ksStatRandom are accumulated per piece in the KS loop above;
#  ksP0 is the single observed sample value.)

WCPks <- t.test(ksStatBias, ksStatRandom, paired = TRUE, alternative = "two.sided")  # Bias vs Random

ksData <- data.frame(
  Category = factor(rep(c("Max Bias", "Random", "Sample"), each = length(ksStatBias)),
                    levels = c("Max Bias", "Random", "Sample")),
  Value = c(ksStatBias,
            ksStatRandom,
            c(ksP0, rep(NA, length(ksStatBias) - 1)))   # Sample is one observed value
)

pKS <- ggplot(ksData, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.8, alpha = 0.9, shape = 1, stroke = 0.7) +
  scale_fill_manual(values = c("Max Bias" = "grey", "Random" = "grey", "Sample" = "grey")) +
  scale_color_manual(values = c("Max Bias" = "black", "Random" = "black", "Sample" = "black")) +
  labs(y = expression(plain(paste("KS Statistic")))) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  ) +
  geom_text(aes(x = 1.5, y = 0.96, label = myAsterics(round(WCPks$p.value, 3))),
            size = 10, fontface = "bold")   # Bias vs Random
print(pKS)
ggsave("KS_PvsN.eps", plot = pKS, device = cairo_ps, width = 10, height = 10, dpi = 1200, family = "Arial")
ggsave("KS_PvsN.pdf", plot = pKS, device = "pdf",      width = 10, height = 10, dpi = 1200)


#####################RDF###################


#################individual ######


myAsterics <- function(myPval){
  if (myPval >= 0.05){return ("ns")}
  else if((myPval < 0.05) & (myPval >= 0.005)){return("*")}
  else if((myPval < 0.005) & (myPval >= 0.0005)){return("**")}
  else {return("***")}
}

WCPG <- t.test(xmaxPCF, xmaxPCFRand, paired = TRUE, alternative = "two.sided" )

library(ggplot2)

# Data with explicit factor ordering

data <- data.frame(Category = factor(rep(c("Max Bias", "Random",  
                                            "Sample"), each = 10),
                                     levels = c("Max Bias", "Random",
                                                 "Sample")),
                   
                   Value = c("Max Bias" = xmaxPCF,  "Random" = xmaxPCFRand, 
                             
                             "Sample" = c(gmaxP0, rep(NA, length(xmaxPCF) - length(gmaxP0)))  # Extend B with NA
                            
                   ))  # Trim C to match length of A

 

# Create the boxplot with asterisks and horizontal lines
ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.8, alpha = 0.9, shape = 1, stroke = 0.7) +
  scale_x_discrete(limits = c("Max Bias", "Random", 
                               "Sample")) +
  scale_fill_manual(values = c("grey", "grey",  "grey")) +  
  scale_color_manual(values = c("black", "black", "black")) +
  labs(y = expression(plain(paste("Maximum RDF")))) +
  coord_cartesian(ylim = c(0, 7.35)) +  # Increase Y-axis limit for annotations
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  ) +
  # Add horizontal comparison lines
  #geom_segment(aes(x = 1, xend = 2, y = 6.35, yend = 6.35), color = "black", size = 0.9) +
  #geom_segment(aes(x = 3, xend = 4, y = 6.35, yend = 6.35), color = "black", size = 0.5) +
  #geom_segment(aes(x = 1, xend = 3, y = 7.18, yend = 7.18), color = "black", size = 0.9) +
  #geom_segment(aes(x = 2, xend = 3, y = 6.71, yend = 6.71), color = "black", size = 0.9) +
  # Add p-value asterisks above the lines
  geom_text(aes(x = (1+2)/2, y = 6.46, label = myAsterics(round(WCPG$p.value, 3))), size = 10, fontface = "bold") #+
  #geom_text(aes(x = (3+4)/2, y = 6.46, label = myAsterics(round(WCPNR$p.value, 3))), size = 4, fontface = "bold") +
  #geom_text(aes(x = (2+3)/2, y = 6.82, label = myAsterics(round(WCPNRR$p.value, 3))), size = 6, fontface = "bold") +
  #geom_text(aes(x = (1+3)/2, y = 7.29, label = myAsterics(round(WCPNR$p.value, 3))), size = 6, fontface = "bold") 

ggsave("gmaxggP.eps", device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("gmaxggP.pdf", device = "pdf",      dpi = 1200, width = 10, height = 10)


RDFP0tot <- RDFP0$obs/RDFP0$mmean

# Compute density for each dataset
rdf_A <- rbind(
  data.frame(x = RDFP1$r, y = RDFPtot,     group = "RDFPtot"),      # Bias
  data.frame(x = RDFP1$r, y = RDFPExp30,   group = "RDFPExp30"),
  data.frame(x = RDFP1$r, y = RDFPExp60,   group = "RDFPExp60"),
  data.frame(x = RDFP1$r, y = RDFPExp90,   group = "RDFPExp90"),
  data.frame(x = RDFP1$r, y = RDFPExp180,  group = "RDFPExp180"),
  data.frame(x = RDFP1$r, y = RDFPExp320,  group = "RDFPExp320"),
  data.frame(x = RDFP1$r, y = RDFPRandtot, group = "RDFPRandtot"),  # Random
  data.frame(x = RDFP0$r, y = RDFP0tot,    group = "RDFP0tot")      # Sample
)
rdf_A$group <- factor(rdf_A$group, levels = c("RDFPtot","RDFPExp30","RDFPExp60","RDFPExp90",
                                          "RDFPExp180","RDFPExp320","RDFPRandtot","RDFP0tot"))
color_mapping <- c("RDFPtot"="#8A0046","RDFPExp30"="#C2006B","RDFPExp60"="#E80074",
                   "RDFPExp90"="#FF4D6D","RDFPExp180"="#FF92C7","RDFPExp320"="#FFC65C",
                   "RDFPRandtot"="black","RDFP0tot"="black")
line_mapping <- c("RDFPtot"="solid","RDFPExp30"="solid","RDFPExp60"="solid","RDFPExp90"="solid",
                  "RDFPExp180"="solid","RDFPExp320"="solid","RDFPRandtot"="dashed","RDFP0tot"="solid")
pRDFa <- ggplot(rdf_A, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = line_mapping) +
  labs(x = expression(plain(paste("Distance from Stroma"," ","(",mu,"m",")"))),
       y = expression(plain(paste("Radial Distribution Function")))) +
  coord_cartesian(ylim = c(0, 5.05), xlim = c(0, 200)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_text(size = 40, face = "plain"),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  )
print(pRDFa)
ggsave("RDFDiffCom0gg.eps", plot = pRDFa, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("RDFDiffCom0gg.pdf", plot = pRDFa, device = "pdf",      dpi = 1200, width = 10, height = 10)


rdf_B <- rbind(
  data.frame(x = RDFP1$r, y = RDFPtot,     group = "RDFPtot"),      # Bias
  data.frame(x = RDFP1$r, y = RDFPExp30,   group = "RDFPExp30"),
  data.frame(x = RDFP1$r, y = RDFPExp60,   group = "RDFPExp60"),
  data.frame(x = RDFP1$r, y = RDFPExp90,   group = "RDFPExp90"),
  data.frame(x = RDFP1$r, y = RDFPExp180,  group = "RDFPExp180"),
  data.frame(x = RDFP1$r, y = RDFPExp320,  group = "RDFPExp320"),
  data.frame(x = RDFP1$r, y = RDFPRandtot, group = "RDFPRandtot"),  # Random
  data.frame(x = RDFP0$r, y = RDFP0tot,    group = "RDFP0tot")      # Sample
)
rdf_B$group <- factor(rdf_B$group, levels = c("RDFPtot","RDFPExp30","RDFPExp60","RDFPExp90",
                                          "RDFPExp180","RDFPExp320","RDFPRandtot","RDFP0tot"))
color_mapping <- c("RDFPtot"="#8A0046","RDFPExp30"="#C2006B","RDFPExp60"="#E80074",
                   "RDFPExp90"="#FF4D6D","RDFPExp180"="#FF92C7","RDFPExp320"="#FFC65C",
                   "RDFPRandtot"="black","RDFP0tot"="black")
line_mapping <- c("RDFPtot"="solid","RDFPExp30"="solid","RDFPExp60"="solid","RDFPExp90"="solid",
                  "RDFPExp180"="solid","RDFPExp320"="solid","RDFPRandtot"="dashed","RDFP0tot"="solid")
pRDFb <- ggplot(rdf_B, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = line_mapping) +
  labs(x = expression(plain(paste("Distance from Stroma"," ","(",mu,"m",")"))),
       y = expression(plain(paste("Radial Distribution Function")))) +
  coord_cartesian(ylim = c(0, 5.05), xlim = c(0, 200)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_text(size = 40, face = "plain"),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  )
print(pRDFb)
ggsave("RDFCombAll.eps", plot = pRDFb, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("RDFCombAll.pdf", plot = pRDFb, device = "pdf",      dpi = 1200, width = 10, height = 10)


##################################

library(ggplot2)
library(dplyr)

# Shared x-axis
Z <- RhoP1$Z[1:410]

# Combine all line data (Exp120/150/240 dropped)
line_df <- data.frame(
  Z = rep(Z, 8),
  rho = c(RhoP1$rho[1:410], RhoPExp30_1$rho[1:410], RhoPExp60_1$rho[1:410],
          RhoPExp90_1$rho[1:410], RhoPExp180_1$rho[1:410], RhoPExp320_1$rho[1:410],
          RhoPRand1$rho[1:410], RhoP0$rho[1:410]),
  group = rep(c("P1","Exp30","Exp60","Exp90","Exp180","Exp320","Rand","P0"), each = 410)
)
line_df$group <- factor(line_df$group,
                        levels = c("P1","Exp30","Exp60","Exp90","Exp180","Exp320","Rand","P0"))

# Confidence bands: Max Bias / Random / Sample
band_df <- bind_rows(
  data.frame(Z = Z, lo = RhoP1$lo[1:410],     hi = RhoP1$hi[1:410],     group = "P1"),
  data.frame(Z = Z, lo = RhoPRand1$lo[1:410], hi = RhoPRand1$hi[1:410], group = "Rand"),
  data.frame(Z = Z, lo = RhoP0$lo[1:410],     hi = RhoP0$hi[1:410],     group = "P0")
)

# palette D1
color_map <- c(
  "P1"     = "#8A0046",   # Max Bias
  "Exp30"  = "#C2006B",
  "Exp60"  = "#E80074",
  "Exp90"  = "#FF4D6D",
  "Exp180" = "#FF92C7",
  "Exp320" = "#FFC65C",
  "Rand"   = "black",     # Random (dashed)
  "P0"     = "black"      # Sample  (solid)
)
linetype_map <- c(
  "P1" = "solid", "Exp30" = "solid", "Exp60" = "solid", "Exp90" = "solid",
  "Exp180" = "solid", "Exp320" = "solid", "Rand" = "dashed", "P0" = "solid"
)

# Horizontal reference lines (ave[20] per condition)
ref_df <- data.frame(
  y = c(RhoP1$ave[20], RhoPExp30_1$ave[20], RhoPExp60_1$ave[20], RhoPExp90_1$ave[20],
        RhoPExp180_1$ave[20], RhoPExp320_1$ave[20], RhoPRand1$ave[20], RhoP0$ave[20]),
  group = c("P1","Exp30","Exp60","Exp90","Exp180","Exp320","Rand","P0")
)

pRSFline <- ggplot() +
  geom_ribbon(data = band_df, aes(x = Z, ymin = lo, ymax = hi, group = group),
              fill = "grey50", alpha = 0.3) +
  geom_line(data = line_df, aes(x = Z, y = rho, color = group, linetype = group), size = 1.5) +
  geom_hline(data = ref_df, aes(yintercept = y),
             color = "#99004D", linetype = "dashed", size = 1.0) +
  scale_color_manual(values = color_map) +
  scale_linetype_manual(values = linetype_map) +
  labs(x = expression(plain(paste("Distance from Stroma"," ","(",mu,"m",")"))),
       y = "Resource Selection Function") +
  coord_cartesian(ylim = c(-0.0005, 0.0035), xlim = c(0, 200)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_text(size = 40, face = "plain"),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  )
print(pRSFline)
ggsave("RSFgg.eps", plot = pRSFline, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("RSFgg.pdf", plot = pRSFline, device = "pdf",      dpi = 1200, width = 10, height = 10)


# ============ RSF box plot: Bias / Random / Sample ============
# RSF statistic  S = integral |rho/ave - 1| dZ   (ave = rho-bar; strength of resource selection)
rsfStat <- function(rhoObj){
  z <- rhoObj$Z
  v <- abs(rhoObj$rho / rhoObj$ave - 1)
  sum(diff(z) * (head(v, -1) + tail(v, -1)) / 2)   # trapezoidal integral
}

rsfBias <- c(); rsfRandom <- c()
for (i in 1:numbTissue){
  rsfBias   <- c(rsfBias,   rsfStat(eval(parse(text = paste0("RhoP", i)))))
  rsfRandom <- c(rsfRandom, rsfStat(eval(parse(text = paste0("RhoPRand", i)))))
}
rsfP0 <- rsfStat(RhoP0)                               # observed sample (single value)

WCPrsf <- t.test(rsfBias, rsfRandom, paired = TRUE, alternative = "two.sided")  # Bias vs Random

rsfData <- data.frame(
  Category = factor(rep(c("Max Bias", "Random", "Sample"), each = length(rsfBias)),
                    levels = c("Max Bias", "Random", "Sample")),
  Value = c(rsfBias,
            rsfRandom,
            c(rsfP0, rep(NA, length(rsfBias) - 1)))   # Sample is one observed value
)

pRSF <- ggplot(rsfData, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.8, alpha = 0.9, shape = 1, stroke = 0.7) +
  scale_fill_manual(values = c("grey", "grey", "grey")) +
  scale_color_manual(values = c("black", "black", "black")) +
  labs(y = expression(plain(paste("Integrated |", rho, "/", bar(rho), "-1| dz")))) +
  coord_cartesian(ylim = c(0, 410)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  ) +
  geom_text(aes(x = 1.5, y = 392, label = myAsterics(round(WCPrsf$p.value, 3))),
            size = 10, fontface = "bold")   # Bias vs Random
print(pRSF)
ggsave("RSF_Boxplot.eps", plot = pRSF, device = cairo_ps, width = 10, height = 10, dpi = 1200, family = "Arial")
ggsave("RSF_Boxplot.pdf", plot = pRSF, device = "pdf",      width = 10, height = 10, dpi = 1200)

myAsterics <- function(myPval){
  if (myPval >= 0.05){return ("ns")}
  else if((myPval < 0.05) & (myPval >= 0.005)){return("*")}
  else if((myPval < 0.005) & (myPval >= 0.0005)){return("**")}
  else {return("***")}
}

assign(paste0("JPtot", 0),((JcrossP0$obs)/(JcrossP0$mmean)))
JPtot0 = c(JPtot0, rep(NA, length(JPtot) - length(JPtot0)))
length(JPtot0)
xmaxJP0 <- JPtot0[120]

JPtot <- JPtot[0:501]
length(JPtot)
WCPJ <- t.test(xmaxJ, xmaxJRand, paired = TRUE, alternative = "two.sided")

JPExp30 <- JPExp30[0:501]
JPExp60 <- JPExp60[0:501]
JPExp90 <- JPExp90[0:501]
JPExp180 <- JPExp180[0:501] 
JPExp320 <-  JPExp320[0:501]
JPRandtot <-  JPRandtot[0:501]


#### J function (Exp120/150/240 dropped; palette D1; each curve paired with its own r) ####
JPtot0 <- (JcrossP0$obs) / (JcrossP0$mmean)        # sample (un-padded; same 501 r-grid)

jdf <- rbind(
  data.frame(x = JcrossP1$r, y = JPtot,     group = "JPtot"),      # Bias
  data.frame(x = JcrossP1$r, y = JPExp30,   group = "JPExp30"),
  data.frame(x = JcrossP1$r, y = JPExp60,   group = "JPExp60"),
  data.frame(x = JcrossP1$r, y = JPExp90,   group = "JPExp90"),
  data.frame(x = JcrossP1$r, y = JPExp180,  group = "JPExp180"),
  data.frame(x = JcrossP1$r, y = JPExp320,  group = "JPExp320"),
  data.frame(x = JcrossP1$r, y = JPRandtot, group = "JPRandtot"),  # Random
  data.frame(x = JcrossP0$r, y = JPtot0,    group = "JPtot0")      # Sample
)
jdf$group <- factor(jdf$group, levels = c("JPtot","JPExp30","JPExp60","JPExp90",
                                          "JPExp180","JPExp320","JPRandtot","JPtot0"))

color_mapping <- c(
  "JPtot"     = "#8A0046",   # Bias
  "JPExp30"   = "#C2006B",
  "JPExp60"   = "#E80074",
  "JPExp90"   = "#FF4D6D",
  "JPExp180"  = "#FF92C7",
  "JPExp320"  = "#FFC65C",
  "JPRandtot" = "black",     # Random (dashed)
  "JPtot0"    = "black"      # Sample  (solid)
)
line_mapping <- c(
  "JPtot" = "solid", "JPExp30" = "solid", "JPExp60" = "solid", "JPExp90" = "solid",
  "JPExp180" = "solid", "JPExp320" = "solid", "JPRandtot" = "dashed", "JPtot0" = "solid"
)

pJ <- ggplot(jdf, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = line_mapping) +
  labs(x = expression(plain(paste("Distance from Stroma", " ", "(", mu, "m", ")"))),
       y = expression(plain(paste("J Function")))) +
  coord_cartesian(ylim = c(0, 1.6), xlim = c(0, 200)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_text(size = 40, face = "plain"),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  )
print(pJ)
ggsave("JfunctCombAll.eps", plot = pJ, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("JfunctCombAll.pdf", plot = pJ, device = "pdf",      dpi = 1200, width = 10, height = 10)

###############
########
########
data <- data.frame(Category = factor(rep(c("Max Bias", "Random",  
                                           "Sample"), each = 10),
                                     levels = c("Max Bias", "Random",
                                                "Sample")),
                   
                   Value = c("Max Bias" = xmaxJ,  "Random" = xmaxJRand, 
                             
                             "Sample" = c(xmaxJP0, rep(NA, length(xmaxJ) - length(xmaxJP0)))  # Extend B with NA
                             
                   ))  # Trim C to match length of A


# Create the boxplot with asterisks and horizontal lines
ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.8, alpha = 0.9, shape = 1, stroke = 0.7) +
  scale_x_discrete(limits = c("Max Bias", "Random", 
                              "Sample")) +
  scale_fill_manual(values = c("grey", "grey",  "grey")) +  
  scale_color_manual(values = c("black", "black", "black")) +
  labs(y = expression(plain(paste("Jmin")))) +
  coord_cartesian(ylim = c(0, 2.35)) +  # Increase Y-axis limit for annotations
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 2.5),
    axis.line.y = element_line(color = "black", size = 2.5),
    axis.ticks.y = element_line(color = "black", size = 2.5),
    axis.ticks.x = element_line(color = "black", size = 2.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 36, face = "plain", color = "black"),
    axis.text.y = element_text(size = 36, face = "plain", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
    legend.position = "none"
  ) +
  # Add horizontal comparison lines
  #geom_segment(aes(x = 1, xend = 2, y = 6.35, yend = 6.35), color = "black", size = 0.9) +
  #geom_segment(aes(x = 3, xend = 4, y = 6.35, yend = 6.35), color = "black", size = 0.5) +
  #geom_segment(aes(x = 1, xend = 3, y = 7.18, yend = 7.18), color = "black", size = 0.9) +
  #geom_segment(aes(x = 2, xend = 3, y = 6.71, yend = 6.71), color = "black", size = 0.9) +
  # Add p-value asterisks above the lines
  geom_text(aes(x = (1+2)/2, y = 2.26, label = myAsterics(round(WCPJ$p.value, 3))), size = 10, fontface = "bold") #+
#geom_text(aes(x = (3+4)/2, y = 6.46, label = myAsterics(round(WCPNR$p.value, 3))), size = 4, fontface = "bold") +
#geom_text(aes(x = (2+3)/2, y = 6.82, label = myAsterics(round(WCPNRR$p.value, 3))), size = 6, fontface = "bold") +
#geom_text(aes(x = (1+3)/2, y = 7.29, label = myAsterics(round(WCPNR$p.value, 3))), size = 6, fontface = "bold") 

ggsave("JmaxggP.eps", device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("JmaxggP.pdf", device = "pdf",      dpi = 1200, width = 10, height = 10)
