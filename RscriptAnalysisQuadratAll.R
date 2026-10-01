rm(list=ls())
library(sp)
library(rmapshaper)
library(ggplot2)
library(magrittr)
library(sf)
library(spatstat)
library(goftest)
library(readr)
library(kSamples)

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))


mynumber <- paste(readLines("number.txt"), collapse=" ")
number <-as.numeric(mynumber)
numb <- as.numeric(paste(readLines("polNumb.txt"), collapse=" "))
linn <-readLines(file("Rparams.txt",open="r"))
close(file("Rparams.txt",open="r"))
quadratSize <- as.numeric(linn[1]) #1000
PCF_r <- as.numeric(linn[2]) # 450
J_r <- as.numeric(linn[3]) #250
Z_r_length <- as.numeric(linn[4]) #500
dens_r_length <- as.numeric(linn[5]) #500
PCF_r_length <- as.integer(PCF_r*2)
J_r_length <- as.integer(J_r*2)

xmin1 <-0
xmax1 <-quadratSize
ymin1 <-0
ymax1 <-quadratSize

#No changes
rPCF <-seq(from = 0, to = PCF_r, by = 0.5)
rJ <-seq(from = 0, to = J_r, by = 0.5)


nncrossp <- c()
nncrossn <- c()
numR <- 0
numb <- 21
for (i in 1:10){
mydataPos <- read.csv(file = paste("PosR30_",i,"_21.csv",sep=""))
mydataNeg <- read.csv(file = paste("NegR30_",i,"_21.csv",sep=""))
mydataStr <- read.csv(file = paste("StrR30_",i,"_21.csv",sep=""))


#mytotalpattern <- ppp(mydata[,1], mydata[,2], c(xmin1,xmax1), c(ymin1,ymax1), marks=factor(mydata[,3]))
mypatternPos <- ppp(mydataPos[,1], mydataPos[,2], c(xmin1,xmax1), c(ymin1,ymax1))# marks=factor(mydataPos[,3]))
mypatternNeg <- ppp(mydataNeg[,1], mydataNeg[,2], c(xmin1,xmax1), c(ymin1,ymax1))# marks=factor(mydataNeg[,3]))
mypatternStr <- ppp(mydataStr[,1], mydataStr[,2], c(xmin1,xmax1), c(ymin1,ymax1))# marks=factor(mydataStr[,3]))



if ((mypatternStr$n + mypatternPos$n + mypatternNeg$n)> 3000){
  numR <- numR + 1

#rho hat 
Z <- distmap(mypatternStr)
myRho <- rhohat(mypatternNeg,Z)
myRhoP <- rhohat(mypatternPos,Z)



write.table(myRhoP, file = paste(i,"_myRhoPRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(myRho, file = paste(i,"_myRhoNRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
#rho hat plot
# jpeg(file=paste0(i, "_rhohatPos.jpeg"),width = 8, height = 8, units = 'in', res = 300)
# plot(myRhoP, col = "red",main="Rho hat BrdU+",
#      xlab= "Distance to Stroma (um)",
#      ylab= "BrdU+ Density (# /unit area)",
#      font.main=3, font.lab=2, legend = FALSE)
# dev.off()
# jpeg(file=paste0(i,"_rhohatNeg.jpeg"))
# plot(myRho, col = "blue", main="Rho hat BrdU-",
#      xlab= "Distance to Stroma (um)",
#      ylab= "BrdU- Density (# /unit area)",
#      font.main=3, font.lab=2, legend = FALSE)
# dev.off()

nncrossp <- c()
nncrossn <- c()

nn <- nncross(mypatternPos, mypatternStr)$dist
mm <- nncross(mypatternNeg, mypatternStr)$dist
for(jj in 1:length(nn)){
  nncrossp <- c(nncrossp, nn[jj])}
for(hh in 1:length(mm)){
  nncrossn <- c(nncrossn, mm[hh])}
write.table(nn, file = paste(i,"_myNNPRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(mm, file = paste(i,"_myNNNRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

## starting PCF,J,G and L
# plot the pcf or rdf
XY <- superimpose(Pos=mypatternPos,Str=mypatternStr)
XYN <- superimpose(Neg=mypatternNeg,Str=mypatternStr)
PN <- superimpose(Pos=mypatternPos,Neg=mypatternNeg)
mysim <- function(ppp1,ppp2){
  superimpose(PosR = split(rlabel(ppp1))$Pos, Str = ppp2)
}

mysimN <- function(ppp1,ppp2){
  superimpose(PosR = split(rlabel(ppp1))$Neg, Str = ppp2)
}


myenv <- envelope(XYN, fun = pcfcross, correction=c("Ripley"),r = rPCF, nsim=4, simulate=expression(mysimN(PN,mypatternStr)))
myenvp <- envelope(XY, fun = pcfcross, correction=c("Ripley"),nsim=39,r= rPCF, simulate=expression(mysim(PN,mypatternStr)))

write.table(myenvp, file = paste(i,"_myenvRDFPRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(myenv, file = paste(i,"_myenvRDFNRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

# ymax = max(c(max(myenvp$obs[5:PCF_r_length]), max(myenv$obs[5:PCF_r_length])))
# jpeg(file=paste0(i,"_PCFcross.jpeg"),width = 8, height = 8, units = 'in', res = 300)
# plot(myenvp[15:PCF_r_length], main = "gcross BrdU+ and BrdU-", xlab = expression(paste("Distance from Stroma"," ", "(", mu,"m",")")), ylab = expression("g"["BrdU+,Str"]),
#      cex.lab = 1.3, font.lab = 4, legend = FALSE, col = "red", lwd = 3)
# lines(myenv[15:PCF_r_length], col = "blue", lty = 1, lwd = 3, legend = FALSE)
# abline(h = 1)
# dev.off()


#Lcross ####################################################
myenvL <- envelope(XYN, fun = Lcross, correction=c("Ripley"),nsim=4,r = rPCF, simulate=expression(mysimN(PN,mypatternStr)))
myenvpL <- envelope(XY, fun = Lcross, correction=c("Ripley"),nsim=39, r= rPCF, simulate=expression(mysim(PN,mypatternStr)))

write.table(myenvpL, file = paste(i,"_myenvLcrossPRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(myenvL, file = paste(i,"_myenvLcrossNRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
# 
# # plot Lcross
# jpeg(file=paste0(i,"_LcrossPos.jpeg"),width = 8, height = 8, units = 'in', res = 300)
# plot(myenvpL, ((.-r)/r) ~ r,  main = "Lcross BrdU+", xlab = expression(paste("Distance from Stroma"," ", "(", mu,"m",")")), 
#      ylab = expression("L"["BrdU+,Str"]),
#      cex.lab = 1.3, font.lab = 4, legend = FALSE, col = "red", lwd = 3)
# dev.off()
# 
# ymaxN = min(myenvL$obs)
# jpeg(file=paste0(i,"_LcrossNeg.jpeg"),width = 8, height = 8, units = 'in', res = 300)
# plot(myenvL, ((.-r)/r) ~ r, main = "Lcross BrdU-", xlab = expression(paste("r"," ", "(", mu,"m",")")), ylab = expression("L"["BrdU+,Str"]),
#      cex.lab = 1.3, font.lab = 4, legend = FALSE, col = "blue", lwd = 3)
# dev.off()


#Jcross #########################################################
myenvJ <- envelope(XYN, fun = Jcross, correction=c("rs"),nsim=39,r =rJ, simulate=expression(mysimN(PN,mypatternStr)))
myenvpJ <- envelope(XY, fun = Jcross, correction=c("rs"),nsim=39,r = rJ, simulate=expression(mysim(PN,mypatternStr)))

write.table(myenvpJ, file = paste(i,"_myenvJcrossPRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(myenvJ, file = paste(i,"_myenvJcrossNRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

# xmaxJ = max(myenvJ$r)
# ymaxJ = max(myenvpJ$obs, myenvJ$obs)
# jpeg(file=paste0(i,"_Jcross.jpeg"),width = 8, height = 8, units = 'in', res = 300)
# plot(myenvpJ, main = "Jcross BrdU+ and BrdU-", xlab = expression(paste("r"," ", "(", mu,"m",")")), ylab = expression("J"["BrdU+,Str"]),
#      cex.lab = 1.3, font.lab = 4, legend = FALSE,  col = "red", lwd = 3)
# lines(myenvJ, col = "blue", lty = 1, lwd = 3, legend = FALSE)
# abline(h = 1)
# dev.off()


#Gcross #############################################################
myenvG <- envelope(XYN, fun = Gcross, correction=c("rs"),nsim=39,r = rJ,  simulate=expression(mysimN(PN,mypatternStr)))
myenvpG <- envelope(XY, fun = Gcross, correction=c("rs"),nsim=39, r = rJ,simulate=expression(mysim(PN,mypatternStr)))

write.table(myenvpG, file = paste(i,"_myenvGcrossPRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(myenvG, file = paste(i,"_myenvGcrossNRand_21",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
# 
# # plot Gcross
# jpeg(file=paste0(i, "_Gcross.jpeg"),width = 8, height = 8, units = 'in', res = 300)
# plot(myenvpG, main = "Gcross BrdU+ and BrdU-", xlab = expression(paste("r"," ", "(", mu,"m",")")),
#      ylab = expression("G"["BrdU+,Str"]),
#      cex.lab = 1.3, font.lab = 4, legend = FALSE, col = "red", lwd = 3)
# lines(myenvG, col = "blue", lty = 1, lwd = 3)
# dev.off()


}}
