
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
  assign(paste("NNN", number, sep = ""),  (read.csv(file = paste(i,"_myNNN3D_", N, ".csv",sep=""))))
  assign(paste("NNP", number, sep = ""),  (read.csv(file = paste(i,"_myNNP3D_", N, ".csv",sep=""))))
  
  assign(paste("NNNRand", number, sep = ""),  (read.csv(file = paste(number,"_myNNNR_", N, ".csv",sep=""))))
  assign(paste("NNPRand", number, sep = ""),  (read.csv(file = paste(number,"_myNNPR_", N, ".csv",sep=""))))
  
  assign(paste("NNNExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp30_", N, ".csv",sep=""))))
  assign(paste("NNPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp30_", N, ".csv",sep=""))))
  
  assign(paste("NNNExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp60_", N, ".csv",sep=""))))
  assign(paste("NNPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp60_", N, ".csv",sep=""))))
  
  assign(paste("NNNExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp90_", N, ".csv",sep=""))))
  assign(paste("NNPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp90_", N, ".csv",sep=""))))
  
  assign(paste("NNNExp120_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp120_", N, ".csv",sep=""))))
  assign(paste("NNPExp120_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp120_", N, ".csv",sep=""))))
  
  # assign(paste("NNNExp150_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp150_", N, ".csv",sep=""))))
  # assign(paste("NNPExp150_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp150_", N, ".csv",sep=""))))
  
  assign(paste("NNNExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp180_", N, ".csv",sep=""))))
  assign(paste("NNPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp180_", N, ".csv",sep=""))))
  
  # assign(paste("NNNExp240_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp240_", N, ".csv",sep=""))))
  # assign(paste("NNPExp240_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp240_", N, ".csv",sep=""))))
  
  assign(paste("NNNExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myNNNExp320_", N, ".csv",sep=""))))
  assign(paste("NNPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myNNPExp320_", N, ".csv",sep=""))))
  
}
########## getting data for distribution to the nearest neighbour
NNPtotal <- c()
NNNtotal <- c()
NNPRandtotal <- c()
NNNRandtotal <- c()
NNPExp60 <- c()
NNNExp60 <- c()
NNPExp30 <- c()
NNNExp30 <- c()
NNPExp90 <- c()
NNNExp90 <- c()
NNPExp120 <- c()
NNNExp120 <- c()
NNPExp180 <- c()
NNNExp180 <- c()
NNPExp320 <- c()
NNNExp320 <- c()

for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("NNP",i)
  pop.nameN =  paste0("NNN",i)
  pop.nameR = paste0("NNPRand",i)
  pop.nameRN = paste0("NNNRand",i)
  pop.namePExp60 = paste0("NNPExp60_",i)
  pop.nameNExp60 = paste0("NNNExp60_",i)
  pop.namePExp30 = paste0("NNPExp30_",i)
  pop.nameNExp30 = paste0("NNNExp30_",i)
  pop.namePExp90 = paste0("NNPExp90_",i)
  pop.nameNExp90 = paste0("NNNExp90_",i)
  pop.namePExp120 = paste0("NNPExp120_",i)
  pop.nameNExp120 = paste0("NNNExp120_",i)
  # pop.namePExp150 = paste0("NNPExp150_",i)
  # pop.nameNExp150 = paste0("NNNExp150_",i)
  pop.namePExp180 = paste0("NNPExp180_",i)
  pop.nameNExp180 = paste0("NNNExp180_",i)
  # pop.namePExp240 = paste0("NNPExp240_",i)
  # pop.nameNExp240 = paste0("NNNExp240_",i)
  pop.namePExp320 = paste0("NNPExp320_",i)
  pop.nameNExp320 = paste0("NNNExp320_",i)
  
  NNPtotal <-  c(NNPtotal , eval(parse(text = pop.name))$x)
  NNNtotal <-  c(NNNtotal , eval(parse(text = pop.nameN))$x)
  NNPRandtotal <-  c(NNPRandtotal , eval(parse(text = pop.nameR))$x)
  NNNRandtotal <-  c(NNNRandtotal , eval(parse(text = pop.nameRN))$x)
  NNPExp30 <-  c(NNPExp30 , eval(parse(text = pop.namePExp30))$x)
  NNNExp30 <-  c(NNNExp30 , eval(parse(text = pop.nameNExp30))$x)
  NNPExp60 <-  c(NNPExp60 , eval(parse(text = pop.namePExp60))$x)
  NNNExp60 <-  c(NNNExp60 , eval(parse(text = pop.nameNExp60))$x)
  NNPExp90 <-  c(NNPExp90 , eval(parse(text = pop.namePExp90))$x)
  NNNExp90 <-  c(NNNExp90 , eval(parse(text = pop.nameNExp90))$x)
  NNPExp120 <-  c(NNPExp120 , eval(parse(text = pop.namePExp120))$x)
  NNNExp120 <-  c(NNNExp120 , eval(parse(text = pop.nameNExp120))$x)
  # NNPExp150 <-  c(NNPExp150 , eval(parse(text = pop.namePExp150))$x)
  # NNNExp150 <-  c(NNNExp150 , eval(parse(text = pop.nameNExp150))$x)
  NNPExp180 <-  c(NNPExp180 , eval(parse(text = pop.namePExp180))$x)
  NNNExp180 <-  c(NNNExp180 , eval(parse(text = pop.nameNExp180))$x)
  # NNPExp240 <-  c(NNPExp240 , eval(parse(text = pop.namePExp240))$x)
  # NNNExp240 <-  c(NNNExp240 , eval(parse(text = pop.nameNExp240))$x)
  NNPExp320 <-  c(NNPExp320 , eval(parse(text = pop.namePExp320))$x)
  NNNExp320 <-  c(NNNExp320 , eval(parse(text = pop.nameNExp320))$x)
}

plot(density(NNPExp320_5$x))
lines(density(NNP5$x), col = "darkgreen")

##################
xmaxP <- c()
xmaxN <- c()
xmaxPRand <- c()
xmaxNRand <- c()

xmaxPExp30 <- c()
xmaxNExp30 <- c()
xmaxPExp60 <- c()
xmaxNExp60 <- c()
xmaxPExp90 <- c()
xmaxNExp90 <- c()
xmaxPExp120 <- c()
xmaxNExp120 <- c()
xmaxPExp180 <- c()
xmaxNExp180 <- c()
xmaxPExp320 <- c()
xmaxNExp320 <- c()

for (i in 1:numbTissue){ 
  # Create new variable names
  pop.name = paste0("NNP",i)
  pop.nameR = paste0("NNPRand",i)
  pop.nameN =  paste0("NNN",i)
  pop.nameRN = paste0("NNNRand",i)
  
  pop.namePExp60 = paste0("NNPExp60_",i)
  pop.nameNExp60 = paste0("NNNExp60_",i)
  pop.namePExp30 = paste0("NNPExp30_",i)
  pop.nameNExp30 = paste0("NNNExp30_",i)
  pop.namePExp90 = paste0("NNPExp90_",i)
  pop.nameNExp90 = paste0("NNNExp90_",i)
  pop.namePExp120 = paste0("NNPExp120_",i)
  pop.nameNExp120 = paste0("NNNExp120_",i)
  # pop.namePExp150 = paste0("NNPExp150_",i)
  # pop.nameNExp150 = paste0("NNNExp150_",i)
  pop.namePExp180 = paste0("NNPExp180_",i)
  pop.nameNExp180 = paste0("NNNExp180_",i)
  # pop.namePExp240 = paste0("NNPExp240_",i)
  # pop.nameNExp240 = paste0("NNNExp240_",i)
  pop.namePExp320 = paste0("NNPExp320_",i)
  pop.nameNExp320 = paste0("NNNExp320_",i)
  
  xmaxP <- c(xmaxP,density(eval(parse(text = pop.name))$x)$x[which.max(density(eval(parse(text = pop.name))$x)$y)])
  xmaxN <- c(xmaxN,density(eval(parse(text = pop.nameN))$x)$x[which.max(density(eval(parse(text = pop.nameN))$x)$y)])
  xmaxPRand <- c(xmaxPRand,density(eval(parse(text = pop.nameR))$x)$x[which.max(density(eval(parse(text = pop.nameR))$x)$y)])
  xmaxNRand <- c(xmaxNRand,density(eval(parse(text = pop.nameRN))$x)$x[which.max(density(eval(parse(text = pop.nameRN))$x)$y)])
  xmaxPExp30 <- c(xmaxPExp30,density(eval(parse(text = pop.namePExp30))$x)$x[which.max(density(eval(parse(text = pop.namePExp30))$x)$y)])
  xmaxNExp30 <- c(xmaxNExp30,density(eval(parse(text = pop.nameNExp30))$x)$x[which.max(density(eval(parse(text = pop.nameNExp30))$x)$y)])
  xmaxPExp60 <- c(xmaxPExp60,density(eval(parse(text = pop.namePExp60))$x)$x[which.max(density(eval(parse(text = pop.namePExp60))$x)$y)])
  xmaxNExp60 <- c(xmaxNExp60,density(eval(parse(text = pop.nameNExp60))$x)$x[which.max(density(eval(parse(text = pop.nameNExp60))$x)$y)])
  xmaxPExp90 <- c(xmaxPExp90,density(eval(parse(text = pop.namePExp90))$x)$x[which.max(density(eval(parse(text = pop.namePExp90))$x)$y)])
  xmaxNExp90 <- c(xmaxNExp90,density(eval(parse(text = pop.nameNExp90))$x)$x[which.max(density(eval(parse(text = pop.nameNExp90))$x)$y)])
  xmaxPExp120 <- c(xmaxPExp120,density(eval(parse(text = pop.namePExp120))$x)$x[which.max(density(eval(parse(text = pop.namePExp120))$x)$y)])
  xmaxNExp120 <- c(xmaxNExp120,density(eval(parse(text = pop.nameNExp120))$x)$x[which.max(density(eval(parse(text = pop.nameNExp120))$x)$y)])
  # xmaxPExp150 <- c(xmaxPExp150,density(eval(parse(text = pop.namePExp150))$x)$x[which.max(density(eval(parse(text = pop.namePExp150))$x)$y)])
  # xmaxNExp150 <- c(xmaxNExp150,density(eval(parse(text = pop.nameNExp150))$x)$x[which.max(density(eval(parse(text = pop.nameNExp150))$x)$y)])
  xmaxPExp180 <- c(xmaxPExp180,density(eval(parse(text = pop.namePExp180))$x)$x[which.max(density(eval(parse(text = pop.namePExp180))$x)$y)])
  xmaxNExp180 <- c(xmaxNExp180,density(eval(parse(text = pop.nameNExp180))$x)$x[which.max(density(eval(parse(text = pop.nameNExp180))$x)$y)])
  # xmaxPExp240 <- c(xmaxPExp240,density(eval(parse(text = pop.namePExp240))$x)$x[which.max(density(eval(parse(text = pop.namePExp240))$x)$y)])
  # xmaxNExp240 <- c(xmaxNExp240,density(eval(parse(text = pop.nameNExp240))$x)$x[which.max(density(eval(parse(text = pop.nameNExp240))$x)$y)])
  xmaxPExp320 <- c(xmaxPExp320,density(eval(parse(text = pop.namePExp320))$x)$x[which.max(density(eval(parse(text = pop.namePExp320))$x)$y)])
  xmaxNExp320 <- c(xmaxNExp320,density(eval(parse(text = pop.nameNExp320))$x)$x[which.max(density(eval(parse(text = pop.nameNExp320))$x)$y)])
  
}




#################

################## BrdU+ plots
## per-variation-scale bandwidth: SAME for curves, peaks, and null band (define before first use)
bwFix <- bw.nrd0(NNPRandtotal) * numbTissue^(1/5)
xmaxn <- which.max(density(NNNtotal)$y)
xmaxp <- which.max(density(NNPtotal, bw = bwFix)$y)
xmaxnRand <- which.max(density(NNNRandtotal)$y)
xmaxpRand <- which.max(density(NNPRandtotal, bw = bwFix)$y)
xmaxnExp30 <- which.max(density(NNNExp30)$y)
xmaxpExp30 <- which.max(density(NNPExp30, bw = bwFix)$y)
xmaxnExp60 <- which.max(density(NNNExp60)$y)
xmaxpExp60 <- which.max(density(NNPExp60)$y)
xmaxnExp90 <- which.max(density(NNNExp90)$y)
xmaxpExp90 <- which.max(density(NNPExp90, bw = bwFix)$y)
xmaxnExp120 <- which.max(density(NNNExp120)$y)
xmaxpExp120 <- which.max(density(NNPExp120)$y)
# xmaxnExp150 <- which.max(density(NNNExp150)$y)
# xmaxpExp150 <- which.max(density(NNPExp150)$y)
xmaxnExp180 <- which.max(density(NNNExp180)$y)
xmaxpExp180 <- which.max(density(NNPExp180, bw = bwFix)$y)
# xmaxnExp240 <- which.max(density(NNNExp240)$y)
# xmaxpExp240 <- which.max(density(NNPExp240)$y)
xmaxnExp320 <- which.max(density(NNNExp320)$y)
xmaxpExp320 <- which.max(density(NNPExp320, bw = bwFix)$y)
maxY <- max(max(density(NNPtotal, bw = bwFix)$y), max(density(NNNtotal)$y),
            max(density(NNPRandtotal, bw = bwFix)$y),max(density(NNNRandtotal)$y))


#test9 <- shapiro.test(xmaxPExp150)
#test11 <- shapiro.test(xmaxPExp240)

#test17 <- shapiro.test(xmaxNExp150)
#test19 <- shapiro.test(xmaxNExp240)


WCPPX <- t.test(xmaxP, xmaxPRand, paired = FALSE,alternative = "two.sided") # 1 3
WCPNPX <- t.test(xmaxP, xmaxN, alternative = "two.sided") #1 5
WCNNRPX <- t.test(xmaxN, xmaxNRand, paired = FALSE, alternative = "two.sided")#5 7
WCPRNPX <- t.test(xmaxPRand, xmaxN, alternative = "two.sided")# 3 5
WCPNRPX <- t.test(xmaxP, xmaxNRand, alternative = "two.sided")# 1 7

maxXmaxPPRand <- max(round(max(xmaxP)), round(max(xmaxPRand))) # 1,3
maxXmaxPRandN <- max(round(max(xmaxPRand)), round(max(xmaxN))) # 3,5
maxXmaxNNRand <- max(round(max(xmaxN)),round(max(xmaxNRand))) # 5,7
maxXmaxPN <- max(round(max(xmaxP)),round(max(xmaxN))) # 1,5
maxXmaxPNRand <- max(round(max(xmaxP)),round(max(xmaxNRand))) # 1,7
maxAllXmax <- max(round(max(xmaxP)),round(max(xmaxPRand)),
                  round(max(xmaxN)), round(max(xmaxNRand)))

WCPPExp30 <- t.test(xmaxPRand, xmaxPExp30, paired = FALSE,alternative = "two.sided") # 1 3
WCPPExp60 <- t.test(xmaxPRand, xmaxPExp60, paired = FALSE,alternative = "two.sided") # 1 3
WCPPExp90 <- t.test(xmaxPRand, xmaxPExp90, paired = FALSE,alternative = "two.sided") # 1 3
WCPPExp120 <- t.test(xmaxPRand, xmaxPExp120, paired = FALSE,alternative = "two.sided") # 1 3
#WCPPExp150 <- t.test(xmaxPRand, xmaxPExp150, paired = FALSE,alternative = "two.sided") # 1 3
WCPPExp180 <- t.test(xmaxPRand, xmaxPExp180, paired = FALSE,alternative = "two.sided") # 1 3
#WCPPExp240 <- t.test(xmaxPRand, xmaxPExp240, paired = FALSE,alternative = "two.sided") # 1 3
WCPPExp320 <- t.test(xmaxPRand, xmaxPExp320, paired = FALSE,alternative = "two.sided") # 1 3

maxXmaxPExp30 <- max(round(max(xmaxPRand)), round(max(xmaxPExp30))) # 1,3
maxXmaxPExp60 <- max(round(max(xmaxPRand)), round(max(xmaxPExp60))) # 1,3
maxXmaxPExp90 <- max(round(max(xmaxPRand)), round(max(xmaxPExp90))) # 1,3
maxXmaxPExp120 <- max(round(max(xmaxPRand)), round(max(xmaxPExp120))) # 1,3
#maxXmaxPExp150 <- max(round(max(xmaxPRand)), round(max(xmaxPExp150))) # 1,3
maxXmaxPExp180 <- max(round(max(xmaxPRand)), round(max(xmaxPExp180))) # 1,3
#maxXmaxPExp240 <- max(round(max(xmaxPRand)), round(max(xmaxPExp240))) # 1,3
maxXmaxPExp320 <- max(round(max(xmaxPRand)), round(max(xmaxPExp320))) # 1,3


# Compute density for each dataset (bwFix defined above, before its first use)
df <- data.frame(
  x = c(density(NNPtotal, bw = bwFix)$x, 
        density(NNPExp30, bw = bwFix)$x, #density(NNPExp60)$x, 
        density(NNPExp90, bw = bwFix)$x, #density(NNPExp120)$x, #density(NNPExp150)$x,
        density(NNPExp180, bw = bwFix)$x, #density(NNPExp240)$x,
        density(NNPExp320, bw = bwFix)$x,
        density(NNPRandtotal, bw = bwFix)$x),#, density(NNP0$x)$x),
  y = c(density(NNPtotal, bw = bwFix)$y,
        density(NNPExp30, bw = bwFix)$y, #density(NNPExp60)$y, 
        density(NNPExp90, bw = bwFix)$y, #density(NNPExp120)$y, #density(NNPExp150)$y,
        density(NNPExp180, bw = bwFix)$y, #density(NNPExp240)$y, 
        density(NNPExp320, bw = bwFix)$y,
        density(NNPRandtotal, bw = bwFix)$y),#, density(NNP0$x)$y)
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
  each = length(density(NNPtotal, bw = bwFix)$x)))


df$group <- factor(df$group, levels = c(
  "BrdU+ Max Bias",
  "BrdU+ Exp30",
  #"NNPExp60",
  "BrdU+ Exp90",
  #"NNPExp120",
  #"NNPExp150",
  "BrdU+ Exp180",
  # "NNPExp240",
  "BrdU+ Exp320",
  "BrdU+ Random"))#, "NNP0" 
#)
#)


pink_gradient <- c(
  "#7B2D8B",  # Medium purple
  "#B5006A",  # Dark magenta-rose
  "#E8003A",  # Bright crimson
  "#FF7800",  # Vivid orange
  "#FF99CC"   # Pale pink
)
names(pink_gradient) <- c("BrdU+ Max Bias",
                          "BrdU+ Exp30",
                          # "NNPExp60",
                          "BrdU+ Exp90",
                          #  "NNPExp120",
                          #"NNPExp150",
                          "BrdU+ Exp180",
                          # "NNPExp240",
                          "BrdU+ Exp320")

color_mapping <- c(
  pink_gradient,
  "BrdU+ Random" = "black"
)

line_mapping <- c(
  "BrdU+ Max Bias" = "solid",
  "BrdU+ Exp30" = "solid",
  "BrdU+ Exp90" = "solid",
  "BrdU+ Exp180" = "solid",
  "BrdU+ Exp320" = "solid",
  "BrdU+ Random" = "solid")
#)
# Horizontal reference lines
# Create the ggplot
## null band for the distance frequency distribution: 2.5-97.5% spread of the Random-
## condition distance density across the 99 variations (Random labelling = the null)
gx <- density(NNPRandtotal, bw = bwFix)$x
randDensMat <- sapply(1:numbTissue, function(k){
  di <- eval(parse(text = paste0("NNPRand", k)))$x
  density(di, bw = bwFix, from = min(gx), to = max(gx), n = length(gx))$y
})
nnFreqLo <- apply(randDensMat, 1, quantile, probs = 0.025, na.rm = TRUE)
nnFreqHi <- apply(randDensMat, 1, quantile, probs = 0.975, na.rm = TRUE)

## dynamic y-range: fit all condition curves + the null band within the x-window (+5% pad)
xLo <- 0; xHi <- myMergeXmax()
inWin <- df$x >= xLo & df$x <= xHi
inG   <- gx >= xLo & gx <= xHi
yAll  <- c(df$y[inWin], nnFreqLo[inG], nnFreqHi[inG]); yAll <- yAll[is.finite(yAll)]
yPad  <- 0.05 * (max(yAll) - min(yAll)); yLo <- min(yAll) - yPad; yHi <- max(yAll) + yPad

ggplot(df, aes(x = x, y = y, color = group, linetype = group)) +
  # random-labelling null band (spread of the Random distance density across variations)
  geom_ribbon(data = data.frame(x = gx, ymin = nnFreqLo, ymax = nnFreqHi),
              aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
              fill = "grey60", alpha = 0.25) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = line_mapping) +
  labs(
    #title = "Distribution of Distances to the Nearest Stroma",
    x = expression(plain(paste("Distance from Stroma", " ", "(", mu, "m", ")"))),
    y = "Probability Density"
  ) +
  coord_cartesian(xlim = c(xLo, xHi), ylim = c(yLo, yHi)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),  # Keep X-axis ticks
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

ggsave(filename = "DistrStr.eps",
       device = cairo_ps,
       dpi = 1200,
       width = 10,
       height = 10,
       family = "Arial")

ggsave(filename = "DistrStr.pdf",
       device = "pdf",
       dpi = 1200,
       width = 10,
       height = 10)

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

# BrdU+ t-tests: each Exp condition vs Random
WCPPExp30  <- t.test(xmaxPRand, xmaxPExp30,  paired = FALSE, alternative = "two.sided")
WCPPExp60  <- t.test(xmaxPRand, xmaxPExp60,  paired = FALSE, alternative = "two.sided")
WCPPExp90  <- t.test(xmaxPRand, xmaxPExp90,  paired = FALSE, alternative = "two.sided")
WCPPExp120 <- t.test(xmaxPRand, xmaxPExp120, paired = FALSE, alternative = "two.sided")
WCPPExp180 <- t.test(xmaxPRand, xmaxPExp180, paired = FALSE, alternative = "two.sided")
WCPPExp320 <- t.test(xmaxPRand, xmaxPExp320, paired = FALSE, alternative = "two.sided")

# Data with explicit factor ordering

data <- data.frame(Category = factor(rep(c("BrdU+ Max Bias",  "BrdU+ Exp30", "BrdU+ Exp90",
                                           "BrdU+ Exp180", "BrdU+ Exp320", "BrdU+ Random"), each = numbTissue),
                                     levels = c("BrdU+ Max Bias",  "BrdU+ Exp30", "BrdU+ Exp90",
                                                "BrdU+ Exp180", "BrdU+ Exp320", "BrdU+ Random")),
                   Value = c("BrdU+ Max Bias" = xmaxP,
                             "BrdU+ Exp30" = xmaxPExp30, "BrdU+ Exp90" = xmaxPExp90,
                             "BrdU+ Exp180" = xmaxPExp180, "BrdU+ Exp320" = xmaxPExp320,
                             "BrdU+ Random" = xmaxPRand
                   )
)  # Trim C to match length of A


#positions <- list(c(1, 3), c(1, 5), c(5, 7), c(3, 5))  # X-axis pairs for comparison
#y_pos <- c(115, 118, 121, 124, 127)  # Y positions for the lines & text

# x positions: 1=Max Bias, 2=Exp30, 3=Exp90, 4=Exp180, 5=Exp320, 6=Random
pP_list <- list(
  "BrdU+ Max Bias" = list(p.value = t.test(xmaxPRand, xmaxP, paired = FALSE, alternative = "two.sided")$p.value),
  "BrdU+ Exp30"    = list(p.value = WCPPExp30$p.value),
  "BrdU+ Exp90"    = list(p.value = WCPPExp90$p.value),
  "BrdU+ Exp180"   = list(p.value = WCPPExp180$p.value),
  "BrdU+ Exp320"   = list(p.value = WCPPExp320$p.value)
)
tested_xP <- c("BrdU+ Max Bias" = 1, "BrdU+ Exp30" = 2, "BrdU+ Exp90" = 3,
               "BrdU+ Exp180" = 4, "BrdU+ Exp320" = 5)

bxV <- data$Value[is.finite(data$Value)]
bxR <- max(bxV) - min(bxV); if (bxR == 0) bxR <- 1
yLoBox <- min(bxV) - 0.05 * bxR    # pad below the lowest point
yHiBox <- max(bxV) + 0.20 * bxR    # headroom above the tallest box for the asterisk + text
p_brduP <- ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.5, alpha = 0.9, shape = 1) +
  scale_x_discrete(limits = c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                              "BrdU+ Exp180", "BrdU+ Exp320", "BrdU+ Random"), labels = function(x) sub("^BrdU[+/-]+ ", "", x)) +
  scale_fill_manual(values = c("grey", "grey", "grey", "grey", "grey", "grey")) +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "black")) +
  labs(y = expression(plain(paste("Peaks", " ", "(", mu, "m", ")")))) +
  coord_cartesian(ylim = c(yLoBox, yHiBox)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = margin(14, 16, 22, 40, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 32, face = "plain"),
    axis.text.y = element_text(size = 32, face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 32, face = "plain"),
    legend.position = "none"
  ) +
  add_sig_stars(pP_list,
                data_list = list("BrdU+ Max Bias" = xmaxP,
                                 "BrdU+ Exp30"    = xmaxPExp30,
                                 "BrdU+ Exp90"    = xmaxPExp90,
                                 "BrdU+ Exp180"   = xmaxPExp180,
                                 "BrdU+ Exp320"   = xmaxPExp320),
                tested_x = tested_xP,
                ylim = c(yLoBox, yHiBox))

print(p_brduP)
ggsave("PeakggP.eps", plot = p_brduP, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("PeakggP.pdf", plot = p_brduP, device = "pdf", dpi = 1200, width = 10, height = 10)

##################
WCNNExp30  <- t.test(xmaxNRand, xmaxNExp30,  paired = FALSE, alternative = "two.sided")
WCNNExp60  <- t.test(xmaxNRand, xmaxNExp60,  paired = FALSE, alternative = "two.sided")
WCNNExp90  <- t.test(xmaxNRand, xmaxNExp90,  paired = FALSE, alternative = "two.sided")
WCNNExp120 <- t.test(xmaxNRand, xmaxNExp120, paired = FALSE, alternative = "two.sided")
WCNNExp180 <- t.test(xmaxNRand, xmaxNExp180, paired = FALSE, alternative = "two.sided")
WCNNExp320 <- t.test(xmaxNRand, xmaxNExp320, paired = FALSE, alternative = "two.sided")

data <- data.frame(Category = factor(rep(c("BrdU- Max Bias",  "BrdU- Exp30", "BrdU- Exp90",
                                           "BrdU- Exp180", "BrdU- Exp320", "BrdU- Random"), each = numbTissue),
                                     levels = c("BrdU- Max Bias", "BrdU- Exp30", "BrdU- Exp90",
                                                "BrdU- Exp180", "BrdU- Exp320", "BrdU- Random")),
                   Value = c("BrdU- Max Bias" = xmaxN,
                             "BrdU- Exp30" = xmaxNExp30, "BrdU- Exp90" = xmaxNExp90,
                             "BrdU- Exp180" = xmaxNExp180, "BrdU- Exp320" = xmaxNExp320,
                             "BrdU- Random" = xmaxNRand
                   )
)

# x positions: 1=Max Bias, 2=Exp30, 3=Exp90, 4=Exp180, 5=Exp320, 6=Random
pN_list <- list(
  "BrdU- Max Bias" = list(p.value = t.test(xmaxNRand, xmaxN, paired = FALSE, alternative = "two.sided")$p.value),
  "BrdU- Exp30"    = list(p.value = WCNNExp30$p.value),
  "BrdU- Exp90"    = list(p.value = WCNNExp90$p.value),
  "BrdU- Exp180"   = list(p.value = WCNNExp180$p.value),
  "BrdU- Exp320"   = list(p.value = WCNNExp320$p.value)
)
tested_xN <- c("BrdU- Max Bias" = 1, "BrdU- Exp30" = 2, "BrdU- Exp90" = 3,
               "BrdU- Exp180" = 4, "BrdU- Exp320" = 5)

bxV <- data$Value[is.finite(data$Value)]
bxR <- max(bxV) - min(bxV); if (bxR == 0) bxR <- 1
yLoBox <- min(bxV) - 0.05 * bxR    # pad below the lowest point
yHiBox <- max(bxV) + 0.20 * bxR    # headroom above the tallest box for the asterisk + text
p_brduN <- ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.5, alpha = 0.9, shape = 1) +
  scale_x_discrete(limits = c("BrdU- Max Bias", "BrdU- Exp30", "BrdU- Exp90",
                              "BrdU- Exp180", "BrdU- Exp320", "BrdU- Random"), labels = function(x) sub("^BrdU[+/-]+ ", "", x)) +
  scale_fill_manual(values = c("grey", "grey", "grey", "grey", "grey", "grey")) +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "black")) +
  labs(y = expression(plain(paste("Peaks", " ", "(", mu, "m", ")")))) +
  coord_cartesian(ylim = c(yLoBox, yHiBox)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin = margin(14, 16, 22, 40, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 32, face = "plain"),
    axis.text.y = element_text(size = 32, face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 32, face = "plain"),
    legend.position = "none"
  ) +
  add_sig_stars(pN_list,
                data_list = list("BrdU- Max Bias" = xmaxN,
                                 "BrdU- Exp30"    = xmaxNExp30,
                                 "BrdU- Exp90"    = xmaxNExp90,
                                 "BrdU- Exp180"   = xmaxNExp180,
                                 "BrdU- Exp320"   = xmaxNExp320),
                tested_x = tested_xN,
                ylim = c(yLoBox, yHiBox))

print(p_brduN)
ggsave("PeakNgg.eps", plot = p_brduN, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("PeakNgg.pdf", plot = p_brduN, device = "pdf", dpi = 1200, width = 10, height = 10)

