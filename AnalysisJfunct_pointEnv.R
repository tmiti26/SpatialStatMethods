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
  assign(paste("JcrossP", number, sep = ""),  (read.csv(file = paste(i,"_myenvJcrossP3D_", N, ".csv",sep=""))))
  #assign(paste("JcrossN", number, sep = ""),  (read.csv(file = paste(i,"_myenvJcrossN3D_", N, ".csv",sep=""))))
  
  assign(paste("JcrossPRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPR_", N, ".csv",sep=""))))
  #assign(paste("JcrossNRand", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNR_", N, ".csv",sep=""))))
  
  assign(paste("JcrossPExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp30_", N, ".csv",sep=""))))
  #assign(paste("JcrossNExp30_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp30_", N, ".csv",sep=""))))
  
  assign(paste("JcrossPExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp60_", N, ".csv",sep=""))))
  #assign(paste("JcrossNExp60_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp60_", N, ".csv",sep=""))))
  
  assign(paste("JcrossPExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp90_", N, ".csv",sep=""))))
  #assign(paste("JcrossNExp90_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp90_", N, ".csv",sep=""))))
  
  #assign(paste("JcrossPExp120_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp120_", N, ".csv",sep=""))))
  #assign(paste("JcrossNExp120_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp120_", N, ".csv",sep=""))))
  
  # assign(paste("JcrossPExp150_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp150_", N, ".csv",sep=""))))
  # assign(paste("JcrossNExp150_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp150_", N, ".csv",sep=""))))
  
  assign(paste("JcrossPExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp180_", N, ".csv",sep=""))))
  #assign(paste("JcrossNExp180_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp180_", N, ".csv",sep=""))))
  
  # assign(paste("JcrossPExp240_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp240_", N, ".csv",sep=""))))
  # assign(paste("JcrossNExp240_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp240_", N, ".csv",sep=""))))
  
  assign(paste("JcrossPExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossPExp320_", N, ".csv",sep=""))))
  #assign(paste("JcrossNExp320_", number, sep = ""),  (read.csv(file = paste(number,"_myenvJcrossNExp320_", N, ".csv",sep=""))))
  
}

#*********************************** J function ***********************
x <- JcrossP1$r

for (i in 1:numbTissue){ 
  pop.name       = paste0("JcrossP",i)
  #pop.nameN      = paste0("JcrossN",i)
  pop.nameR      = paste0("JcrossPRand",i)
  #pop.nameRN     = paste0("JcrossNRand",i)
  pop.namePExp30  = paste0("JcrossPExp30_",i)
  #pop.nameNExp30  = paste0("JcrossNExp30_",i)
  pop.namePExp60  = paste0("JcrossPExp60_",i)
  #pop.nameNExp60  = paste0("JcrossNExp60_",i)
  pop.namePExp90  = paste0("JcrossPExp90_",i)
  #pop.nameNExp90  = paste0("JcrossNExp90_",i)
  #pop.namePExp120 = paste0("JcrossPExp120_",i)
  #pop.nameNExp120 = paste0("JcrossNExp120_",i)
  pop.namePExp180 = paste0("JcrossPExp180_",i)
  #pop.nameNExp180 = paste0("JcrossNExp180_",i)
  pop.namePExp320 = paste0("JcrossPExp320_",i)
  #pop.nameNExp320 = paste0("JcrossNExp320_",i)
  
  assign(paste0("Jdiff",      i), (eval(parse(text = pop.name))$obs)  / (eval(parse(text = pop.name))$mmean))
  #assign(paste0("JdiffN",     i), (eval(parse(text = pop.nameN))$obs) / (eval(parse(text = pop.nameN))$mmean))
  assign(paste0("JdiffRand",  i), (eval(parse(text = pop.nameR))$obs) / (eval(parse(text = pop.nameR))$mmean))
  #assign(paste0("JdiffNRand", i), (eval(parse(text = pop.nameRN))$obs)/ (eval(parse(text = pop.nameRN))$mmean))
  assign(paste0("JdiffL",      i), (eval(parse(text = pop.name))$lo)  / (eval(parse(text = pop.name))$mmean))
  assign(paste0("JdiffH",      i), (eval(parse(text = pop.name))$hi)  / (eval(parse(text = pop.name))$mmean))
  assign(paste0("JdiffRandL",  i), (eval(parse(text = pop.nameR))$lo) / (eval(parse(text = pop.nameR))$mmean))
  assign(paste0("JdiffRandH",  i), (eval(parse(text = pop.nameR))$hi) / (eval(parse(text = pop.nameR))$mmean))
  assign(paste0("JdiffPExp30_",  i), (eval(parse(text = pop.namePExp30))$obs)  / (eval(parse(text = pop.namePExp30))$mmean))
  #assign(paste0("JdiffNExp30_",  i), (eval(parse(text = pop.nameNExp30))$obs)  / (eval(parse(text = pop.nameNExp30))$mmean))
  assign(paste0("JdiffPExp60_",  i), (eval(parse(text = pop.namePExp60))$obs)  / (eval(parse(text = pop.namePExp60))$mmean))
  #assign(paste0("JdiffNExp60_",  i), (eval(parse(text = pop.nameNExp60))$obs)  / (eval(parse(text = pop.nameNExp60))$mmean))
  assign(paste0("JdiffPExp90_",  i), (eval(parse(text = pop.namePExp90))$obs)  / (eval(parse(text = pop.namePExp90))$mmean))
  #assign(paste0("JdiffNExp90_",  i), (eval(parse(text = pop.nameNExp90))$obs)  / (eval(parse(text = pop.nameNExp90))$mmean))
  #assign(paste0("JdiffPExp120_", i), (eval(parse(text = pop.namePExp120))$obs) / (eval(parse(text = pop.namePExp120))$mmean))
  #assign(paste0("JdiffNExp120_", i), (eval(parse(text = pop.nameNExp120))$obs) / (eval(parse(text = pop.nameNExp120))$mmean))
  assign(paste0("JdiffPExp180_", i), (eval(parse(text = pop.namePExp180))$obs) / (eval(parse(text = pop.namePExp180))$mmean))
  #assign(paste0("JdiffNExp180_", i), (eval(parse(text = pop.nameNExp180))$obs) / (eval(parse(text = pop.nameNExp180))$mmean))
  assign(paste0("JdiffPExp320_", i), (eval(parse(text = pop.namePExp320))$obs) / (eval(parse(text = pop.namePExp320))$mmean))
  #assign(paste0("JdiffNExp320_", i), (eval(parse(text = pop.nameNExp320))$obs) / (eval(parse(text = pop.nameNExp320))$mmean))
}

JPtot     <- rep(0, length(Jdiff1))
JNtot     <- rep(0, length(Jdiff1))
JPRandtot <- rep(0, length(Jdiff1))
JNRandtot <- rep(0, length(Jdiff1))
JPtotL     <- rep(0, length(Jdiff1))
JPtotH     <- rep(0, length(Jdiff1))
JPRandtotL <- rep(0, length(Jdiff1))
JPRandtotH <- rep(0, length(Jdiff1))
JPExp30   <- rep(0, length(Jdiff1))
JNExp30   <- rep(0, length(Jdiff1))
JPExp60   <- rep(0, length(Jdiff1))
JNExp60   <- rep(0, length(Jdiff1))
JPExp90   <- rep(0, length(Jdiff1))
JNExp90   <- rep(0, length(Jdiff1))
JPExp120  <- rep(0, length(Jdiff1))
JNExp120  <- rep(0, length(Jdiff1))
JPExp180  <- rep(0, length(Jdiff1))
JNExp180  <- rep(0, length(Jdiff1))
JPExp320  <- rep(0, length(Jdiff1))
JNExp320  <- rep(0, length(Jdiff1))

for (i in 1:numbTissue){ 
  pop.name        = paste0("Jdiff",i)
  #pop.nameN       = paste0("JdiffN",i)
  pop.nameR       = paste0("JdiffRand",i)
  #pop.nameRN      = paste0("JdiffNRand",i)
  pop.namePExp30  = paste0("JdiffPExp30_",i)
  #pop.nameNExp30  = paste0("JdiffNExp30_",i)
  pop.namePExp60  = paste0("JdiffPExp60_",i)
  #pop.nameNExp60  = paste0("JdiffNExp60_",i)
  pop.namePExp90  = paste0("JdiffPExp90_",i)
  #pop.nameNExp90  = paste0("JdiffNExp90_",i)
  #pop.namePExp120 = paste0("JdiffPExp120_",i)
  #pop.nameNExp120 = paste0("JdiffNExp120_",i)
  pop.namePExp180 = paste0("JdiffPExp180_",i)
  #pop.nameNExp180 = paste0("JdiffNExp180_",i)
  pop.namePExp320 = paste0("JdiffPExp320_",i)
  #pop.nameNExp320 = paste0("JdiffNExp320_",i)
  
  JPtot     <- JPtot     + eval(parse(text = pop.name))       / numbTissue
  #JNtot     <- JNtot     + eval(parse(text = pop.nameN))      / numbTissue
  JPRandtot <- JPRandtot + eval(parse(text = pop.nameR))      / numbTissue
  #JNRandtot <- JNRandtot + eval(parse(text = pop.nameRN))     / numbTissue
  JPtotL     <- JPtotL     + eval(parse(text = paste0("JdiffL",i)))     / numbTissue
  JPtotH     <- JPtotH     + eval(parse(text = paste0("JdiffH",i)))     / numbTissue
  JPRandtotL <- JPRandtotL + eval(parse(text = paste0("JdiffRandL",i))) / numbTissue
  JPRandtotH <- JPRandtotH + eval(parse(text = paste0("JdiffRandH",i))) / numbTissue
  JPExp30   <- JPExp30   + eval(parse(text = pop.namePExp30)) / numbTissue
  #JNExp30   <- JNExp30   + eval(parse(text = pop.nameNExp30)) / numbTissue
  JPExp60   <- JPExp60   + eval(parse(text = pop.namePExp60)) / numbTissue
  #JNExp60   <- JNExp60   + eval(parse(text = pop.nameNExp60)) / numbTissue
  JPExp90   <- JPExp90   + eval(parse(text = pop.namePExp90)) / numbTissue
  #JNExp90   <- JNExp90   + eval(parse(text = pop.nameNExp90)) / numbTissue
  #JPExp120  <- JPExp120  + eval(parse(text = pop.namePExp120))/ numbTissue
  #JNExp120  <- JNExp120  + eval(parse(text = pop.nameNExp120))/ numbTissue
  JPExp180  <- JPExp180  + eval(parse(text = pop.namePExp180))/ numbTissue
  #JNExp180  <- JNExp180  + eval(parse(text = pop.nameNExp180))/ numbTissue
  JPExp320  <- JPExp320  + eval(parse(text = pop.namePExp320))/ numbTissue
  #JNExp320  <- JNExp320  + eval(parse(text = pop.nameNExp320))/ numbTissue
}


########## J values at r=120
xmaxJ      <- c()
xmaxJRand  <- c()
xmaxJN     <- c()
xmaxJNRand <- c()
xmaxJPExp30  <- c()
xmaxJNExp30  <- c()
xmaxJPExp60  <- c()
xmaxJNExp60  <- c()
xmaxJPExp90  <- c()
xmaxJNExp90  <- c()
xmaxJPExp120 <- c()
xmaxJNExp120 <- c()
xmaxJPExp180 <- c()
xmaxJNExp180 <- c()
xmaxJPExp320 <- c()
xmaxJNExp320 <- c()

for (i in 1:numbTissue){ 
  pop.name        = paste0("Jdiff",i)
  pop.nameR       = paste0("JdiffRand",i)
  #pop.nameN       = paste0("JdiffN",i)
  #pop.nameRN      = paste0("JdiffNRand",i)
  pop.namePExp30  = paste0("JdiffPExp30_",i)
  #pop.nameNExp30  = paste0("JdiffNExp30_",i)
  pop.namePExp60  = paste0("JdiffPExp60_",i)
  #pop.nameNExp60  = paste0("JdiffNExp60_",i)
  pop.namePExp90  = paste0("JdiffPExp90_",i)
  #pop.nameNExp90  = paste0("JdiffNExp90_",i)
  #pop.namePExp120 = paste0("JdiffPExp120_",i)
  #pop.nameNExp120 = paste0("JdiffNExp120_",i)
  pop.namePExp180 = paste0("JdiffPExp180_",i)
  #pop.nameNExp180 = paste0("JdiffNExp180_",i)
  pop.namePExp320 = paste0("JdiffPExp320_",i)
  #pop.nameNExp320 = paste0("JdiffNExp320_",i)
  
  xmaxJ       <- c(xmaxJ,       eval(parse(text = pop.name))[120])
  xmaxJRand   <- c(xmaxJRand,   eval(parse(text = pop.nameR))[120])
  #xmaxJN      <- c(xmaxJN,      eval(parse(text = pop.nameN))[120])
  #xmaxJNRand  <- c(xmaxJNRand,  eval(parse(text = pop.nameRN))[120])
  xmaxJPExp30  <- c(xmaxJPExp30,  eval(parse(text = pop.namePExp30))[120])
  #xmaxJNExp30  <- c(xmaxJNExp30,  eval(parse(text = pop.nameNExp30))[120])
  xmaxJPExp60  <- c(xmaxJPExp60,  eval(parse(text = pop.namePExp60))[120])
  #xmaxJNExp60  <- c(xmaxJNExp60,  eval(parse(text = pop.nameNExp60))[120])
  xmaxJPExp90  <- c(xmaxJPExp90,  eval(parse(text = pop.namePExp90))[120])
  #xmaxJNExp90  <- c(xmaxJNExp90,  eval(parse(text = pop.nameNExp90))[120])
  #xmaxJPExp120 <- c(xmaxJPExp120, eval(parse(text = pop.namePExp120))[120])
  #xmaxJNExp120 <- c(xmaxJNExp120, eval(parse(text = pop.nameNExp120))[120])
  xmaxJPExp180 <- c(xmaxJPExp180, eval(parse(text = pop.namePExp180))[120])
  #xmaxJNExp180 <- c(xmaxJNExp180, eval(parse(text = pop.nameNExp180))[120])
  xmaxJPExp320 <- c(xmaxJPExp320, eval(parse(text = pop.namePExp320))[120])
  #xmaxJNExp320 <- c(xmaxJNExp320, eval(parse(text = pop.nameNExp320))[120])
}


####### significance helper functions
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

####### BrdU+ boxplot (Jmin)
WCPJ     <- t.test(xmaxJ,     xmaxJRand,    paired = TRUE, alternative = "two.sided")
WCPJR30  <- t.test(xmaxJRand, xmaxJPExp30,  paired = TRUE, alternative = "two.sided")
WCPJR90  <- t.test(xmaxJRand, xmaxJPExp90,  paired = TRUE, alternative = "two.sided")
WCPJR180 <- t.test(xmaxJRand, xmaxJPExp180, paired = TRUE, alternative = "two.sided")
WCPJR320 <- t.test(xmaxJRand, xmaxJPExp320, paired = TRUE, alternative = "two.sided")

data <- data.frame(
  Category = factor(rep(c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                          "BrdU+ Exp180",  "BrdU+ Exp320", "BrdU+ Random"), each = numbTissue),
                    levels = c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                               "BrdU+ Exp180",  "BrdU+ Exp320", "BrdU+ Random")),
  Value = c("BrdU+ Max Bias" = xmaxJ,
            "BrdU+ Exp30"   = xmaxJPExp30,
            "BrdU+ Exp90"   = xmaxJPExp90,
            "BrdU+ Exp180"  = xmaxJPExp180,
            "BrdU+ Exp320"  = xmaxJPExp320,
            "BrdU+ Random"  = xmaxJRand)
)

# x positions: 1=Max Bias, 2=Exp30, 3=Exp90, 4=Exp180, 5=Exp320, 6=Random
pP_list <- list(
  "BrdU+ Max Bias" = list(p.value = WCPJ$p.value),
  "BrdU+ Exp30"    = list(p.value = WCPJR30$p.value),
  "BrdU+ Exp90"    = list(p.value = WCPJR90$p.value),
  "BrdU+ Exp180"   = list(p.value = WCPJR180$p.value),
  "BrdU+ Exp320"   = list(p.value = WCPJR320$p.value)
)
tested_xP <- c("BrdU+ Max Bias" = 1, "BrdU+ Exp30" = 2, "BrdU+ Exp90" = 3,
               "BrdU+ Exp180" = 4, "BrdU+ Exp320" = 5)

bxV <- data$Value[is.finite(data$Value)]
bxR <- max(bxV) - min(bxV); if (bxR == 0) bxR <- 1
yLoBox <- min(bxV) - 0.05 * bxR    # pad below the lowest point
yHiBox <- max(bxV) + 0.20 * bxR    # headroom above the tallest box for the asterisk + text
p_jP <- ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Category), width = 0.2, size = 1.5, alpha = 0.9, shape = 1) +
  scale_x_discrete(limits = c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                              "BrdU+ Exp180",  "BrdU+ Exp320", "BrdU+ Random"), labels = function(x) sub("^BrdU[+/-]+ ", "", x)) +
  scale_fill_manual(values = c("grey", "grey", "grey", "grey", "grey", "grey")) +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "black")) +
  labs(y = expression(plain(paste("Jmin")))) +
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
                data_list = list("BrdU+ Max Bias" = xmaxJ,
                                 "BrdU+ Exp30"    = xmaxJPExp30,
                                 "BrdU+ Exp90"    = xmaxJPExp90,
                                 "BrdU+ Exp180"   = xmaxJPExp180,
                                 "BrdU+ Exp320"   = xmaxJPExp320),
                tested_x = tested_xP,
                ylim = c(yLoBox, yHiBox))

print(p_jP)
ggsave("JXmaxggP.eps", plot = p_jP, device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("JXmaxggP.pdf", plot = p_jP, device = "pdf",      dpi = 1200, width = 10, height = 10)

######### J function line plot (BrdU+)
df <- data.frame(
  x = rep(JcrossP1$r, 6),
  y = c(JPtot, JPExp30, JPExp90, JPExp180, JPExp320, JPRandtot),
  group = rep(c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                "BrdU+ Exp180",  "BrdU+ Exp320", "BrdU+ Random"),
              each = length(JPtot))
)

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
  "BrdU+ Random"   = "solid")

## dynamic y-range: fit all condition curves + the null envelope within the x-window (+5% pad)
xLo <- 0; xHi <- myMergeXmax()
inWin <- df$x >= xLo & df$x <= xHi
inR   <- JcrossP1$r >= xLo & JcrossP1$r <= xHi
yAll  <- c(df$y[inWin], JPRandtotL[inR], JPRandtotH[inR]); yAll <- yAll[is.finite(yAll)]
yPad  <- 0.05 * (max(yAll) - min(yAll)); yLo <- min(yAll) - yPad; yHi <- max(yAll) + yPad

ggplot(df, aes(x = x, y = y, color = group, linetype = group)) +
  # random-labelling null envelope (identical for every condition -> shown once, in grey)
  geom_ribbon(data = data.frame(x = JcrossP1$r, ymin = JPRandtotL, ymax = JPRandtotH),
              aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
              fill = "grey60", alpha = 0.25) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = line_mapping) +
  labs(
    x = expression(plain(paste("Distance from Stroma", " ", "(", mu, "m", ")"))),
    y = expression(plain(paste("J Function")))
  ) +
  coord_cartesian(ylim = c(yLo, yHi), xlim = c(xLo, xHi)) +
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

ggsave("JCombAll.eps", device = cairo_ps, dpi = 1200, width = 10, height = 10, family = "Arial")
ggsave("JCombAll.pdf", device = "pdf",      dpi = 1200, width = 10, height = 10)
