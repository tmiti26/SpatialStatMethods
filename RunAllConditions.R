rm(list=ls())
library(sp)
library(rmapshaper)
library(ggplot2)
library(magrittr)
library(sf)
library(spatstat)
library(GET)          # global rank-envelope test (Myllymaki et al. 2017)
library(goftest)
library(readr)
library(kSamples)

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))

## ==================== DRIVER SETTINGS ====================================
## Runs every condition below. It is RESUMABLE: for each (cond, i) it skips
## the work if that variation's p-value file already exists, and each variation
## is wrapped in tryCatch so a failure is logged (analysis_log.txt) and the run
## continues. Re-running the script picks up exactly where it stopped.
condList <- c("R", "3D", "Exp30", "Exp60", "Exp90", "Exp120", "Exp180", "Exp320")
N <- 3              # specimen (quadrat) number; must match chosenN in the Python generators
nruns <- 99          # number of generated variations produced by the Python generators
nsimEnv <- 199        # random-labelling simulations per variation. MUST be >= 19 for a 5% global test
                     # (GET needs nsimEnv+1 functions with (nsimEnv+1)*alpha >= 1). Use 199/999 for real runs.
pwRank <- max(1, min(floor((nsimEnv - 1)/2), round((nsimEnv + 1) * 0.025)))   # pointwise rank; ~5% two-sided
rsfBw <- 15          # fixed RSF smoothing bandwidth (~1 cell diameter, um); same for observed and all nulls

## guard: GET's global envelope test needs enough functions to resolve alpha.
## With s = nsimEnv + 1 functions, the smallest achievable p is 1/s, so a 5% test
## requires s*0.05 >= 1  ->  nsimEnv >= 19. Stop early with a clear message otherwise.
if ((nsimEnv + 1) * 0.05 < 1) {
  stop(sprintf("nsimEnv = %d is too small for a 5%% global test: need nsimEnv >= 19 (got s = %d functions).",
               nsimEnv, nsimEnv + 1))
}
## ========================================================================

#mynumber <- paste(readLines("number.txt"), collapse=" ")
number <- 3#as.numeric(mynumber)
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

xmin1 <-0
xmax1 <-quadratSize
ymin1 <-0
ymax1 <-quadratSize

rPCF <-seq(from = 0, to = PCF_r, by = 0.5)
rJ <-seq(from = 0, to = J_r, by = 0.5)

## ---- NA handling shared by every summary function (fill sparse small-r NAs by
## up-to-2 non-NA neighbours per side; truncate the unsupported large-r tail) ----
myFillNeighbours <- function(v){
  bad <- which(!is.finite(v))
  for (k in bad){
    left  <- rev(v[seq_len(k-1)]); left  <- left[is.finite(left)]
    right <- v[seq.int(k+1, length(v))]; right <- right[is.finite(right)]
    nb <- c(head(left, 2), head(right, 2))
    if (length(nb) > 0) v[k] <- mean(nb)
  }
  v
}

myCleanCurves <- function(myObs, mySim){
  allMat <- rbind(myObs, mySim)
  finiteCol <- apply(allMat, 2, function(col) all(is.finite(col)))
  nc <- length(myObs)
  tailStart <- nc + 1
  for (cc in nc:1){ if (!finiteCol[cc]) tailStart <- cc else break }
  lastGood <- tailStart - 1
  keep <- seq_len(max(lastGood, 2))
  list(keep = keep,
       obs = myFillNeighbours(myObs[keep]),
       sim = t(apply(mySim[, keep, drop = FALSE], 1, myFillNeighbours)))
}

myGlobalFromEnv <- function(env){
  simDf <- as.data.frame(attr(env, "simfuns"))
  rr <- simDf[[1]]
  simMat <- t(as.matrix(simDf[, -1, drop = FALSE]))
  cl <- myCleanCurves(as.numeric(env$obs), simMat)
  cset <- create_curve_set(list(r = rr[cl$keep], obs = cl$obs, sim_m = t(cl$sim)))
  global_envelope_test(cset, type = "erl")
}

myGlobalBand <- function(myObs, mySim, myGrid){
  cl <- myCleanCurves(myObs, mySim)
  myObs <- cl$obs; mySim <- cl$sim; myGrid <- myGrid[cl$keep]
  myPwLo <- apply(mySim, 2, quantile, probs = 0.025, na.rm = TRUE)
  myPwHi <- apply(mySim, 2, quantile, probs = 0.975, na.rm = TRUE)
  myCset <- create_curve_set(list(r = myGrid, obs = myObs, sim_m = t(mySim)))
  myGe <- global_envelope_test(myCset, type = "erl")
  data.frame(r = myGrid, obs = myObs,
             pwLo = myPwLo, pwHi = myPwHi,
             globLo = myGe$lo, globHi = myGe$hi,
             central = myGe$central,
             pglobal = attr(myGe, "p"))
}

## ---- the analysis for ONE variation (cond, i): identical statistics to the
## original script; returns nothing, writes this variation's output files ------
myRunVariation <- function(i, cond){

  mydataPos <- read.csv(file = paste0("Pos",cond,"_",i,"_",N,".csv"))
  mydataNeg <- read.csv(file = paste0("Neg",cond,"_",i,"_",N,".csv"))
  mydataStr <- read.csv(file = paste0("Str",cond,"_",i,"_",N,".csv"))

  mypatternPos <- ppp(mydataPos[,1], mydataPos[,2], c(xmin1,xmax1), c(ymin1,ymax1))
  mypatternNeg <- ppp(mydataNeg[,1], mydataNeg[,2], c(xmin1,xmax1), c(ymin1,ymax1))
  mypatternStr <- ppp(mydataStr[,1], mydataStr[,2], c(xmin1,xmax1), c(ymin1,ymax1))

  if ((mypatternStr$n + mypatternPos$n + mypatternNeg$n) > 3000){

    myPvalsI <- data.frame()

    ## null used everywhere: fixed stroma + fixed cell positions, permuted BrdU labels
    PN <- superimpose(Pos = mypatternPos, Neg = mypatternNeg)
    mysim <- function(ppp1,ppp2){ superimpose(PosR = split(rlabel(ppp1))$Pos, Str = ppp2) }
    mysimN <- function(ppp1,ppp2){ superimpose(PosR = split(rlabel(ppp1))$Neg, Str = ppp2) }

    ## nearest-stroma distances
    nn <- nncross(mypatternPos, mypatternStr)$dist
    mm <- nncross(mypatternNeg, mypatternStr)$dist
    write.table(nn, file = paste0(i,"_myNNP",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(mm, file = paste0(i,"_myNNN",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")

    myAllD <- nncross(PN, mypatternStr)$dist
    Zrsf <- distmap(mypatternStr)
    myMaxD <- max(max(myAllD), max(Zrsf)) * 1.001
    Z <- Zrsf

    ## availability baseline = intensity of ALL tumour cells (fixed spatial sigma = rsfBw = 15 um)
    myBaseline <- density(unmark(PN), sigma = rsfBw)

    ## ---- RSF, FOUR variants (BrdU+), each with raw curve + pointwise/global band ----
    ##   v1 = original : no grid, floating bw, NO baseline
    ##   v2 = baseline : no grid, floating bw, baseline
    ##   v3 = grid     : fixed grid + fixed bw, NO baseline
    ##   v4 = grid+base: fixed grid + fixed bw + baseline
    ## helper: build the band for one variant. simFun(X) returns the rhohat for pattern X
    ## with that variant's settings. Sim curves are stacked onto the observed grid (or
    ## interpolated only if a floating grid returns a different length), then handed to
    ## myGlobalBand (which does the NA handling / tail truncation at the GET boundary).
    myRsfVariant <- function(rhoObs, simFun){
      dfo <- as.data.frame(rhoObs); zg <- dfo[[1]]; obs <- dfo$rho
      simMat <- matrix(NA_real_, nrow = nsimEnv, ncol = length(zg))
      for (bb in 1:nsimEnv){
        dfs <- as.data.frame(simFun(split(rlabel(PN))$Pos))
        sr <- dfs$rho
        if (length(sr) == length(zg)) simMat[bb, ] <- sr
        else simMat[bb, ] <- approx(x = dfs[[1]], y = sr, xout = zg, rule = 2)$y
      }
      myGlobalBand(obs, simMat, zg)
    }

    ## BrdU- raw curve (single, fixed-grid style; kept for reference)
    myRho <- rhohat(mypatternNeg, Z, from = 0, to = myMaxD, n = dens_r_length, bw = rsfBw)
    write.table(myRho, file = paste0(i,"_myRhoN",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")

    ## v1 : original (no grid, floating bw, no baseline)
    myRhoPv1 <- rhohat(mypatternPos, Z)
    myBandV1 <- myRsfVariant(myRhoPv1, function(X) rhohat(X, Z))
    write.table(myRhoPv1, file = paste0(i,"_myRhoPv1",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myBandV1, file = paste0(i,"_myRhoPbandv1",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")

    ## v2 : + baseline (no grid, floating bw)
    myRhoPv2 <- rhohat(mypatternPos, Z, baseline = myBaseline)
    myBandV2 <- myRsfVariant(myRhoPv2, function(X) rhohat(X, Z, baseline = myBaseline))
    write.table(myRhoPv2, file = paste0(i,"_myRhoPv2",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myBandV2, file = paste0(i,"_myRhoPbandv2",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")

    ## v3 : grid + fixed bw, no baseline
    myRhoPv3 <- rhohat(mypatternPos, Z, from = 0, to = myMaxD, n = dens_r_length, bw = rsfBw)
    myBandV3 <- myRsfVariant(myRhoPv3, function(X) rhohat(X, Z, from = 0, to = myMaxD, n = dens_r_length, bw = rsfBw))
    write.table(myRhoPv3, file = paste0(i,"_myRhoPv3",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myBandV3, file = paste0(i,"_myRhoPbandv3",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")

    ## v4 : grid + fixed bw + baseline
    myRhoPv4 <- rhohat(mypatternPos, Z, baseline = myBaseline, from = 0, to = myMaxD, n = dens_r_length, bw = rsfBw)
    myBandV4 <- myRsfVariant(myRhoPv4, function(X) rhohat(X, Z, baseline = myBaseline, from = 0, to = myMaxD, n = dens_r_length, bw = rsfBw))
    write.table(myRhoPv4, file = paste0(i,"_myRhoPv4",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myBandV4, file = paste0(i,"_myRhoPbandv4",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")

    myPvalsI <- rbind(myPvalsI,
                      data.frame(i = i, stat = "RSF_Pos_v1", p = myBandV1$pglobal[1]),
                      data.frame(i = i, stat = "RSF_Pos_v2", p = myBandV2$pglobal[1]),
                      data.frame(i = i, stat = "RSF_Pos_v3", p = myBandV3$pglobal[1]),
                      data.frame(i = i, stat = "RSF_Pos_v4", p = myBandV4$pglobal[1]))

    ## ---- nearest-stroma distance: ECDF pointwise/global envelope + KS p-value ----
    myDeval <- seq(0, myMaxD, length.out = J_r_length)
    myEcdf <- function(myD, myGrid) sapply(myGrid, function(t) mean(myD <= t))
    myFobs <- myEcdf(nn, myDeval)
    myFsim <- matrix(NA_real_, nrow = nsimEnv, ncol = length(myDeval))
    for (bb in 1:nsimEnv){
      myNullPos <- split(rlabel(PN))$Pos
      myFsim[bb, ] <- myEcdf(nncross(myNullPos, mypatternStr)$dist, myDeval)
    }
    myNNband <- myGlobalBand(myFobs, myFsim, myDeval)
    write.table(myNNband, file = paste0(i,"_myNNPband",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")

    myDall <- myAllD
    myKSobs <- as.numeric(ks.test(nn, myDall)$statistic)
    myKSsim <- numeric(nsimEnv)
    for (bb in 1:nsimEnv){
      myKSsim[bb] <- as.numeric(ks.test(nncross(split(rlabel(PN))$Pos, mypatternStr)$dist, myDall)$statistic)
    }
    myKSp <- (1 + sum(myKSsim >= myKSobs)) / (1 + nsimEnv)
    myPvalsI <- rbind(myPvalsI,
                      data.frame(i = i, stat = "NN_Pos_global", p = myNNband$pglobal[1]),
                      data.frame(i = i, stat = "NN_Pos_KS",     p = myKSp))

    ## ---- cross summary functions: pointwise (envelope) + global (rank envelope) ----
    XY  <- superimpose(Pos=mypatternPos,Str=mypatternStr)
    XYN <- superimpose(Neg=mypatternNeg,Str=mypatternStr)

    ## RDF (pcfcross)
   # myenv  <- envelope(XYN, fun = pcfcross, correction=c("Ripley"), r = rPCF, nsim=nsimEnv, nrank=pwRank, savefuns=TRUE, simulate=expression(mysimN(PN,mypatternStr)))
    myenvp <- envelope(XY,  fun = pcfcross, correction=c("Ripley"), r = rPCF, nsim=nsimEnv, nrank=pwRank, savefuns=TRUE, simulate=expression(mysim(PN,mypatternStr)))
    write.table(myenvp, file = paste0(i,"_myenvRDFP",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    #write.table(myenv,  file = paste0(i,"_myenvRDFN",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    myenvpGlob <- myGlobalFromEnv(myenvp); #myenvGlob <- myGlobalFromEnv(myenv)
    myRDFPg <- as.data.frame(myenvpGlob); myRDFPg$pglobal <- attr(myenvpGlob,"p")
    #myRDFNg <- as.data.frame(myenvGlob);  myRDFNg$pglobal <- attr(myenvGlob,"p")
    write.table(myRDFPg, file = paste0(i,"_myenvRDFPglob",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    #write.table(myRDFNg, file = paste0(i,"_myenvRDFNglob",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    myPvalsI <- rbind(myPvalsI,
                      data.frame(i = i, stat = "RDF_Pos", p = attr(myenvpGlob,"p")))
                      #data.frame(i = i, stat = "RDF_Neg", p = attr(myenvGlob,"p")))


    ## Jcross
    #myenvJ  <- envelope(XYN, fun = Jcross, correction=c("rs"), r = rJ, nsim=nsimEnv, nrank=pwRank, savefuns=TRUE, simulate=expression(mysimN(PN,mypatternStr)))
    myenvpJ <- envelope(XY,  fun = Jcross, correction=c("rs"), r = rJ, nsim=nsimEnv, nrank=pwRank, savefuns=TRUE, simulate=expression(mysim(PN,mypatternStr)))
    write.table(myenvpJ, file = paste0(i,"_myenvJcrossP",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    #write.table(myenvJ,  file = paste0(i,"_myenvJcrossN",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    myenvpJGlob <- myGlobalFromEnv(myenvpJ); #myenvJGlob <- myGlobalFromEnv(myenvJ)
    myJPg <- as.data.frame(myenvpJGlob); myJPg$pglobal <- attr(myenvpJGlob,"p")
    #myJNg <- as.data.frame(myenvJGlob);  myJNg$pglobal <- attr(myenvJGlob,"p")
    write.table(myJPg, file = paste0(i,"_myenvJcrossPglob",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    #write.table(myJNg, file = paste0(i,"_myenvJcrossNglob",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
    myPvalsI <- rbind(myPvalsI,
                      data.frame(i = i, stat = "Jcross_Pos", p = attr(myenvpJGlob,"p")))
                      #data.frame(i = i, stat = "Jcross_Neg", p = attr(myenvJGlob,"p")))


    ## LAST write = the per-variation p-value file. Its existence is the "done" marker
    ## used for resuming, so it is written only after everything else has succeeded.
    write.table(myPvalsI, file = paste0(i,"_pvals_",cond,"_",N,".csv"),row.names=FALSE,col.names=TRUE,sep=",")
  }
}

## ==================== DRIVER LOOP (resumable) ============================
for (cond in condList){
  message("==== condition: ", cond, " ====")
  for (i in 1:nruns){

    doneMarker <- paste0(i,"_pvals_",cond,"_",N,".csv")
    if (file.exists(doneMarker)) next          # RESUME: this variation already finished

    if (!file.exists(paste0("Pos",cond,"_",i,"_",N,".csv"))){
      cat(sprintf("SKIP (input not generated) cond=%s i=%s\n", cond, i),
          file = "analysis_log.txt", append = TRUE)
      next
    }

    tryCatch(
      {
        myRunVariation(i, cond)
        message("  cond=", cond, " i=", i, " done")
      },
      error = function(e){
        cat(sprintf("ERROR cond=%s i=%s : %s\n", cond, i, conditionMessage(e)),
            file = "analysis_log.txt", append = TRUE)
        message("  cond=", cond, " i=", i, " ERROR (logged) -- continuing")
      }
    )
  }

  ## aggregate this condition's per-variation p-value files into pvalues_<cond>_<N>.csv
  pvFiles <- list.files(pattern = paste0("^[0-9]+_pvals_", cond, "_", N, "\\.csv$"))
  if (length(pvFiles) > 0){
    myPvalsAll <- do.call(rbind, lapply(pvFiles, read.csv))
    myPvalsAll <- myPvalsAll[order(myPvalsAll$i), ]
    write.table(myPvalsAll, file = paste0("pvalues_",cond,"_",N,".csv"),
                row.names=FALSE, col.names=TRUE, sep=",")
    message("  wrote pvalues_", cond, "_", N, ".csv  (", length(pvFiles), " variations)")
  }
}
message("ALL DONE")
