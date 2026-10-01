
rm(list = ls())
library(ggplot2)
library(grid)

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

## ================= settings =================
numbTissue <- 99
N <- 12

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

condList  <- c("3D", "Exp30", "Exp90", "Exp180", "Exp320", "R")
condLabel <- c("3D"="BrdU+ Max Bias", "Exp30"="BrdU+ Exp30", "Exp90"="BrdU+ Exp90",
               "Exp180"="BrdU+ Exp180", "Exp320"="BrdU+ Exp320", "R"="BrdU+ Random")
variantTitle <- c("v1: no grid, floating bw, no baseline",
                  "v2: no grid, floating bw, with baseline",
                  "v3: grid, fixed 15um bw, no baseline",
                  "v4: grid, fixed 15um bw, with baseline")

## colours (match the other figures)
pink_gradient <- c("#7B2D8B", "#B5006A", "#E8003A", "#FF7800", "#FF99CC")
names(pink_gradient) <- c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                          "BrdU+ Exp180", "BrdU+ Exp320")
color_mapping <- c(pink_gradient, "BrdU+ Random" = "black")

## ---- significance-star helpers (shared with the other scripts) ----
myAsterics <- function(myPval){
  if (myPval >= 0.05){return("ns")}
  else if((myPval <= 0.05) & (myPval > 0.005)){return("*")}
  else if((myPval <= 0.005) & (myPval > 0.0005)){return("**")}
  else {return("***")}
}
add_sig_stars <- function(p_list, data_list, tested_x, ylim, frac = 0.05) {
  offset <- frac * (ylim[2] - ylim[1])
  layers <- list()
  for (nm in names(p_list)) {
    pval <- p_list[[nm]]$p.value
    if (pval < 0.05) {
      y_star <- max(data_list[[nm]], na.rm = TRUE) + offset
      layers <- c(layers, list(
        annotate("text", x = tested_x[[nm]], y = y_star,
                 label = myAsterics(pval), size = 8, fontface = "bold")))
    }
  }
  layers
}

## ---- pull the populated legend grob out of a built ggplot (version-robust) ----
myGetLegend <- function(p){
  g <- ggplotGrob(p)
  nm <- vapply(g$grobs, function(x) x$name, character(1))
  idx <- which(grepl("guide-box", nm))
  if (!length(idx)) return(NULL)
  ng <- vapply(idx, function(k){ gb <- g$grobs[[k]]
    if (is.null(gb$grobs)) 0L else length(gb$grobs) }, integer(1))
  g$grobs[[ idx[ which.max(ng) ] ]]
}

## ---- write a grob (e.g. the legend) to its own file; device from the extension ----
mySaveGrob <- function(grob, file, width, height){
  ext <- tolower(tools::file_ext(file))
  if      (ext %in% c("eps","ps"))   cairo_ps(file, width = width, height = height)
  else if (ext == "png")             png(file, width = width, height = height, units = "in", res = 1200)
  else if (ext %in% c("tif","tiff")) tiff(file, width = width, height = height, units = "in", res = 1200)
  else                               pdf(file, width = width, height = height)
  on.exit(dev.off())
  grid.newpage(); grid.draw(grob)
}

## ---- read one (variant, condition) band file for variation i ----
myBandFile <- function(v, cond, i) read.csv(paste0(i, "_myRhoPbandv", v, cond, "_", N, ".csv"))

## ---- average the RSF band (obs + pointwise lo/hi) across variations ----
myReadMeanBand <- function(v, cond){
  myObs <- NULL; myLo <- NULL; myHi <- NULL; myR <- NULL
  for (i in 1:numbTissue){
    myD <- myBandFile(v, cond, i)
    if (is.null(myObs)){ myObs <- rep(0, nrow(myD)); myLo <- myObs; myHi <- myObs; myR <- myD$r }
    myObs <- myObs + myD$obs  / numbTissue
    myLo  <- myLo  + myD$pwLo / numbTissue
    myHi  <- myHi  + myD$pwHi / numbTissue
  }
  data.frame(r = myR, obs = myObs, lo = myLo, hi = myHi)
}

## ---- RSF index = integral of |rho/rhoBar - 1| dz for one variation (trapezoidal) ----
myRsfIndexOne <- function(rho, zz){
  ok <- is.finite(rho) & is.finite(zz); rho <- rho[ok]; zz <- zz[ok]
  if (length(rho) < 3) return(NA_real_)
  rhoBar <- mean(rho)
  val <- abs(rho / rhoBar - 1)
  dz  <- diff(zz)
  sum((val[-1] + val[-length(val)]) / 2 * dz)   # unitless
}
## ---- RSF index for every variation of a (variant, condition) ----
myRsfIndexVec <- function(v, cond){
  vapply(1:numbTissue, function(i){
    d <- myBandFile(v, cond, i); myRsfIndexOne(d$obs, d$r)
  }, numeric(1))
}

## ================= per-variant: curve figure + index boxplot =================
for (v in 1:4){
  
  ## ----- (A) RSF curves + random null band -----
  myCurves <- lapply(condList, function(cc) myReadMeanBand(v, cc)); names(myCurves) <- condList
  myGrid <- myCurves[["R"]]$r
  df <- do.call(rbind, lapply(condList, function(cc)
    data.frame(x = myCurves[[cc]]$r, y = myCurves[[cc]]$obs, group = condLabel[[cc]])))
  df$group <- factor(df$group, levels = unname(condLabel[condList]))
  myNullLo <- myCurves[["R"]]$lo; myNullHi <- myCurves[["R"]]$hi
  
  xLo <- 0; xHi <- myMergeXmax()
  inX <- df$x >= xLo & df$x <= xHi; inR <- myGrid >= xLo & myGrid <= xHi
  yAll <- c(df$y[inX], myNullLo[inR], myNullHi[inR]); yAll <- yAll[is.finite(yAll)]
  yPad <- 0.05 * (max(yAll) - min(yAll)); yLo <- min(yAll) - yPad; yHi <- max(yAll) + yPad
  
  pCurve <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_ribbon(data = data.frame(x = myGrid, ymin = myNullLo, ymax = myNullHi),
                aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
                fill = "grey60", alpha = 0.25) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_mapping,
                       labels = function(x) sub("^BrdU[+/-]+ ", "", x)) +
    labs(x = expression(plain(paste("Distance from Stroma", " ", "(", mu, "m", ")"))),
         y = expression(plain(paste("Resource Selection Function")))) +
    coord_cartesian(ylim = c(yLo, yHi), xlim = c(xLo, xHi)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 1),
          axis.line.y = element_line(color = "black", size = 1),
          axis.ticks.x = element_line(color = "black", size = 1),
          axis.ticks.y = element_line(color = "black", size = 1),
          axis.ticks.length = unit(0.25, "cm"),
          text = element_text(size = 32),                       # tick labels 32, horizontal
          axis.title = element_text(face = "plain", size = 32),  # BOTH x and y titles shown
          plot.margin = margin(14, 16, 22, 40, "pt"),
          legend.title = element_blank(),
          legend.text = element_text(size = 26, face = "plain"),
          legend.position = "right")

  ## ---- (a) main RSF curve WITHOUT the legend (10x10, pdf + eps) ----
  pCurveMain <- pCurve + theme(legend.position = "none")
  ggsave(paste0("RSF_v", v, ".pdf"), pCurveMain, width = 10, height = 10, dpi = 1200)
  ggsave(paste0("RSF_v", v, ".eps"), pCurveMain, device = cairo_ps,
         width = 10, height = 10, dpi = 1200, family = "Arial")
  print(pCurveMain)

  ## ---- (b) legend as its own image (pdf + eps), auto-sized ----
  leg <- myGetLegend(pCurve)
  if (!is.null(leg)){
    lw <- tryCatch(convertWidth(sum(leg$widths),   "in", valueOnly = TRUE), error = function(e) NA)
    lh <- tryCatch(convertHeight(sum(leg$heights), "in", valueOnly = TRUE), error = function(e) NA)
    if (!is.finite(lw) || lw <= 0) lw <- 6
    if (!is.finite(lh) || lh <= 0) lh <- 6
    mySaveGrob(leg, paste0("RSF_v", v, "_legend.pdf"), lw + 0.2, lh + 0.2)
    mySaveGrob(leg, paste0("RSF_v", v, "_legend.eps"), lw + 0.2, lh + 0.2)
    message("wrote RSF_v", v, "_legend.pdf and .eps")
  }
  
  ## ----- (B) RSF index boxplot across conditions -----
  idxList <- lapply(condList, function(cc) myRsfIndexVec(v, cc)); names(idxList) <- condList
  dataI <- data.frame(
    Category = factor(rep(unname(condLabel[condList]), each = numbTissue),
                      levels = unname(condLabel[condList])),
    Value    = unlist(idxList, use.names = FALSE))
  
  ## test each condition's index vs Random (paired across variations)
  testConds <- setdiff(condList, "R")
  pI_list <- setNames(lapply(testConds, function(cc)
    list(p.value = t.test(idxList[[cc]], idxList[["R"]], paired = TRUE)$p.value)),
    unname(condLabel[testConds]))
  tested_xI <- setNames(match(unname(condLabel[testConds]), levels(dataI$Category)),
                        unname(condLabel[testConds]))
  dataI_list <- setNames(lapply(testConds, function(cc) idxList[[cc]]), unname(condLabel[testConds]))
  
  bxV <- dataI$Value[is.finite(dataI$Value)]
  bxR <- max(bxV) - min(bxV); if (bxR == 0) bxR <- 1
  yLoBox <- min(bxV) - 0.05 * bxR
  yHiBox <- max(bxV) + 0.20 * bxR
  
  pIdx <- ggplot(dataI, aes(x = Category, y = Value, fill = Category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = Category), width = 0.2, size = 1.5, alpha = 0.9, shape = 1) +
    scale_x_discrete(limits = unname(condLabel[condList]),
                     labels = function(x) sub("^BrdU[+/-]+ ", "", x)) +
    scale_fill_manual(values = rep("grey", length(condList))) +
    scale_color_manual(values = rep("black", length(condList))) +
    labs(y = expression(plain(paste("RSF index")))) +
    coord_cartesian(ylim = c(yLoBox, yHiBox)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 1),
          axis.line.y = element_line(color = "black", size = 1),
          axis.ticks.x = element_line(color = "black", size = 1),
          axis.ticks.y = element_line(color = "black", size = 1),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 32, face = "plain"),  # boxplot: 45 deg
          axis.text.y = element_text(size = 32, face = "plain"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 32, face = "plain"),
          plot.margin = margin(14, 16, 22, 40, "pt"),
          legend.position = "none") +
    add_sig_stars(pI_list, data_list = dataI_list, tested_x = tested_xI,
                  ylim = c(yLoBox, yHiBox))
  ggsave(paste0("RSFindex_v", v, ".pdf"), pIdx, width = 10, height = 10, dpi = 1200)
  ggsave(paste0("RSFindex_v", v, ".eps"), pIdx, device = cairo_ps,
         width = 10, height = 10, dpi = 1200, family = "Arial")
  print(pIdx)
  
  message("wrote RSF_v", v, ".pdf and RSFindex_v", v, ".pdf  (", variantTitle[v], ")")
}