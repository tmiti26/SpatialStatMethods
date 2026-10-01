rm(list = ls())
library(ggplot2)
library(grid)

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

## ================= settings (shared with the other curve plots) =================
numbTissue <- 99
N <- 14

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

pink_gradient <- c("#7B2D8B", "#B5006A", "#E8003A", "#FF7800", "#FF99CC")
names(pink_gradient) <- c("BrdU+ Max Bias", "BrdU+ Exp30", "BrdU+ Exp90",
                          "BrdU+ Exp180", "BrdU+ Exp320")
color_mapping <- c(pink_gradient, "BrdU+ Random" = "black")

xLo <- 0; xHi <- myMergeXmax()
fmtP <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

## median global p-value for one condition, from pvalues_<cond>_<N>.csv.
## Self-diagnosing: if the stat name is not present it prints what IS available.
myMedianP <- function(cond, statName){
  pf <- paste0("pvalues_", cond, "_", N, ".csv")
  if (!file.exists(pf)){ message("  [p] missing file: ", pf); return(NA_real_) }
  pv <- read.csv(pf, stringsAsFactors = FALSE)
  pv$stat <- trimws(pv$stat)
  sel <- pv$stat == statName
  if (!any(sel)){
    message("  [p] '", statName, "' not found in ", pf,
            " -- available stats: ", paste(unique(pv$stat), collapse = ", "))
    return(NA_real_)
  }
  median(pv$p[sel], na.rm = TRUE)
}

## shared theme (legend ON, carries the p-values)
myCurveTheme <- function(){
  theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 1),
          axis.line.y = element_line(color = "black", size = 1),
          axis.ticks.x = element_line(color = "black", size = 1),
          axis.ticks.y = element_line(color = "black", size = 1),
          axis.ticks.length = unit(0.25, "cm"),
          text = element_text(size = 32),
          axis.title = element_text(face = "plain", size = 32),
          plot.margin = margin(14, 16, 22, 40, "pt"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 16, face = "bold"))
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

## assemble + draw one figure given per-condition curves (r, obs, lo, hi) and p-values
## split = TRUE  ->  save the graph without its legend, and write the legend to its own file
myDrawGlobal <- function(cv, pvals, yLabel, fileName, split = FALSE){
  rr <- cv[["R"]]$r
  df <- do.call(rbind, lapply(condList, function(cc)
    data.frame(x = cv[[cc]]$r, y = cv[[cc]]$obs, group = condLabel[[cc]])))
  df$group <- factor(df$group, levels = unname(condLabel[condList]))
  nullLo <- cv[["R"]]$lo; nullHi <- cv[["R"]]$hi
  labs_p <- setNames(paste0(sub("^BrdU[+/-]+ ", "", condLabel[condList]),
                            "  (p=", fmtP(pvals[condList]), ")"), condLabel[condList])

  inX <- df$x >= xLo & df$x <= xHi; inR <- rr >= xLo & rr <= xHi
  yAll <- c(df$y[inX], nullLo[inR], nullHi[inR]); yAll <- yAll[is.finite(yAll)]
  yPad <- 0.05 * (max(yAll) - min(yAll)); yLo <- min(yAll) - yPad; yHi <- max(yAll) + yPad

  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_ribbon(data = data.frame(x = rr, ymin = nullLo, ymax = nullHi),
                aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
                fill = "grey60", alpha = 0.25) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_mapping, labels = labs_p) +
    labs(x = expression(plain(paste("Distance from Stroma", " ", "(", mu, "m", ")"))), y = yLabel) +
    coord_cartesian(ylim = c(yLo, yHi), xlim = c(xLo, xHi)) +
    myCurveTheme()
  epsName <- sub("(\\.[^.]+)$", ".eps", fileName)
  if (split){
    ## ---- main graph WITHOUT the legend ----
    pMain <- p + theme(legend.position = "none")
    ggsave(fileName, pMain, width = 10, height = 10, dpi = 1200)
    ggsave(epsName,  pMain, device = cairo_ps, width = 10, height = 10, dpi = 1200, family = "Arial")
    print(pMain)
    ## ---- legend as its own image (pdf + eps), auto-sized ----
    leg <- myGetLegend(p)
    if (!is.null(leg)){
      lw <- tryCatch(convertWidth(sum(leg$widths),   "in", valueOnly = TRUE), error = function(e) NA)
      lh <- tryCatch(convertHeight(sum(leg$heights), "in", valueOnly = TRUE), error = function(e) NA)
      if (!is.finite(lw) || lw <= 0) lw <- 6
      if (!is.finite(lh) || lh <= 0) lh <- 6
      legPdf <- sub("(\\.[^.]+)$", "_legend.pdf", fileName)
      legEps <- sub("(\\.[^.]+)$", "_legend.eps", fileName)
      mySaveGrob(leg, legPdf, lw + 0.2, lh + 0.2)
      mySaveGrob(leg, legEps, lw + 0.2, lh + 0.2)
      message("wrote ", legPdf, " and ", legEps)
    }
  } else {
    ggsave(fileName, p, width = 10, height = 10, dpi = 1200)
    ggsave(epsName,  p, device = cairo_ps, width = 10, height = 10, dpi = 1200, family = "Arial")
    print(p)
  }
  message("wrote ", fileName, " and ", epsName)
}

## ---------- (1) cross statistics: RDF, J  (global envelope files, p from pvalues) ----------
## normalize obs and band by the global 'central'; interpolate onto variation-1's grid.
myReadMeanGlobCross <- function(stat, cond){
  ref <- read.csv(paste0(1, "_myenv", stat, "Pglob", cond, "_", N, ".csv"))
  rr  <- ref$r; obs <- rep(0, length(rr)); cen <- obs; lo <- obs; hi <- obs
  hasP <- "pglobal" %in% names(ref); ps <- numeric(numbTissue)
  for (i in 1:numbTissue){
    d <- read.csv(paste0(i, "_myenv", stat, "Pglob", cond, "_", N, ".csv"))
    obs <- obs + approx(d$r, d$obs,     rr, rule = 2)$y / numbTissue
    cen <- cen + approx(d$r, d$central, rr, rule = 2)$y / numbTissue
    lo  <- lo  + approx(d$r, d$lo,      rr, rule = 2)$y / numbTissue
    hi  <- hi  + approx(d$r, d$hi,      rr, rule = 2)$y / numbTissue
    if (hasP) ps[i] <- d$pglobal[1]
  }
  ## normalized by central; p taken from the glob file if that column exists
  list(r = rr, obs = obs / cen, lo = lo / cen, hi = hi / cen,
       p = if (hasP) median(ps, na.rm = TRUE) else NA_real_)
}
plotGlobCross <- function(stat, statP, yLabel, fileName){
  cv <- lapply(condList, function(cc) myReadMeanGlobCross(stat, cc)); names(cv) <- condList
  ## p from the glob file's pglobal column if present, else fall back to pvalues_<cond>.csv
  pv <- sapply(condList, function(cc)
    if (!is.na(cv[[cc]]$p)) cv[[cc]]$p else myMedianP(cc, statP)); names(pv) <- condList
  myDrawGlobal(cv, pv, yLabel, fileName)
}
plotGlobCross("RDF",    "RDF_Pos",    expression(plain("Radial Distribution Function")), "RDF_global.pdf")
plotGlobCross("Jcross", "Jcross_Pos", expression(plain("Cross J-function")),             "Jcross_global.pdf")

## ---------- (2) band statistics: NN (freq. distribution) and RSF v1-v4 ----------
## band files carry their own global band (globLo/globHi) AND p-value (pglobal).
myReadMeanBandGlob <- function(fileStem, cond){
  ref <- read.csv(paste0(1, "_", fileStem, cond, "_", N, ".csv"))
  rr <- ref$r; obs <- rep(0, length(rr)); lo <- obs; hi <- obs; ps <- numeric(numbTissue)
  for (i in 1:numbTissue){
    d <- read.csv(paste0(i, "_", fileStem, cond, "_", N, ".csv"))
    obs <- obs + approx(d$r, d$obs,    rr, rule = 2)$y / numbTissue
    lo  <- lo  + approx(d$r, d$globLo, rr, rule = 2)$y / numbTissue
    hi  <- hi  + approx(d$r, d$globHi, rr, rule = 2)$y / numbTissue
    ps[i] <- d$pglobal[1]
  }
  list(r = rr, obs = obs, lo = lo, hi = hi, p = median(ps, na.rm = TRUE))
}
plotGlobBand <- function(fileStem, yLabel, fileName, split = FALSE){
  cv <- lapply(condList, function(cc) myReadMeanBandGlob(fileStem, cc)); names(cv) <- condList
  pv <- sapply(condList, function(cc) cv[[cc]]$p); names(pv) <- condList
  myDrawGlobal(cv, pv, yLabel, fileName, split = split)
}
## nearest-stroma distance frequency distribution (ECDF), global band + p
plotGlobBand("myNNPband", expression(plain(paste("Nearest-stroma distance  F(d)"))), "NN_global.pdf")
## RSF, all four variants
for (v in 1:4)
  plotGlobBand(paste0("myRhoPbandv", v),
               expression(plain("Resource Selection Function")),
               paste0("RSF_global_v", v, ".pdf"), split = TRUE)
