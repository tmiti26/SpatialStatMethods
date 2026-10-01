# ============================================================================
# Sample figure set — reads the per-tissue CSVs produced by the analysis
# pipeline and makes the standard plots: for each summary function a CURVE
# across conditions (+ the Sample), plus its scalar BOX plot (Max Bias /
# Random / Sample). BrdU+ only. No BrdU-, no G/L-cross.
#
# Figures written (PDF + EPS, 10x10, Arial):
#   RSF_sample        + RSFindex_sample     (Resource Selection Function)
#   RDF_sample        + gmax_sample         (Radial Distribution Function)
#   Jfunct_sample     + J60_sample          (cross J function)
#   NNdist_sample     + Medians_sample + KS_sample  (nearest-stroma distance)
#
# House style: plain black titles 40 / tick numbers 36, thick axes+ticks 2.5,
# open-circle jitter, significance labels size 10, y-title margin.
# ============================================================================
rm(list = ls())
library(ggplot2)
library(grid)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

## ============================ CONFIG ======================================
## Edit these to match your files on disk.
numbTissue <- 99      # simulated variations per condition (i = 1..numbTissue)
N          <- 5       # specimen-number suffix in the CONDITION filenames
imyNumb    <- 5       # the Sample specimen (bare filenames <imyNumb>_<stat>.csv)

## condition display name -> filename token (the leading "P" is added by the
## reader; e.g. Max Bias RDF = <i>_myenvRDFP3D30_5.csv)
condTok <- c("Max Bias" = "3D",
             "Exp30"    = "Exp30",
             "Exp60"    = "Exp60",
             "Exp90"    = "Exp90",
             "Exp180"   = "Exp180",
             "Exp320"   = "Exp320",
             "Random"   = "R")
rsf_variant <- 4      # RSF reads the band-variant files: <i>_myRhoPbandv<rsf_variant><cond>_<N>.csv
condCurve <- names(condTok)        # order shown on the curve legend
refCond   <- "Random"              # the random reference (used for the null band + KS)

## manuscript palette (Bias magenta -> gradient -> pink; Random/Sample black)
col_map <- c("Max Bias" = "#8A0046", "Exp30" = "#C2006B", "Exp60" = "#E80074",
             "Exp90" = "#FF4D6D", "Exp180" = "#FF92C7", "Exp320" = "#FFC65C",
             "Random" = "black", "Sample" = "black")
lty_map <- setNames(rep("solid", length(col_map)), names(col_map))
lty_map["Random"] <- "dashed"      # Random drawn dashed; Sample stays solid

xMax <- 200          # x-axis cut for the curves (um); gmax also assessed <= this
J60row <- 120        # row index corresponding to r = 60 um on the J grid

## ======================= filename + io helpers ============================
fCond   <- function(pfx, cond, i) paste0(i, "_", pfx, condTok[[cond]], "_", N, ".csv")
fSample <- function(pfx)          paste0(imyNumb, "_", pfx, ".csv")
fRsf    <- function(cond, i) paste0(i, "_myRhoPbandv", rsf_variant, condTok[[cond]], "_", N, ".csv")
fRsfSmp <- function()        paste0(imyNumb, "_myRhoPbandv", rsf_variant, ".csv")
rd      <- function(f) tryCatch(read.csv(f), error = function(e) NULL)

## mean curve across tissues: gridcol = x column name, valfn(d) = y series
mean_curve <- function(pfx, cond, gridcol, valfn) {
  acc <- NULL; grid <- NULL; n <- 0
  for (i in 1:numbTissue) {
    d <- rd(fCond(pfx, cond, i)); if (is.null(d)) next
    y <- valfn(d)
    if (is.null(acc)) { grid <- d[[gridcol]]; acc <- rep(0, length(y)) }
    if (length(y) == length(acc)) { acc <- acc + y; n <- n + 1 }
  }
  if (n == 0) return(NULL)
  data.frame(x = grid, y = acc / n)
}
sample_curve <- function(pfx, gridcol, valfn) {
  d <- rd(fSample(pfx)); if (is.null(d)) return(NULL)
  data.frame(x = d[[gridcol]], y = valfn(d))
}
## per-tissue scalar vector for a condition
per_tissue <- function(pfx, cond, scalarfn)
  vapply(1:numbTissue, function(i) { d <- rd(fCond(pfx, cond, i))
  if (is.null(d)) NA_real_ else scalarfn(d) }, numeric(1))
sample_scalar <- function(pfx, scalarfn) { d <- rd(fSample(pfx))
if (is.null(d)) NA_real_ else scalarfn(d) }

## ---- startup check: do the expected files exist? (fail loudly, show reality) ----
.probe <- c(RSF        = fRsf(condCurve[1], 1),
            RDF        = fCond("myenvRDFP",    condCurve[1], 1),
            J          = fCond("myenvJcrossP", condCurve[1], 1),
            NN         = fCond("myNNP",        condCurve[1], 1),
            Sample_RSF = fRsfSmp())
.miss <- .probe[!file.exists(.probe)]
if (length(.miss)) {
  message("WARNING: these expected files were NOT found in:\n   ", getwd())
  for (nm in names(.miss)) message("   [", nm, "] looked for: ", .miss[[nm]])
  message("\nFiles actually present that look relevant:")
  print(head(sort(list.files(pattern = "_my(NNP|RhoP|envRDFP|envJcrossP)")), 30))
  stop("Filename mismatch. Copy a real name from the list above and fix ",
       "condTok / N / imyNumb in the CONFIG block, then re-run.")
}

## ======================= shared house style ===============================
base_theme <- function(x_title = TRUE) theme_minimal() + theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x  = element_line(color = "black", size = 2.5),
  axis.line.y  = element_line(color = "black", size = 2.5),
  axis.ticks.x = element_line(color = "black", size = 2.5),
  axis.ticks.y = element_line(color = "black", size = 2.5),
  axis.ticks.length = unit(0.2, "cm"),
  axis.text.x  = element_text(size = 36, face = "plain", color = "black"),
  axis.text.y  = element_text(size = 36, face = "plain", color = "black"),
  axis.title.x = if (x_title) element_text(size = 40, face = "plain") else element_blank(),
  axis.title.y = element_text(size = 40, face = "plain", margin = margin(r = 15)),
  plot.margin  = margin(14, 16, 22, 40, "pt"),
  legend.position = "none")

save_fig <- function(p, name) {
  ggsave(paste0(name, ".pdf"), p, width = 10, height = 10, dpi = 1200)
  ggsave(paste0(name, ".eps"), p, device = cairo_ps, width = 10, height = 10,
         dpi = 1200, family = "Arial")
  invisible(p)
}
myAsterics <- function(p) if (is.na(p) || p >= 0.05) "ns" else
  if (p > 0.005) "*" else if (p > 0.0005) "**" else "***"

## ---- styled curve plot (conditions + Sample, optional grey null band) ----
curve_plot <- function(df, band, ylab, name, xlim = c(0, xMax), ylim = NULL) {
  df$group <- factor(df$group, levels = c(condCurve, "Sample"))
  if (is.null(ylim)) {
    m <- df$x >= xlim[1] & df$x <= xlim[2]; v <- df$y[m & is.finite(df$y)]
    if (!is.null(band)) { mb <- band$x >= xlim[1] & band$x <= xlim[2]
    v <- c(v, band$lo[mb], band$hi[mb]) }
    v <- v[is.finite(v)]; r <- diff(range(v)); if (r <= 0) r <- max(abs(max(v)), 1)
    ylim <- c(min(v) - 0.05 * r, max(v) + 0.05 * r)
  }
  p <- ggplot()
  if (!is.null(band))
    p <- p + geom_ribbon(data = band, aes(x = x, ymin = lo, ymax = hi),
                         inherit.aes = FALSE, fill = "grey60", alpha = 0.25)
  p <- p +
    geom_line(data = df, aes(x = x, y = y, color = group, linetype = group), size = 1.5) +
    scale_color_manual(values = col_map) + scale_linetype_manual(values = lty_map) +
    labs(x = expression(plain(paste("Distance from Stroma (", mu, "m)"))), y = ylab) +
    coord_cartesian(xlim = xlim, ylim = ylim) + base_theme(TRUE)
  save_fig(p, name)
}

## ---- styled box plot: Max Bias / Random / Sample, star = Bias vs Random ----
box_plot <- function(biasV, randV, sampleV, ylab, name, ylim = NULL) {
  df <- rbind(
    data.frame(Category = "Max Bias", Value = biasV),
    data.frame(Category = "Random",   Value = randV),
    data.frame(Category = "Sample",   Value = c(sampleV,
                                                rep(NA_real_, max(0, length(biasV) - length(sampleV))))))
  df$Category <- factor(df$Category, levels = c("Max Bias", "Random", "Sample"))
  v <- df$Value[is.finite(df$Value)]
  if (is.null(ylim)) { r <- diff(range(v)); if (r <= 0) r <- max(abs(max(v)), 1)
  ylim <- c(min(v) - 0.05 * r, max(v) + 0.15 * r) }
  bg <- biasV[is.finite(biasV)]; rg <- randV[is.finite(randV)]
  pv <- tryCatch(wilcox.test(bg, rg, paired = length(bg) == length(rg),
                             exact = FALSE)$p.value, error = function(e) NA_real_)
  p <- ggplot(df, aes(x = Category, y = Value, fill = Category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = Category), width = 0.2, size = 1.8, alpha = 0.9,
                shape = 1, stroke = 0.7) +
    scale_fill_manual(values = rep("grey", 3)) +
    scale_color_manual(values = rep("black", 3)) +
    labs(y = ylab) + coord_cartesian(ylim = ylim) + base_theme(FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 36,
                                     face = "plain", color = "black")) +
    annotate("text", x = 1.5, y = ylim[2] - 0.06 * diff(ylim),
             label = myAsterics(pv), size = 10, fontface = "bold")
  save_fig(p, name)
}

## build a conditions+Sample curve frame for one (pfx, gridcol, valfn)
build_curve_df <- function(pfx, gridcol, valfn) {
  df <- do.call(rbind, lapply(condCurve, function(cc) {
    m <- mean_curve(pfx, cc, gridcol, valfn); if (is.null(m)) return(NULL)
    data.frame(x = m$x, y = m$y, group = cc) }))
  s <- sample_curve(pfx, gridcol, valfn)
  if (!is.null(s)) df <- rbind(df, data.frame(x = s$x, y = s$y, group = "Sample"))
  df
}

## ============================================================================
## 1) RSF — Resource Selection Function (obs vs distance, from band variant v)
## ============================================================================
rsf_mean <- function(cond, valcol) {           # mean of a band column across tissues
  acc <- NULL; grid <- NULL; n <- 0
  for (i in 1:numbTissue) { d <- rd(fRsf(cond, i)); if (is.null(d)) next
  y <- d[[valcol]]; if (is.null(acc)) { grid <- d$r; acc <- rep(0, length(y)) }
  if (length(y) == length(acc)) { acc <- acc + y; n <- n + 1 } }
  if (n == 0) return(NULL); data.frame(x = grid, y = acc / n)
}
rsf_df <- do.call(rbind, lapply(condCurve, function(cc) {
  m <- rsf_mean(cc, "obs"); if (is.null(m)) return(NULL)
  data.frame(x = m$x, y = m$y, group = cc) }))
sdd <- rd(fRsfSmp())
if (!is.null(sdd)) rsf_df <- rbind(rsf_df, data.frame(x = sdd$r, y = sdd$obs, group = "Sample"))
rlo <- rsf_mean(refCond, "pwLo"); rhi <- rsf_mean(refCond, "pwHi")
rsf_band <- if (!is.null(rlo)) data.frame(x = rlo$x, lo = rlo$y, hi = rhi$y) else NULL
curve_plot(rsf_df, rsf_band, expression(plain("Resource Selection Function")), "RSF_sample")

## RSF index = integral |obs/mean(obs) - 1| dr  (mean baseline, sorted grid)
rsf_index <- function(d) {
  rho <- d$obs; zz <- d$r; ok <- is.finite(rho) & is.finite(zz); rho <- rho[ok]; zz <- zz[ok]
  if (length(rho) < 3) return(NA_real_)
  o <- order(zz); zz <- zz[o]; rho <- rho[o]; rhoBar <- mean(rho)
  v <- abs(rho / rhoBar - 1); sum((v[-1] + v[-length(v)]) / 2 * diff(zz))
}
rsf_tissue <- function(cond)
  vapply(1:numbTissue, function(i) { d <- rd(fRsf(cond, i))
  if (is.null(d)) NA_real_ else rsf_index(d) }, numeric(1))
box_plot(rsf_tissue("Max Bias"), rsf_tissue("Random"),
         { d <- rd(fRsfSmp()); if (is.null(d)) NA_real_ else rsf_index(d) },
         expression(plain("RSF index")), "RSFindex_sample")

## ============================================================================
## 2) RDF — Radial Distribution Function (obs/mmean vs r)
## ============================================================================
rdf_df <- build_curve_df("myenvRDFP", "r", function(d) d$obs / d$mmean)
curve_plot(rdf_df, NULL, expression(plain("Radial Distribution Function")), "RDF_sample")

## gmax = max(obs/mmean) over r <= xMax
gmax_fn <- function(d) { in_r <- is.finite(d$r) & d$r <= xMax
r <- (d$obs / d$mmean)[in_r]; max(r[is.finite(r)]) }
box_plot(per_tissue("myenvRDFP", "Max Bias", gmax_fn),
         per_tissue("myenvRDFP", "Random", gmax_fn),
         sample_scalar("myenvRDFP", gmax_fn),
         expression(plain("Maximum RDF")), "gmax_sample")

## ============================================================================
## 3) J — cross J function (obs/mmean vs r)
## ============================================================================
j_df <- build_curve_df("myenvJcrossP", "r", function(d) d$obs / d$mmean)
curve_plot(j_df, NULL, expression(plain("J Function")), "Jfunct_sample", ylim = c(0, 1.6))

## J at 60 um = value at row J60row
j60_fn <- function(d) (d$obs / d$mmean)[J60row]
box_plot(per_tissue("myenvJcrossP", "Max Bias", j60_fn),
         per_tissue("myenvJcrossP", "Random", j60_fn),
         sample_scalar("myenvJcrossP", j60_fn),
         expression(plain(paste("J (60 ", mu, "m)"))), "J60_sample")

## ============================================================================
## 4) Nearest-stroma distance — distribution curve + medians + KS
## ============================================================================
## pooled distances per condition (for the density curve, the KS reference, medians)
nn_pool <- function(cond) { x <- numeric(0)
for (i in 1:numbTissue) { d <- rd(fCond("myNNP", cond, i))
if (!is.null(d)) x <- c(x, d[[1]]) }; x[is.finite(x)] }
nn_sample <- function() { d <- rd(fSample("myNNP"))
if (is.null(d)) numeric(0) else { v <- d[[1]]; v[is.finite(v)] } }

## distribution curve: density of pooled distances per condition (+ Sample)
nn_curve <- do.call(rbind, lapply(condCurve, function(cc) {
  x <- nn_pool(cc); if (length(x) < 3) return(NULL)
  de <- density(x, from = 0, to = xMax, n = 512)
  data.frame(x = de$x, y = de$y, group = cc) }))
sx <- nn_sample()
if (length(sx) >= 3) { de <- density(sx, from = 0, to = xMax, n = 512)
nn_curve <- rbind(nn_curve, data.frame(x = de$x, y = de$y, group = "Sample")) }
curve_plot(nn_curve, NULL, expression(plain("Frequency")), "NNdist_sample")

## medians box (per tissue)
med_fn <- function(d) median(d[[1]], na.rm = TRUE)
box_plot(per_tissue("myNNP", "Max Bias", med_fn),
         per_tissue("myNNP", "Random", med_fn),
         median(nn_sample()),
         expression(plain(paste("Median distance (", mu, "m)"))), "Medians_sample")

## KS box: KS(BrdU+ , pooled random) per tissue  (BrdU+ vs its random reference)
NNPRandtotal <- nn_pool(refCond)
ks_fn <- function(d) as.numeric(ks.test(d[[1]], NNPRandtotal)$statistic)
box_plot(per_tissue("myNNP", "Max Bias", ks_fn),
         per_tissue("myNNP", "Random", ks_fn),
         { sx <- nn_sample(); if (length(sx) < 1) NA_real_ else
           as.numeric(ks.test(sx, NNPRandtotal)$statistic) },
         expression(plain("KS statistic")), "KS_sample")

message("Done. Wrote RSF/RDF/Jfunct/NNdist curves and RSFindex/gmax/J60/Medians/KS boxes.")