# ============================================================
# Grouped box plots across subsampling PERCENTAGES (1,5,20,50,75)
# x = Exp condition groups; within each group, one box per percentage.
# Metrics: Maximum RDF (gmax), Jmin, KS statistic, Medians, RSF index.
# Style matches the previous figures (thick axes, big fonts, EPS+PDF, Arial).
# NO per-plot legend -- a separate percentage legend is written at the end.
# Self-contained: reads the raw per-tissue/condition/percent CSVs directly.
# Set  type <- "N"  to produce the BrdU- versions instead.
# ============================================================
library(ggplot2)

# Set the working directory to this script's folder (where the CSVs live),
# exactly as the other Analysis*.R scripts do. Without this, read.csv() finds
# nothing and every box comes out empty.
current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

type <- "P"                       # "P" = BrdU+, "N" = BrdU-
tissues  <- 1:99
percents <- c(1, 5, 20, 50, 75)
pctLevels <- as.character(percents)
refPct <- "75"                    # post-hoc reference (the most-points estimate)

# Significance marking:
#   FALSE = label ONLY the non-significant comparisons with "ns" (declutters:
#           since most % differ from the 75% reference, marking the exceptions
#           is clearer than a star on nearly every box).
#   TRUE  = also draw */**/*** on the significant comparisons.
mark_significant <- FALSE
fig_w <- 12; fig_h <- 6  # grouped box-plot figure size (inches; more rectangular = wider:taller)
label_size <- 5          # text size of the "ns" / significance labels

# Exp conditions: display name -> file token
condTokens <- c("Max Bias" = "3D", "Exp30" = "Exp30", "Exp90" = "Exp90",
                "Exp180" = "Exp180", "Exp320" = "Exp320", "Random" = "R")
condLevels   <- names(condTokens)            # 6 groups (gmax/Jmin/Medians/RSF)
condLevelsKS <- names(condTokens)[1:5]       # KS has no Random box

# percentage palette (box fills shown at alpha 0.7, as in the example)
pct_pal <- c("1" = "black", "5" = "darkgreen", "20" = "orange",
             "50" = "blue", "75" = "red")

# ---- per-sample statistic readers (depend on type P/N) -------------------
safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)

gmax_stat <- function(i, tok, p) safe({
  d <- read.csv(paste0(i, "_myenvRDF", type, tok, "_", p, ".csv"))
  in200 <- is.finite(d$r) & d$r <= 200          # assess the RDF peak only up to 200 um
  r <- (d$obs / d$mmean)[in200]
  max(r[is.finite(r)])
})
jmin_stat <- function(i, tok, p) safe({
  d <- read.csv(paste0(i, "_myenvJcross", type, tok, "_", p, ".csv"))
  (d$obs / d$mmean)[120]
})
median_stat <- function(i, tok, p) safe({
  median(read.csv(paste0(i, "_myNN", type, tok, "_", p, ".csv"))[[1]])
})
ks_stat <- function(i, tok, p) safe({
  x <- read.csv(paste0(i, "_myNN", type, tok, "_", p, ".csv"))[[1]]
  r <- read.csv(paste0(i, "_myNN", type, "R_", p, ".csv"))[[1]]
  as.numeric(ks.test(x, r)$statistic)
})
# RSF index, matching AnalysisRho's myRsfIndexOne exactly: read the band file's
# curve (obs) over r, baseline = mean(obs), trapezoidal integral of
# |obs/mean(obs) - 1| dz. rsf_stat_v(v) returns the reader for band variant v.
rsf_stat_v <- function(v) function(i, tok, p) safe({
  d   <- read.csv(paste0(i, "_myRhoPbandv", v, tok, "_", p, ".csv"))  # band files are BrdU+ only
  rho <- d$obs; zz <- d$r
  ok  <- is.finite(rho) & is.finite(zz); rho <- rho[ok]; zz <- zz[ok]
  if (length(rho) < 3) return(NA_real_)
  ord <- order(zz); zz <- zz[ord]; rho <- rho[ord]     # ensure sorted grid
  rhoBar <- mean(rho)                                   # baseline = mean of the curve
  val <- abs(rho / rhoBar - 1)
  sum((val[-1] + val[-length(val)]) / 2 * diff(zz))     # trapezoidal, unitless
})

# ---- assemble a long data frame for one statistic ------------------------
build_long <- function(statfn, useConds) {
  out <- list()
  for (cn in useConds) {
    tok <- condTokens[[cn]]
    for (p in percents) {
      vals <- sapply(tissues, statfn, tok = tok, p = p)
      out[[length(out) + 1]] <- data.frame(
        Tissue = tissues, Condition = cn, Percent = as.character(p), Value = as.numeric(vals))
    }
  }
  df <- do.call(rbind, out)
  df$Condition <- factor(df$Condition, levels = useConds)
  df$Percent   <- factor(df$Percent,   levels = pctLevels)
  df
}

# ---- data-driven y-limits: span the actual data (all points are drawn),
#      pad the bottom a little and leave headroom on top for the labels.
#      This is what stops boxes/points/"ns" from being clipped.
dyn_ylim <- function(df, top_pad = 0.05, bottom_pad = 0.05) {
  # Simple, exact: use the data's own min and max, add 5% headroom each side.
  # (Runs on the real per-tissue values when you execute the script.) If one
  # metric ever has a genuine blow-up you don't want on-axis, pass an explicit
  # ylim = c(lo, hi) to plot_grouped() to override this.
  v <- df$Value[is.finite(df$Value)]
  if (!length(v)) return(c(0, 1))
  lo <- min(v); hi <- max(v)
  r <- hi - lo; if (r <= 0) r <- max(abs(hi), 1)
  c(lo - bottom_pad * r, hi + top_pad * r)
}

# ---- post-hoc significance: each % vs the reference % (paired Wilcoxon,
#      Holm-corrected within each Exp group).
#      Default: a single "ns" above each NON-significant box (mark_significant
#      = FALSE). If mark_significant = TRUE, significant boxes get */**/*** too.
myAsterics <- function(p) {
  if (is.na(p) || p >= 0.05) return("ns")
  else if (p > 0.005)  return("*")
  else if (p > 0.0005) return("**")
  else                 return("***")
}
make_stars <- function(df, ylim, frac = 0.05) {
  off <- frac * (ylim[2] - ylim[1])
  conds <- levels(df$Condition); pcts <- levels(df$Percent); n <- length(pcts)
  rows <- list()
  for (ci in seq_along(conds)) {
    cn   <- conds[ci]
    refdf <- df[df$Condition == cn & df$Percent == refPct, c("Tissue", "Value")]
    tests <- setdiff(pcts, refPct); pv <- c()
    for (pc in tests) {
      cur <- df[df$Condition == cn & df$Percent == pc, c("Tissue", "Value")]
      m   <- merge(cur, refdf, by = "Tissue")             # pair by Tissue ID: only tissues in BOTH
      ok  <- is.finite(m$Value.x) & is.finite(m$Value.y)   # keep complete, finite pairs
      pv <- c(pv, tryCatch(
        if (sum(ok) >= 2) wilcox.test(m$Value.x[ok], m$Value.y[ok],
                                      paired = TRUE, exact = FALSE)$p.value else NA_real_,
        error = function(e) NA_real_))
    }
    pv <- p.adjust(pv, method = "holm")            # correct within the group
    for (j in seq_along(tests)) {
      if (is.na(pv[j])) next                        # untestable -> no label
      sig <- pv[j] < 0.05
      if (sig && !mark_significant) next            # declutter: skip significant
      lab <- if (sig) myAsterics(pv[j]) else "ns"
      k  <- match(tests[j], pcts)
      xx <- ci + 0.8 * ((k - (n + 1) / 2) / n)     # dodged x-position of the box
      v  <- df$Value[df$Condition == cn & df$Percent == tests[j]]
      v  <- v[is.finite(v)]
      if (length(v) < 1) next
      top <- min(max(v), ylim[2])                       # top of THIS box's visible data
      yy  <- min(top + off, ylim[2] - 0.5 * off)        # place "ns" just above it, inside the panel
      rows[[length(rows) + 1]] <- data.frame(x = xx, y = yy, lab = lab)
    }
  }
  if (length(rows)) do.call(rbind, rows) else
    data.frame(x = numeric(0), y = numeric(0), lab = character(0))
}

# ---- grouped box-plot in the house style (no legend) ---------------------
plot_grouped <- function(df, ylab_expr, fname, ylim = NULL) {
  if (is.null(ylim)) ylim <- dyn_ylim(df)   # pass ylim = c(lo, hi) to override                    # dynamic axis from the data
  stars <- make_stars(df, ylim)
  p <- ggplot(df, aes(x = Condition, y = Value, fill = Percent)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black",
                 size = 0.7, position = position_dodge(width = 0.8), width = 0.7) +
    geom_point(aes(color = Percent), show.legend = FALSE, shape = 1, size = 1.8, stroke = 0.7, alpha = 0.9,
               position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8)) +
    geom_text(data = stars, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, size = label_size, fontface = "bold") +
    scale_fill_manual(values = pct_pal) +
    scale_color_manual(values = pct_pal) +
    labs(y = ylab_expr) +
    coord_cartesian(ylim = ylim) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 1),
      axis.line.y = element_line(color = "black", size = 1),
      axis.ticks.y = element_line(color = "black", size = 1),
      axis.ticks.x = element_line(color = "black", size = 1),
      plot.margin = margin(5.5, 5.5, 5.5, 25, "pt"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 44, face = "plain", color = "black"),
      axis.text.y = element_text(size = 44, face = "plain", color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 44, face = "plain"),
      legend.position = "none"
    )
  ggsave(paste0(fname, "_", type, ".eps"), p, device = cairo_ps,
         width = fig_w, height = fig_h, dpi = 1200, family = "Arial")
  ggsave(paste0(fname, "_", type, ".pdf"), p, device = "pdf",
         width = fig_w, height = fig_h, dpi = 1200)
  p
}

# ---- build + save each metric (y-limits now computed from the data) ------
plot_grouped(build_long(gmax_stat,   condLevels),
             expression(plain("Maximum RDF")),                "gmaxGrouped")
plot_grouped(build_long(jmin_stat,   condLevels),
             expression(plain("Jmin")),                       "JminGrouped")
plot_grouped(build_long(median_stat, condLevels),
             expression(plain(paste("Medians (", mu, "m)"))), "MediansGrouped")
plot_grouped(build_long(ks_stat,     condLevelsKS),
             expression(plain("KS Statistic")),               "KSGrouped")
# ---- RSF index: all four band variants (v1..v4), each its own grouped panel ----
variantTitle <- c("v1: no grid, floating bw, no baseline",
                  "v2: no grid, floating bw, with baseline",
                  "v3: grid, fixed 15um bw, no baseline",
                  "v4: grid, fixed 15um bw, with baseline")
for (v in 1:4) {
  rsf_index_df  <- build_long(rsf_stat_v(v), condLevels)   # index for band variant v
  vv            <- rsf_index_df$Value[is.finite(rsf_index_df$Value)]
  rsf_index_max <- if (length(vv)) max(vv) else 1
  rsf_index_ylim <- c(0, rsf_index_max + rsf_index_max * 0.05)   # 0 .. max + 5%
  plot_grouped(rsf_index_df,
               expression(plain("RSF index")),
               paste0("RSFindexGrouped_v", v),
               ylim = rsf_index_ylim)
  message("wrote RSFindexGrouped_v", v, "_", type, "  (", variantTitle[v], ")")
}

# ============================================================
# Separate percentage legend (its own file, used once for the whole figure)
# ============================================================
pctLab <- paste0("Proliferation ", pctLevels, "%")
n <- length(pctLevels)
ys <- rev(seq_len(n))                                  # 1% on top, 75% at bottom
legdf <- data.frame(y = ys, col = unname(pct_pal[pctLevels]), lab = pctLab)
pleg <- ggplot(legdf) +
  geom_point(aes(x = 1, y = y, fill = col), shape = 22, size = 13,
             color = "black", stroke = 0.9, alpha = 0.7) +
  scale_fill_identity() +
  geom_text(aes(x = 1.25, y = y, label = lab), hjust = 0, size = 7.5) +
  annotate("text", x = 0.78, y = n + 1, label = "Dataset",
           fontface = "bold", size = 9, hjust = 0) +
  coord_cartesian(xlim = c(0.7, 4.4), ylim = c(0.4, n + 1.6), clip = "off") +
  theme_void()
ggsave("PercentLegend.eps", pleg, device = cairo_ps, width = 4.6, height = 4, dpi = 1200, family = "Arial")
ggsave("PercentLegend.pdf", pleg, device = "pdf",      width = 4.6, height = 4, dpi = 1200)