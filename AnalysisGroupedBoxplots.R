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
  r <- d$obs / d$mmean; max(r[is.finite(r)])
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
rsf_stat <- function(i, tok, p) safe({
  d  <- read.csv(paste0(i, "_myRho", type, tok, "_", p, ".csv"))
  ok <- is.finite(d$rho); Z <- d$Z[ok]
  rel <- abs(d$rho[ok] / d$ave[1] - 1)
  sum(diff(Z) * (head(rel, -1) + tail(rel, -1)) / 2)
})

# ---- assemble a long data frame for one statistic ------------------------
build_long <- function(statfn, useConds) {
  out <- list()
  for (cn in useConds) {
    tok <- condTokens[[cn]]
    for (p in percents) {
      vals <- sapply(tissues, statfn, tok = tok, p = p)
      out[[length(out) + 1]] <- data.frame(
        Condition = cn, Percent = as.character(p), Value = as.numeric(vals))
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
robust_range <- function(v, max_frac = 0.05, gap_ratio = 3) {
  # Range of the REAL data, trimming only genuinely isolated extremes: a small
  # fraction of points (<= max_frac) sitting beyond a gap larger than gap_ratio
  # times the bulk spread. This shows all of a well-behaved metric (gmax never
  # has such gaps, so nothing is cut), removes true blow-ups (RSF index), and
  # never trims a whole high condition (Max Bias is ~1/6 of the data, far more
  # than max_frac, so it is kept in full).
  v <- sort(v[is.finite(v)]); n <- length(v)
  if (n < 10) return(range(v))
  spread <- diff(quantile(v, c(0.10, 0.90), names = FALSE))   # scale of the bulk
  if (!is.finite(spread) || spread <= 0) spread <- diff(range(v))
  if (spread <= 0) spread <- 1
  thr <- gap_ratio * spread
  hi <- v[n]; ktop <- ceiling(n * (1 - max_frac))
  if (ktop < n) for (k in n:(ktop + 1)) if (v[k] - v[k - 1] > thr) hi <- v[k - 1]
  lo <- v[1]; kbot <- floor(n * max_frac)
  if (kbot >= 1) for (k in 1:kbot) if (v[k + 1] - v[k] > thr) lo <- v[k + 1]
  c(lo, hi)
}
dyn_ylim <- function(df, bottom_pad = 0.05, top_pad = 0.20) {
  v <- df$Value[is.finite(df$Value)]
  if (length(v) < 2) return(if (length(v)) c(v[1] - 1, v[1] + 1) else c(0, 1))
  rr <- robust_range(v); lo <- rr[1]; hi <- rr[2]
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
    refv <- df$Value[df$Condition == cn & df$Percent == refPct]
    tests <- setdiff(pcts, refPct); pv <- c()
    for (pc in tests) {
      v <- df$Value[df$Condition == cn & df$Percent == pc]
      pv <- c(pv, tryCatch(wilcox.test(v, refv, paired = TRUE, exact = FALSE)$p.value,
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
         width = 10, height = 8, dpi = 1200, family = "Arial")
  ggsave(paste0(fname, "_", type, ".pdf"), p, device = "pdf",
         width = 10, height = 8, dpi = 1200)
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
plot_grouped(build_long(rsf_stat,    condLevels),
             expression(plain(paste("RSF index"))),           "RSFindexGrouped")
# If the RSF panel still looks stretched/squashed, set its range explicitly, e.g.:
# plot_grouped(build_long(rsf_stat, condLevels),
#              expression(plain(paste("RSF index"))), "RSFindexGrouped", ylim = c(0, 80))

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
