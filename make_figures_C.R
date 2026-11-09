# ============================================================================
# make_figures_C.R : figures for the Block C (misspecification) Results
# Reads grid_final_C/cell_summary.csv; writes PDF + PNG to ./figures.
#   fig9   mean excess error by distribution: estimated-shape PE vs normal-based
#   fig10  the shape the PE model adopted (median and 5-95% range of beta-hat)
#   fig11a decision map, common covariance (approach 1 vs LDA)
#   fig11b decision map, class-specific covariance (approach 1 vs QDA)
# Needs ggplot2 (>= 3.4). Same colours and conventions as make_figures_AB.R.
# ============================================================================

suppressPackageStartupMessages(library(ggplot2))

IN  <- "grid_final_C/cell_summary.csv"
OUT <- "figures"
TAU <- 0.5
dir.create(OUT, showWarnings = FALSE)

d <- read.csv(IN, stringsAsFactors = FALSE)
stopifnot(nrow(d) > 0)

ink <- "#0b0b0b"; ink2 <- "#52514e"; faint <- "#e6e5e1"
blue <- "#2a78d6"; orange <- "#eb6834"

# Heaviest tails first; the skew-normal (asymmetry only) last.
dist_lv  <- c("Student t, 3 d.f.", "Contaminated normal", "Skew-t, 5 d.f.",
              "Student t, 5 d.f.", "Student t, 10 d.f.", "Skew-normal")
dist_map <- c(t3 = dist_lv[1], cn = dist_lv[2], st = dist_lv[3],
              t5 = dist_lv[4], t10 = dist_lv[5], sn = dist_lv[6])
d$dist_lab <- factor(unname(dist_map[d$dist]), levels = dist_lv)
d$cov_lab  <- factor(ifelse(d$cov_type == "equal", "Equal\ncovariance", "Unequal\ncovariance"),
                     levels = c("Equal\ncovariance", "Unequal\ncovariance"))
d$n_lab    <- factor(paste0("n = ", d$n), levels = c("n = 50", "n = 250"))

theme_pe <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = faint, linewidth = 0.3),
      axis.text        = element_text(colour = ink2, size = 11),
      axis.title       = element_text(colour = ink, size = 12),
      strip.text       = element_text(colour = ink, hjust = 0, size = 12, face = "bold"),
      strip.text.y       = element_text(colour = ink, hjust = 0.5, size = 12, face = "bold"),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(colour = ink, size = 11),
      panel.spacing    = grid::unit(0.8, "lines")
    )
}
# Larger text (base 12) needs more room: every figure is drawn 14% wider and 20% taller than
# the sizes given in the save_fig() calls. Set both to 1 to return to the old layout.
W_SCALE <- 8 / 7; H_SCALE <- 1.2
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w * W_SCALE, height = h * H_SCALE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w * W_SCALE, height = h * H_SCALE, dpi = 300)
  cat("saved", name, "\n")
}

## ---- Figure 9: mean excess error by distribution ----------------------------------
rule_lv <- c("PE, fallback",
             "PE, boundary accepted",
             "LDA / QDA")
mk <- function(family, cols) {
  a <- aggregate(d[unname(cols)], by = list(dist = d$dist_lab, n = d$n_lab), FUN = mean)
  do.call(rbind, lapply(seq_along(cols), function(i) {
    data.frame(dist = a$dist, n = a$n, family = family, rule = names(cols)[i],
               value = a[[unname(cols)[i]]])
  }))
}
f9 <- rbind(
  mk("Common covariance (vs LDA)",
     c("PE, fallback" = "excess_pp_PE1_pool_fallback",
       "PE, boundary accepted" = "excess_pp_PE1_pool_accept",
       "LDA / QDA" = "excess_pp_LDA")),
  mk("Class-specific covariance (vs QDA)",
     c("PE, fallback" = "excess_pp_PE1_sep_fallback",
       "PE, boundary accepted" = "excess_pp_PE1_sep_accept",
       "LDA / QDA" = "excess_pp_QDA")))
f9$rule   <- factor(f9$rule, levels = rule_lv)
f9$family <- factor(f9$family, levels = c("Common covariance (vs LDA)",
                                           "Class-specific covariance (vs QDA)"))
f9$dist   <- factor(f9$dist, levels = rev(dist_lv))

fig9 <- ggplot(f9, aes(value, dist, colour = rule, shape = rule)) +
  geom_point(size = 2.2, stroke = 0.7) +
  facet_grid(n ~ family) +
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 4)) +
  scale_colour_manual(values = setNames(c(blue, blue, orange), rule_lv)) +
  scale_shape_manual(values = setNames(c(16, 1, 17), rule_lv)) +
  labs(x = "Mean excess error over the Bayes rule (pp, log scale)", y = NULL) +
  theme_pe() + theme(legend.position = "top")
save_fig(fig9, "fig9_blockC_excess", 7, 4.4)

## ---- Figure 10: the shape the PE model adopted ------------------------------------
f10 <- aggregate(cbind(med = median_bhat_pool, lo = q05_bhat_pool, hi = q95_bhat_pool) ~
                   dist_lab + n_lab, data = d, FUN = mean)
f10$dist_lab <- factor(f10$dist_lab, levels = rev(dist_lv))

fig10 <- ggplot(f10, aes(med, dist_lab, shape = n_lab)) +
  geom_vline(xintercept = 1, colour = ink2, linewidth = 0.3, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = blue, linewidth = 0.5,
                 position = position_dodge(width = 0.5)) +
  geom_point(colour = blue, size = 2.2, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c("n = 50" = 16, "n = 250" = 17)) +
  scale_x_continuous(limits = c(0, 1.6), breaks = seq(0, 1.5, 0.25)) +
  labs(x = "Estimated shape", y = NULL) +
  theme_pe() + theme(axis.title.x = element_text(colour = ink, size = 12, margin = margin(t = 8)))
save_fig(fig10, "fig10_blockC_shape", 7, 3)

## ---- Figure 11: decision maps ----------------------------------------------------------
cls <- function(diff, mcse, tau = TAU) {
  lo <- diff - 1.96 * mcse; hi <- diff + 1.96 * mcse
  out <- ifelse(lo > -tau & hi < tau, "Comparable",
         ifelse(hi < -tau, "Approach 1 better",
         ifelse(lo > tau, "Approach 1 worse", "Inconclusive")))
  factor(out, levels = c("Approach 1 better", "Comparable",
                         "Approach 1 worse", "Inconclusive"))
}
fills <- c("Approach 1 better" = blue, "Comparable" = "#e4e3de",
           "Approach 1 worse" = "#e34948", "Inconclusive" = "#ffffff")

decision_map <- function(df, diff, mcse) {
  df$cls <- cls(df[[diff]], df[[mcse]])
  lv <- as.vector(t(outer(paste0("p = ", sort(unique(df$p))),
                          c("high", "moderate", "low"), paste, sep = ", ")))
  df$row <- factor(paste(paste0("p = ", df$p), df$separation, sep = ", "), levels = rev(lv))
  df$nf  <- factor(df$n)
  # The legend keys come from an invisible dummy layer, so that every category is
  # shown with its colour box even when no cell falls in it (e.g. "worse").
  key_df <- data.frame(cls = factor(names(fills), levels = names(fills)),
                       nf = levels(df$nf)[1], row = levels(df$row)[1])
  # Two-line panel titles, so that they are not clipped ("Contaminated normal").
  two_line <- function(x) sub("Contaminated normal", "Contaminated\nnormal",
                              sub(", ", "\n", x, fixed = TRUE), fixed = TRUE)
  ggplot(df, aes(nf, row, fill = cls)) +
    geom_tile(colour = "#d0cfc9", linewidth = 0.25, show.legend = FALSE) +
    geom_point(data = key_df, shape = 22, size = 4.5, colour = "#d0cfc9", alpha = 0) +
    geom_text(data = subset(df, cls == "Inconclusive"), aes(label = "?"),
              size = 3.6, colour = ink, show.legend = FALSE) +
    facet_grid(cov_lab ~ dist_lab, labeller = labeller(dist_lab = two_line)) +
    scale_fill_manual(values = fills, drop = FALSE, name = NULL) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = "Training sample size per class (n)", y = "Dimension and class separation") +
    theme_pe() + theme(panel.grid = element_blank())
}
save_fig(decision_map(d, "RQ2_pool_diff_pp", "RQ2_pool_mcse_pp"),
         "fig11a_blockC_map_common", 7.5, 4.6)
save_fig(decision_map(d, "RQ2_sep_diff_pp", "RQ2_sep_mcse_pp"),
         "fig11b_blockC_map_class_specific", 7.5, 4.6)

cat("\nAll figures written to", normalizePath(OUT), "\n")
