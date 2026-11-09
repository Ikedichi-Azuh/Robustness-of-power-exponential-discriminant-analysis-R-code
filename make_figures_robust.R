# ============================================================================
# make_figures_robust.R : figures for the robust-comparator Results
# Reads robust_results/robust_summary.csv; writes PDF + PNG to ./figures.
#   fig12  mean excess error: PE, normal-based, MCD, multivariate t
#          (Block A power exponential scenarios and the six Block C distributions)
#   fig13  cell-level difference PE minus t (pp); dashed lines: +/- 0.5 pp margin
# Needs ggplot2 (>= 3.4). Same colours and conventions as make_figures_C.R.
# ============================================================================

suppressPackageStartupMessages(library(ggplot2))

IN  <- "robust_results/robust_summary.csv"
OUT <- "figures"
TAU <- 0.5
dir.create(OUT, showWarnings = FALSE)

d <- read.csv(IN, stringsAsFactors = FALSE)
stopifnot(nrow(d) > 0)

ink <- "#0b0b0b"; ink2 <- "#52514e"; faint <- "#e6e5e1"
blue <- "#2a78d6"; orange <- "#eb6834"; aqua <- "#1baf7a"; plum <- "#8a4fb3"

grp_lv <- c("Power exponential, shape 0.5", "Power exponential, shape 1",
            "Power exponential, shape 2",
            "Student t, 3 d.f.", "Contaminated normal", "Skew-t, 5 d.f.",
            "Student t, 5 d.f.", "Student t, 10 d.f.", "Skew-normal")
grp_map <- c(A_beta0.5 = grp_lv[1], A_beta1 = grp_lv[2], A_beta2 = grp_lv[3],
             t3 = grp_lv[4], cn = grp_lv[5], st = grp_lv[6],
             t5 = grp_lv[7], t10 = grp_lv[8], sn = grp_lv[9])
key <- ifelse(d$source == "C", d$dist, paste0("A_beta", d$beta1))
d$grp   <- factor(unname(grp_map[key]), levels = rev(grp_lv))
d$n_lab <- factor(paste0("n = ", d$n), levels = c("n = 50", "n = 250"))
stopifnot(!anyNA(d$grp))

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

## ---- Figure 12: mean excess error of the four rules ------------------------------------
rule_lv <- c("PE, estimated shape", "LDA / QDA",
             "MCD", "Multivariate t")
mk <- function(family, cols) {
  a <- aggregate(d[unname(cols)], by = list(grp = d$grp, n = d$n_lab), FUN = mean)
  do.call(rbind, lapply(seq_along(cols), function(i) {
    data.frame(grp = a$grp, n = a$n, family = family, rule = names(cols)[i],
               value = a[[unname(cols)[i]]])
  }))
}
f12 <- rbind(
  mk("Common covariance (vs LDA)",
     c("PE, estimated shape"          = "excess_pp_PE1_pool_fallback",
       "LDA / QDA"    = "excess_pp_LDA",
       "MCD"        = "excess_pp_MCD_pool",
       "Multivariate t"               = "excess_pp_T_pool")),
  mk("Class-specific covariance (vs QDA)",
     c("PE, estimated shape"          = "excess_pp_PE1_sep_fallback",
       "LDA / QDA"    = "excess_pp_QDA",
       "MCD"        = "excess_pp_MCD_sep",
       "Multivariate t"               = "excess_pp_T_sep")))
f12$rule   <- factor(f12$rule, levels = rule_lv)
f12$family <- factor(f12$family, levels = c("Common covariance (vs LDA)",
                                             "Class-specific covariance (vs QDA)"))

fig12 <- ggplot(f12, aes(value, grp, colour = rule, shape = rule)) +
  geom_point(size = 2.2, stroke = 0.7) +
  facet_grid(n ~ family) +
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 4)) +
  scale_colour_manual(values = setNames(c(blue, orange, ink2, plum), rule_lv)) +
  scale_shape_manual(values = setNames(c(16, 17, 15, 18), rule_lv)) +
  labs(x = "Mean excess error over the Bayes rule (pp, log scale)", y = NULL) +
  theme_pe()
save_fig(fig12, "fig12_robust_excess", 7.5, 5.2)

## ---- Figure 13: PE minus multivariate t, cell by cell ---------------------------------------
f13 <- rbind(
  data.frame(grp = d$grp, n = d$n_lab, family = "Common covariance",
             diff = d$PE_vs_T_pool_diff_pp),
  data.frame(grp = d$grp, n = d$n_lab, family = "Class-specific covariance",
             diff = d$PE_vs_T_sep_diff_pp))
f13$family <- factor(f13$family, levels = c("Common covariance", "Class-specific covariance"))

fig13 <- ggplot(f13, aes(diff, grp, colour = n, shape = n)) +
  geom_vline(xintercept = c(-TAU, TAU), colour = ink2, linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_point(size = 1.4, alpha = 0.7, position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  facet_wrap(~ family, nrow = 1) +
  scale_colour_manual(values = c("n = 50" = blue, "n = 250" = orange)) +
  scale_shape_manual(values = c("n = 50" = 16, "n = 250" = 17)) +
  labs(x = "Difference in excess error, PE minus multivariate t (pp; negative: PE better)",
       y = NULL) +
  theme_pe()
save_fig(fig13, "fig13_robust_PE_minus_t", 7.5, 4)

cat("\nAll figures written to", normalizePath(OUT), "\n")
