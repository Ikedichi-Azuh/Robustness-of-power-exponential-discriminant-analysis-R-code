# ============================================================================
# make_figures_AB.R : figures for the Results of Blocks A and B
# Reads grid_final/cell_summary.csv and writes PDF (vector, for LaTeX) and PNG
# files to ./figures. Needs ggplot2 (>= 3.4, for `linewidth`). Base R otherwise.
#
# Design rules used (so the figures read at a glance and survive black & white):
#   * one colour per series, same colour for the same series in every figure;
#     blue / orange / aqua pass a colour-vision-deficiency check for lines that
#     can overlap, and every series also has its own linetype and marker;
#   * pp = percentage points; negative paired differences favour the PE rule;
#   * the zero line and the +/-0.5 pp practical-importance margin are drawn in
#     every difference plot; grids and axes are kept faint; no figure titles
#     (use LaTeX captions).
# ============================================================================

suppressPackageStartupMessages(library(ggplot2))

IN  <- "grid_final/cell_summary.csv"
OUT <- "figures"
TAU <- 0.5
dir.create(OUT, showWarnings = FALSE)

d <- read.csv(IN, stringsAsFactors = FALSE)
stopifnot(nrow(d) > 0)

## ---- labels and palette -------------------------------------------------------
ink  <- "#0b0b0b"; ink2 <- "#52514e"; faint <- "#e6e5e1"
blue <- "#2a78d6"; orange <- "#eb6834"; aqua <- "#1baf7a"

d$setting <- ifelse(d$block == "A",
                    paste0("Shape ", d$beta1),
                    paste0("Shapes (", d$beta1, ", ", d$beta2, ")"))
d$setting <- factor(d$setting, levels = c(
  "Shape 0.5", "Shape 1", "Shape 2",
  "Shapes (0.5, 1)", "Shapes (1, 2)", "Shapes (0.5, 2)"))
d$cov_lab <- factor(ifelse(d$cov_type == "equal", "Equal\ncovariance", "Unequal\ncovariance"),
                    levels = c("Equal\ncovariance", "Unequal\ncovariance"))

theme_pe <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = faint, linewidth = 0.3),
      axis.text          = element_text(colour = ink2, size = 11),
      axis.title         = element_text(colour = ink, size = 12),
      strip.text         = element_text(colour = ink, hjust = 0, size = 12, face = "bold"),
      strip.text.y       = element_text(colour = ink, hjust = 0.5, size = 12, face = "bold"),
      legend.position    = "top",
      legend.title       = element_blank(),
      legend.text        = element_text(colour = ink, size = 11),
      legend.key.width   = grid::unit(1.0, "cm"),
      panel.spacing      = grid::unit(0.8, "lines")
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

# Mean over cells of a paired difference, with its Monte Carlo standard error
# (cells are independent, so SE(mean) = sqrt(sum(mcse^2)) / k).
agg_diff <- function(df, by, diff, mcse) {
  key <- do.call(paste, c(df[by], sep = "|"))
  out <- do.call(rbind, lapply(split(df, key), function(g) {
    data.frame(g[1, by, drop = FALSE],
               diff = mean(g[[diff]]),
               se   = sqrt(sum(g[[mcse]]^2)) / nrow(g),
               k    = nrow(g))
  }))
  rownames(out) <- NULL
  out
}

p_cols  <- c("2" = blue, "5" = orange, "10" = aqua)
p_shape <- c("2" = 16,   "5" = 17,     "10" = 15)
p_line  <- c("2" = "solid", "5" = "dashed", "10" = "dotted")
dodge   <- position_dodge(width = 0.45)

## ---- Figure 1: common-covariance rules, equal covariance, Block A -------------
# Mean excess error against n. Shows that estimating the shape helps a little
# when the data are not normal, costs a little at beta = 1 and small n, and that
# knowing the true shape adds almost nothing.
A_eq <- subset(d, block == "A" & cov_type == "equal")
f1 <- aggregate(cbind(PE1 = excess_pp_PE1_pool_fallback,
                      LDA = excess_pp_LDA,
                      TRU = excess_pp_PE4_pool) ~ setting + n,
                data = A_eq, FUN = mean)
rule_lv <- c("PE, estimated shape", "LDA (shape fixed at 1)",
             "PE, true shape (benchmark)")
f1l <- rbind(
  data.frame(f1[c("setting", "n")], rule = rule_lv[1], value = f1$PE1),
  data.frame(f1[c("setting", "n")], rule = rule_lv[2], value = f1$LDA),
  data.frame(f1[c("setting", "n")], rule = rule_lv[3], value = f1$TRU))
f1l$rule <- factor(f1l$rule, levels = rule_lv)

fig1 <- ggplot(f1l, aes(n, value, colour = rule, linetype = rule, shape = rule)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ setting, nrow = 1) +
  scale_x_log10(breaks = c(50, 100, 250, 500)) +
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2)) +
  scale_colour_manual(values = setNames(c(blue, orange, aqua), rule_lv)) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "dotted"), rule_lv)) +
  scale_shape_manual(values = setNames(c(16, 17, 15), rule_lv)) +
  labs(x = "Training sample size per class (n)",
       y = "Mean excess error over\nthe Bayes rule (pp)") +
  theme_pe()
save_fig(fig1, "fig1_pooled_excess", 7, 3)

## ---- Figure 2: class-specific PE (estimated shape) minus QDA, Block A ---------
A <- subset(d, block == "A")
g2 <- agg_diff(A, c("setting", "n", "p"), "RQ2_sep_diff_pp", "RQ2_sep_mcse_pp")
g2$n <- factor(g2$n); g2$p <- factor(g2$p)

fig2 <- ggplot(g2, aes(n, diff, colour = p, shape = p, linetype = p, group = p)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_hline(yintercept = c(-TAU, TAU), colour = ink2, linewidth = 0.3, linetype = "dotted") +
  geom_errorbar(aes(ymin = diff - 1.96 * se, ymax = diff + 1.96 * se),
                width = 0, linewidth = 0.4, position = dodge) +
  geom_line(linewidth = 0.55, position = dodge) +
  geom_point(size = 1.8, position = dodge) +
  facet_wrap(~ setting, nrow = 1) +
  scale_colour_manual(values = p_cols, name = "Dimension p") +
  scale_shape_manual(values = p_shape, name = "Dimension p") +
  scale_linetype_manual(values = p_line, name = "Dimension p") +
  labs(x = "Training sample size per class (n)",
       y = "PE (estimated shape) minus QDA (pp)\nnegative favours PE") +
  theme_pe() + theme(legend.title = element_text(colour = ink))
save_fig(fig2, "fig2_blockA_class_specific", 7, 3)

## ---- Figure 3: same contrast in Block B (unequal shapes) ----------------------
B <- subset(d, block == "B")
g3 <- agg_diff(B, c("setting", "n", "p"), "RQ2_sep_diff_pp", "RQ2_sep_mcse_pp")
g3$n <- factor(g3$n); g3$p <- factor(g3$p)

fig3 <- ggplot(g3, aes(n, diff, colour = p, shape = p, linetype = p, group = p)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_hline(yintercept = c(-TAU, TAU), colour = ink2, linewidth = 0.3, linetype = "dotted") +
  geom_errorbar(aes(ymin = diff - 1.96 * se, ymax = diff + 1.96 * se),
                width = 0, linewidth = 0.4, position = dodge) +
  geom_line(linewidth = 0.55, position = dodge) +
  geom_point(size = 1.8, position = dodge) +
  facet_wrap(~ setting, nrow = 1) +
  scale_colour_manual(values = p_cols, name = "Dimension p") +
  scale_shape_manual(values = p_shape, name = "Dimension p") +
  scale_linetype_manual(values = p_line, name = "Dimension p") +
  labs(x = "Training sample size per class (n)",
       y = "PE (estimated shapes) minus QDA (pp)\nnegative favours PE") +
  theme_pe() + theme(legend.title = element_text(colour = ink))
save_fig(fig3, "fig3_blockB_class_specific", 7, 3)

## ---- Figure 4: decision maps (cell level) --------------------------------------
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
  ggplot(df, aes(nf, row, fill = cls)) +
    geom_tile(colour = "#d0cfc9", linewidth = 0.25, show.legend = FALSE) +
    geom_point(data = key_df, shape = 22, size = 4.5, colour = "#d0cfc9", alpha = 0) +
    geom_text(data = subset(df, cls == "Inconclusive"), aes(label = "?"),
              size = 3.6, colour = ink, show.legend = FALSE) +
    facet_grid(cov_lab ~ setting) +
    scale_fill_manual(values = fills, drop = FALSE, name = NULL) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = "Training sample size per class (n)", y = "Dimension and class separation") +
    theme_pe() + theme(panel.grid = element_blank())
}

fig4a <- decision_map(subset(d, block == "A"), "RQ2_sep_diff_pp", "RQ2_sep_mcse_pp")
save_fig(fig4a, "fig4a_decision_map_blockA", 7, 5)
fig4b <- decision_map(subset(d, block == "B"), "RQ2_sep_diff_pp", "RQ2_sep_mcse_pp")
save_fig(fig4b, "fig4b_decision_map_blockB", 7, 3.6)

## ---- Figure 5: accuracy of the shape estimate (Block A, pooled estimate) ------
s5 <- aggregate(cbind(Bias = bias_pool, RMSE = rmse_pool) ~ setting + n + p,
                data = subset(d, block == "A"), FUN = mean)
s5l <- rbind(
  data.frame(s5[c("setting", "n", "p")], measure = "Bias", value = s5$Bias),
  data.frame(s5[c("setting", "n", "p")], measure = "RMSE", value = s5$RMSE))
s5l$measure <- factor(s5l$measure, levels = c("Bias",
                                              "RMSE"))
s5l$n <- factor(s5l$n); s5l$p <- factor(s5l$p)

fig5 <- ggplot(s5l, aes(n, value, colour = p, shape = p, linetype = p, group = p)) +
  geom_hline(data = data.frame(measure = factor("Bias",
                                                levels = levels(s5l$measure)), y = 0),
             aes(yintercept = y), colour = ink2, linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8) +
  facet_grid(measure ~ setting, scales = "free_y") +
  scale_colour_manual(values = p_cols, name = "Dimension p") +
  scale_shape_manual(values = p_shape, name = "Dimension p") +
  scale_linetype_manual(values = p_line, name = "Dimension p") +
  labs(x = "Training sample size per class (n)", y = NULL) +
  theme_pe() + theme(legend.title = element_text(colour = ink))
save_fig(fig5, "fig5_shape_estimation", 7, 4.2)

## ---- Figure 6: how often the class-specific PE fit falls back to QDA -----------
f6 <- aggregate(fallback_sep ~ setting + n + p, data = d, FUN = mean)
f6$nf <- factor(f6$n); f6$pf <- factor(f6$p)
f6$lab <- sprintf("%.0f", 100 * f6$fallback_sep)
f6$txt <- ifelse(f6$fallback_sep > 0.5, "white", "black")

fig6 <- ggplot(f6, aes(nf, pf, fill = fallback_sep)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab, colour = txt), size = 3.6) +
  scale_colour_identity() +
  facet_wrap(~ setting, nrow = 2) +
  scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = c(0, 1),
                      labels = function(x) paste0(100 * x, "%"),
                      name = "Fallback to the normal rule") +
  labs(x = "Training sample size per class (n)", y = "Dimension p") +
  theme_pe() + theme(panel.grid = element_blank(),
                     legend.title = element_text(colour = ink),
                     legend.key.width = grid::unit(1.6, "cm"))
save_fig(fig6, "fig6_fallback_rates", 7, 4)

cat("\nAll figures written to", normalizePath(OUT), "\n")
