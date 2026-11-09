# ============================================================================
# make_figures_distributions.R : figures of the data-generating distributions
# Uses the simulation's own generators (version 6 script), so the pictures show
# exactly what was simulated. Writes PDF + PNG to ./figures.
#
#   fig7a  PE family: marginal density of one variable for shape 0.5, 1, 2
#          (all with variance 1), linear and log scale
#   fig7b  Block C distributions: marginal density against the normal (linear)
#   fig7c  Block C distributions: the same on a log scale (tails)
#   fig8   two classes at the three calibrated separations (p = 2, equal
#          covariance), for shape 0.5, 1 and 2
# Needs ggplot2 (>= 3.4) and, for fig8, grid_final/calibration.csv.
# ============================================================================

suppressPackageStartupMessages(library(ggplot2))

SCRIPT <- "pe_discriminant_analysis_profile_optimized_6.R"
CALIB  <- "grid_final/calibration.csv"
OUT    <- "figures"
M      <- 1e6          # draws used for each marginal density
dir.create(OUT, showWarnings = FALSE)

env <- new.env()
sys.source(SCRIPT, envir = env)
set.seed(20261105)

ink <- "#0b0b0b"; ink2 <- "#52514e"; faint <- "#e6e5e1"
blue <- "#2a78d6"; orange <- "#eb6834"; aqua <- "#1baf7a"

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
      legend.key.width = grid::unit(1.2, "cm"),
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

# Kernel density of the first variable on a fixed grid (variance is 1 for every
# distribution, because each is scaled to the same covariance matrix).
grid_x <- seq(-6, 6, length.out = 601)
marg <- function(x, label) {
  dd <- density(x, bw = 0.04, from = -6, to = 6, n = length(grid_x))
  data.frame(x = dd$x, dens = dd$y, dist = label)
}
# log10 density, blanked below 1e-5: a kernel estimate from 10^6 draws is too
# noisy there to show (those parts of the line simply end).
log_fun <- function(v) {
  l <- log10(v)
  ifelse(is.finite(l) & l > -5, l, NA_real_)
}
long_scales <- function(df) {
  rbind(
    data.frame(df, scale = "Density", value = df$dens),
    data.frame(df, scale = "log10 density (tails)", value = log_fun(df$dens))
  )
}

## ---- Figure 7a: the power exponential family ------------------------------------
pe_lab <- c("0.5" = "Shape 0.5 (heavy tails, sharp peak)",
            "1"   = "Shape 1 (normal)",
            "2"   = "Shape 2 (light tails, flat top)")
pe <- do.call(rbind, lapply(c(0.5, 1, 2), function(b) {
  s <- env$make_setup(2, b, b, "equal")
  marg(env$draw_centered(s, 1L, M)[, 1], pe_lab[as.character(b)])
}))
pe$dist <- factor(pe$dist, levels = unname(pe_lab))
pe_l <- long_scales(pe)
pe_l$scale <- factor(pe_l$scale, levels = c("Density", "log10 density (tails)"))

fig7a <- ggplot(pe_l, aes(x, value, colour = dist, linetype = dist)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~ scale, nrow = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-5, 5)) +
  scale_colour_manual(values = setNames(c(blue, orange, aqua), levels(pe$dist))) +
  scale_linetype_manual(values = setNames(c("solid", "dashed", "dotted"), levels(pe$dist))) +
  labs(x = "First variable (variance 1)", y = NULL) +
  theme_pe()
save_fig(fig7a, "fig7a_pe_family", 7, 3)

## ---- Figures 7b and 7c: Block C distributions -----------------------------------
cdist <- c(t3 = "Student t, 3 d.f.", t5 = "Student t, 5 d.f.", t10 = "Student t, 10 d.f.",
           cn = "Contaminated normal", sn = "Skew-normal", st = "Skew-t, 5 d.f.")
cc <- do.call(rbind, lapply(names(cdist), function(dn) {
  s <- env$make_setup_C(2, dn, "equal")
  marg(env$draw_centered(s, 1L, M)[, 1], cdist[[dn]])
}))
cc$dist <- factor(cc$dist, levels = unname(cdist))
nrm <- data.frame(x = grid_x, dens = dnorm(grid_x))
nrm_all <- do.call(rbind, lapply(levels(cc$dist), function(l) data.frame(nrm, dist = l)))
nrm_all$dist <- factor(nrm_all$dist, levels = levels(cc$dist))

block_c_plot <- function(value_fun, ylab, ylim = NULL) {
  cc$value <- value_fun(cc$dens); nrm_all$value <- value_fun(nrm_all$dens)
  ggplot() +
    geom_line(data = nrm_all, aes(x, value), colour = ink2, linewidth = 0.4,
              linetype = "dashed", na.rm = TRUE) +
    geom_line(data = cc, aes(x, value), colour = blue, linewidth = 0.7, na.rm = TRUE) +
    facet_wrap(~ dist, nrow = 2) +
    coord_cartesian(xlim = c(-5, 5), ylim = ylim) +
    labs(x = "First variable (variance 1); dashed: normal with the same variance",
         y = ylab) +
    theme_pe() + theme(legend.position = "none")
}
save_fig(block_c_plot(identity, "Density"), "fig7b_blockC_densities", 7, 4)
save_fig(block_c_plot(log_fun, "log10 density", ylim = c(-5, 0)),
         "fig7c_blockC_tails", 7, 4)

## ---- Figure 8: two classes at the calibrated separations ------------------------
cal <- read.csv(CALIB, stringsAsFactors = FALSE)
get_theta <- function(b, tgt) {
  th <- cal$theta[cal$key == paste("2", b, b, "equal", tgt, sep = "|")]
  stopifnot(length(th) == 1L, is.finite(th))
  th
}
sep_lab <- c("0.05" = "5% Bayes error\n(high separation)",
             "0.15" = "15% Bayes error\n(moderate separation)",
             "0.3"  = "30% Bayes error\n(low separation)")
sc <- do.call(rbind, unlist(lapply(c(0.5, 1, 2), function(b) {
  lapply(c(0.05, 0.15, 0.30), function(tg) {
    s  <- env$make_setup(2, b, b, "equal")
    tr <- env$gen_data(s, 300, 300, get_theta(b, tg))
    data.frame(x1 = tr$X[, 1], x2 = tr$X[, 2],
               class = ifelse(tr$y == 1, "Class 1", "Class 2"),
               shape = paste("Shape", b), sep = sep_lab[as.character(tg)])
  })
}), recursive = FALSE))
sc$shape <- factor(sc$shape, levels = c("Shape 0.5", "Shape 1", "Shape 2"))
sc$sep   <- factor(sc$sep, levels = unname(sep_lab))

fig8 <- ggplot(sc, aes(x1, x2, colour = class, shape = class)) +
  geom_point(size = 0.9, alpha = 0.55) +
  facet_grid(sep ~ shape) +
  coord_fixed(xlim = c(-5, 8), ylim = c(-5, 8)) +
  scale_colour_manual(values = c("Class 1" = blue, "Class 2" = orange)) +
  scale_shape_manual(values = c("Class 1" = 16, "Class 2" = 17)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Variable 1", y = "Variable 2") +
  theme_pe()
save_fig(fig8, "fig8_two_classes", 7, 6.2)

cat("\nAll figures written to", normalizePath(OUT), "\n")
