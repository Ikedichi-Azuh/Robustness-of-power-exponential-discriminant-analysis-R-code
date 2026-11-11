# ============================================================================
# analyze_fallback_cor.R : Does the association between shape-estimation error 
# and the paired difference (estimated shape vs fixed beta = 1, approach 1 minus
# approach 2) depend on the replicates in which the estimated-shape fit was
# replaced by the normal-based fit (fallback)?
#
# Reads the stored replicate files grid_final/cells/*B1000.rds (run from the
# folder that contains grid_final/). Changes nothing. Paste everything printed.
#
#  PART 1  within-cell Spearman correlation across replicates, medians by group:
#          (a) all replicates (should reproduce the numbers now in the paper)
#          (b) non-fallback replicates only
#  PART 2  cell-level Spearman correlation (Block A, 54 cells per n) between the
#          RMSE of the pooled shape estimate and the cost of estimating the shape
#          (approach 1 minus approach 4): all replicates vs non-fallback only
# ============================================================================
options(width = 170, digits = 3)
MIN_REP <- 30                       # minimum replicates for a within-cell correlation
files <- list.files("grid_final/cells", pattern = "B1000\\.rds$", full.names = TRUE)
cat("cell files:", length(files), "\n")
stopifnot(length(files) > 0)

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < MIN_REP || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

rows <- lapply(files, function(f) {
  x <- readRDS(f); z <- x$spec; r <- x$res
  # approach 1 minus approach 2, in pp
  d_pool <- 100 * (r$PE1_pool_fallback - r$PE2_pool)
  d_sep  <- 100 * (r$PE1_sep_fallback  - r$PE2_sep)
  # absolute shape-estimation error (raw estimates, before any fallback)
  e_pool <- abs(r$bhat_pool - z$beta1)
  e_sep  <- (abs(r$bhat1_sep - z$beta1) + abs(r$bhat2_sep - z$beta2)) / 2
  fb_pool <- r$fallback_pool_fallback == 1
  fb_sep  <- r$fallback_sep_fallback  == 1
  ok_pool <- !is.na(fb_pool) & !fb_pool
  ok_sep  <- !is.na(fb_sep)  & !fb_sep
  same <- z$beta1 == z$beta2
  # cost of estimating the shape (approach 1 minus approach 4), pooled, equal shapes only
  cost_all <- if (same) 100 * mean(r$PE1_pool_fallback - r$PE4_pool, na.rm = TRUE) else NA_real_
  cost_nf  <- if (same && any(ok_pool))
                100 * mean((r$PE1_pool_fallback - r$PE4_pool)[ok_pool], na.rm = TRUE) else NA_real_
  data.frame(
    block = z$block, n = z$n, p = z$p, beta1 = z$beta1, beta2 = z$beta2,
    pool_all = if (same) sp(e_pool, d_pool) else NA_real_,
    pool_nf  = if (same) sp(e_pool[ok_pool], d_pool[ok_pool]) else NA_real_,
    sep_all  = sp(e_sep, d_sep),
    sep_nf   = sp(e_sep[ok_sep], d_sep[ok_sep]),
    # number of replicates entering each correlation
    nrep_pool_all = sum(is.finite(e_pool) & is.finite(d_pool)),
    nrep_pool_nf  = sum(is.finite(e_pool[ok_pool]) & is.finite(d_pool[ok_pool])),
    nrep_sep_all  = sum(is.finite(e_sep) & is.finite(d_sep)),
    nrep_sep_nf   = sum(is.finite(e_sep[ok_sep]) & is.finite(d_sep[ok_sep])),
    share_fb_pool = mean(fb_pool, na.rm = TRUE),
    share_fb_sep  = mean(fb_sep,  na.rm = TRUE),
    rmse_all = if (same) sqrt(mean((r$bhat_pool - z$beta1)^2)) else NA_real_,
    rmse_nf  = if (same && any(ok_pool)) sqrt(mean((r$bhat_pool[ok_pool] - z$beta1)^2)) else NA_real_,
    cost_all = cost_all, cost_nf = cost_nf)
})
d <- do.call(rbind, rows)
d$setting <- ifelse(d$block == "A", paste0("A beta=", d$beta1),
                    paste0("B (", d$beta1, ",", d$beta2, ")"))

med <- function(v) if (all(is.na(v))) NA_real_ else median(v, na.rm = TRUE)
pos <- function(v) if (all(is.na(v))) NA_real_ else mean(v[!is.na(v)] > 0)

cat("\n=== PART 1. Within-cell Spearman (shape error vs paired difference), median over cells ===\n")
cat("columns: all = every replicate; nf = non-fallback replicates only; ncell = cells with a usable correlation\n\n")
for (rule in c("pool", "sep")) {
  cat("---", ifelse(rule == "pool", "COMMON covariance (PE-pooled vs beta=1)",
                    "CLASS-SPECIFIC (PE-separate vs beta=1)"), "---\n")
  a <- d[[paste0(rule, "_all")]]; b <- d[[paste0(rule, "_nf")]]
  fb <- d[[paste0("share_fb_", rule)]]
  nr_nf <- d[[paste0("nrep_", rule, "_nf")]]
  g <- data.frame(setting = d$setting, n = d$n, a = a, b = b, fb = fb, nr = nr_nf)
  out <- do.call(rbind, lapply(split(g, list(g$setting, g$n), drop = TRUE), function(s)
    data.frame(setting = s$setting[1], n = s$n[1],
               med_all = med(s$a), med_nf = med(s$b),
               ncell_all = sum(!is.na(s$a)), ncell_nf = sum(!is.na(s$b)),
               mean_fallback = mean(s$fb),
               med_nrep_nf = median(s$nr), min_nrep_nf = min(s$nr))))
  print(out[order(out$setting, out$n), ], row.names = FALSE)
  cat("share of cells with positive correlation: all =", round(pos(a), 2),
      "; non-fallback =", round(pos(b), 2), "\n")
  cat("cells with a usable correlation: all =", sum(!is.na(a)), "; non-fallback =", sum(!is.na(b)),
      "of", nrow(d), "cells (cells with fewer than", MIN_REP, "usable replicates are not counted)\n\n")
}

cat("=== PART 2. Cell-level Spearman (Block A): RMSE of pooled shape estimate vs cost of estimating the shape ===\n")
A <- d[d$block == "A", ]
res2 <- do.call(rbind, lapply(split(A, A$n), function(s) {
  data.frame(n = s$n[1], cells = nrow(s),
             rho_all = suppressWarnings(cor(s$rmse_all, s$cost_all, method = "spearman", use = "complete.obs")),
             rho_nf  = suppressWarnings(cor(s$rmse_nf,  s$cost_nf,  method = "spearman", use = "complete.obs")),
             cells_nf = sum(is.finite(s$rmse_nf) & is.finite(s$cost_nf)))
}))
print(res2, row.names = FALSE)
cat("(this is the PE-pooled, common-covariance version; the class-specific cell-level figures in the paper\n",
    " would need the same recomputation if you want them in the table as well)\n")

saveRDS(d, "fallback_cor_cells.rds")
cat("\nDone. Paste everything above this line.\n")
