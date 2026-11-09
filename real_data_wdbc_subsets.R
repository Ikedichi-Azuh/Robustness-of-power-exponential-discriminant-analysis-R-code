# ============================================================================
# real_data_wdbc_subsets.R : sensitivity of the real-data illustration to the choice
# of variables. The 8 "mean" measurements used in real_data_wdbc.R are the pool.
#   p = 2: ALL 28 pairs;   p = 5: ALL 56 five-variable subsets   (complete enumeration)
# For each variable set and training size (25, 50, 100 per class): B = 100 stratified
# random splits, 100 test patients per class, equal priors. Rules: PE (estimated shape,
# fallback policy), LDA, QDA, multivariate-t (common and class-specific). The MCD rule
# is left out (slow, and clearly worse in the main analysis).
# Output: how the paired differences vary across variable sets.
# Dry run:  options(pe.dry_run = TRUE); source("real_data_wdbc_subsets.R")
#           options(pe.dry_run = NULL)
# Needs the same files and packages as real_data_wdbc.R.
# ============================================================================

suppressPackageStartupMessages({ library(MASS); library(parallel) })
source("pe_discriminant_analysis_profile_optimized_6.R")
source("robust_comparators.R")

DRY   <- isTRUE(getOption("pe.dry_run", FALSE))
B     <- if (DRY) 5L else 100L
NTR   <- c(25L, 50L, 100L)
NTE   <- 100L
SEED0 <- 20261107L
OUT   <- if (DRY) "real_data_wdbc_subsets_dry" else "real_data_wdbc_subsets"
CORES <- max(1L, detectCores() - 1L)
TAU   <- 0.5
dir.create(OUT, showWarnings = FALSE)

## ---- data (standardised, as in real_data_wdbc.R) --------------------------------------
data("brca", package = "dslabs")
X0 <- scale(brca$x); y0 <- ifelse(brca$y == "B", 1L, 2L)
cn <- colnames(X0)
pick <- function(pat) { v <- cn[grep(pat, cn)]; stopifnot(length(v) >= 1); v[1] }
POOL <- c(pick("^radius_mean$"), pick("^texture_mean$"), pick("^smoothness_mean$"),
          pick("^compactness_mean$"), pick("^symmetry_mean$"), pick("^concavity_mean$"),
          pick("^concave.*mean$"), pick("^fractal.*mean$"))
cat("pool of variables:", paste(POOL, collapse = ", "), "\n")
sets <- c(combn(8, 2, simplify = FALSE), combn(8, 5, simplify = FALSE))
if (DRY) sets <- c(sets[1:3], sets[29:31])
cat("variable sets:", length(sets), " (p = 2:", sum(lengths(sets) == 2),
    ", p = 5:", sum(lengths(sets) == 5), ")\n")

## ---- one split ---------------------------------------------------------------------------
raw_beta <- function(m, j = 1L) {
  b <- m$estimated_b
  if (is.null(b) || length(b) < j) NA_real_ else as.numeric(b[j])
}
one <- function(vars, ntr, seed) {
  set.seed(seed)
  s1 <- sample(which(y0 == 1L), ntr + NTE); s2 <- sample(which(y0 == 2L), ntr + NTE)
  tr <- c(s1[1:ntr], s2[1:ntr]); te <- c(s1[(ntr + 1):(ntr + NTE)], s2[(ntr + 1):(ntr + NTE)])
  Xtr <- X0[tr, vars, drop = FALSE]; ytr <- y0[tr]
  Xte <- X0[te, vars, drop = FALSE]; yte <- y0[te]
  tryf <- function(expr) tryCatch(expr, error = function(e) NULL)
  pe_err <- function(m) {
    if (is.null(m) || isTRUE(m$fallback_failed) || is.null(m$mus) || is.null(m$S) || is.null(m$b))
      return(NA_real_)
    tryCatch(err_rate(predict_pe(m, Xte), yte), error = function(e) NA_real_)
  }
  pool_acc <- tryf(fit_pe_rule(Xtr, ytr, pooled = TRUE,  est_beta = TRUE, boundary_policy = "accept"))
  pool_fb  <- tryf(fit_pe_rule(Xtr, ytr, pooled = TRUE,  est_beta = TRUE, boundary_policy = "fallback"))
  sep_fb   <- tryf(fit_pe_rule(Xtr, ytr, pooled = FALSE, est_beta = TRUE, boundary_policy = "fallback"))
  yf <- factor(ytr, levels = c(1, 2))
  lda_fit <- tryf(lda(Xtr, yf, prior = c(.5, .5))); qda_fit <- tryf(qda(Xtr, yf, prior = c(.5, .5)))
  err_class <- function(fit) if (is.null(fit)) NA_real_ else
    err_rate(as.integer(predict(fit, Xte)$class), yte)
  tt <- tryf(fit_t(Xtr, ytr))
  e_t <- function(k) if (is.null(tt)) NA_real_ else
    tryCatch(err_rate(predict_t(tt[[k]], Xte), yte), error = function(e) NA_real_)
  c(PE_pool = pe_err(pool_fb), PE_sep = pe_err(sep_fb),
    LDA = err_class(lda_fit), QDA = err_class(qda_fit),
    T_pool = e_t("pool"), T_sep = e_t("sep"),
    bhat = if (is.null(pool_acc)) NA_real_ else raw_beta(pool_acc),
    fb_sep = if (is.null(sep_fb)) NA_real_ else as.numeric(isTRUE(sep_fb$fallback)))
}

## ---- one (variable set, n) task --------------------------------------------------------------
tasks <- expand.grid(set = seq_along(sets), n = NTR)
dm <- function(a, b) {
  v <- 100 * (a - b); v <- v[is.finite(v)]
  c(mean(v), sd(v) / sqrt(length(v)))
}
run_task <- function(k) {
  tk <- tasks[k, ]; idx <- sets[[tk$set]]; vars <- POOL[idx]
  M <- do.call(rbind, lapply(seq_len(B), function(i)
    one(vars, tk$n, SEED0 + 100000L * tk$set + 100L * tk$n + i)))
  z <- list(PE_LDA = dm(M[, "PE_pool"], M[, "LDA"]), PE_QDA = dm(M[, "PE_sep"], M[, "QDA"]),
            PE_Tpool = dm(M[, "PE_pool"], M[, "T_pool"]), PE_Tsep = dm(M[, "PE_sep"], M[, "T_sep"]))
  data.frame(set = tk$set, p = length(idx), n = tk$n, vars = paste(idx, collapse = "-"),
             err_PE_pool = 100 * mean(M[, "PE_pool"], na.rm = TRUE),
             err_LDA = 100 * mean(M[, "LDA"], na.rm = TRUE),
             err_PE_sep = 100 * mean(M[, "PE_sep"], na.rm = TRUE),
             err_QDA = 100 * mean(M[, "QDA"], na.rm = TRUE),
             d_PE_LDA = z$PE_LDA[1], se_PE_LDA = z$PE_LDA[2],
             d_PE_QDA = z$PE_QDA[1], se_PE_QDA = z$PE_QDA[2],
             d_PE_Tpool = z$PE_Tpool[1], se_PE_Tpool = z$PE_Tpool[2],
             d_PE_Tsep = z$PE_Tsep[1], se_PE_Tsep = z$PE_Tsep[2],
             bhat_med = median(M[, "bhat"], na.rm = TRUE),
             fallback_sep = mean(M[, "fb_sep"], na.rm = TRUE))
}
t0 <- proc.time()[["elapsed"]]
res <- do.call(rbind, mclapply(seq_len(nrow(tasks)), run_task, mc.cores = CORES))
cat(sprintf("%d variable-set x n tasks (%d splits) in %.1f minutes\n",
            nrow(res), nrow(res) * B, (proc.time()[["elapsed"]] - t0) / 60))
write.csv(res, file.path(OUT, "set_level_results.csv"), row.names = FALSE)

## ---- summary across variable sets ------------------------------------------------------------------
comp <- list(`PE common - LDA` = c("d_PE_LDA", "se_PE_LDA"),
             `PE class-specific - QDA` = c("d_PE_QDA", "se_PE_QDA"),
             `PE common - t common` = c("d_PE_Tpool", "se_PE_Tpool"),
             `PE class-specific - t class-specific` = c("d_PE_Tsep", "se_PE_Tsep"))
summ <- do.call(rbind, lapply(split(res, list(res$p, res$n)), function(s) {
  do.call(rbind, lapply(names(comp), function(nm) {
    d <- s[[comp[[nm]][1]]]; se <- s[[comp[[nm]][2]]]
    data.frame(p = s$p[1], n = s$n[1], comparison = nm, sets = nrow(s),
               median_diff = median(d), min_diff = min(d), max_diff = max(d),
               share_PE_lower = mean(d < 0),
               share_clearly_lower = mean(d + 1.96 * se < 0),
               share_clearly_higher = mean(d - 1.96 * se > 0),
               share_lower_by_margin = mean(d < -TAU))
  }))
}))
summ <- summ[order(summ$p, summ$n, match(summ$comparison, names(comp))), ]
cat("\n=== Paired difference in test error (pp; negative favours PE) across variable sets ===\n")
cat("share_PE_lower: difference below 0; clearly: 95% interval excludes 0;",
    "lower_by_margin: PE lower by more than 0.5 pp\n")
print(summ, row.names = FALSE, digits = 3)
write.csv(summ, file.path(OUT, "summary_across_sets.csv"), row.names = FALSE)

sh <- do.call(rbind, lapply(split(res, list(res$p, res$n)), function(s)
  data.frame(p = s$p[1], n = s$n[1], median_of_median_bhat = median(s$bhat_med),
             min_bhat = min(s$bhat_med), max_bhat = max(s$bhat_med),
             median_fallback_sep = median(s$fallback_sep),
             median_error_LDA = median(s$err_LDA), min_error_LDA = min(s$err_LDA),
             max_error_LDA = max(s$err_LDA))))
sh <- sh[order(sh$p, sh$n), ]
cat("\n=== Estimated shape (median over splits) and error level across variable sets ===\n")
print(sh, row.names = FALSE, digits = 3)

## ---- figure -------------------------------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  f <- rbind(data.frame(res[c("p", "n")], comparison = "Common covariance (PE - LDA)", diff = res$d_PE_LDA),
             data.frame(res[c("p", "n")], comparison = "Class-specific (PE - QDA)", diff = res$d_PE_QDA))
  f$p <- factor(paste("p =", f$p)); f$n <- factor(paste("n =", f$n), levels = paste("n =", NTR))
  g <- ggplot(f, aes(diff, comparison)) +
    geom_vline(xintercept = 0, colour = "#52514e", linewidth = 0.3) +
    geom_vline(xintercept = c(-TAU, TAU), colour = "#52514e", linewidth = 0.3, linetype = "dashed") +
    geom_point(colour = "#2a78d6", alpha = 0.6, size = 1.3,
               position = position_jitter(height = 0.18, width = 0, seed = 1)) +
    facet_grid(p ~ n) +
    labs(x = "Difference in test error across variable sets (pp; negative: PE better)", y = NULL) +
    theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank())
  dir.create("figures", showWarnings = FALSE)
  nm <- if (DRY) "fig15_wdbc_subsets_dry" else "fig15_wdbc_subsets"
  ggsave(file.path("figures", paste0(nm, ".pdf")), g, width = 7.5, height = 3.8)
  ggsave(file.path("figures", paste0(nm, ".png")), g, width = 7.5, height = 3.8, dpi = 300)
}
cat("\nDone. Paste everything printed above.\n")
