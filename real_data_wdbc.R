# ============================================================================
# real_data_wdbc.R : real-data illustration with the Wisconsin diagnostic breast
# cancer data (569 patients: 357 benign, 212 malignant; dslabs::brca).
#
# The design is fixed here, before any result is looked at:
#   * variable sets (a priori, no selection on performance), all "mean" measurements
#     without perimeter and area, which are almost collinear with radius:
#       p = 2: radius, texture
#       p = 5: + smoothness, compactness, symmetry
#       p = 8: + concavity, concave points, fractal dimension
#   * training sizes per class: 25, 50, 100; test set: 100 per class, drawn from the
#     remaining patients of each class; equal priors (as in the simulation)
#   * B = 500 random stratified splits for each combination of p and n
#   * rules: PE with estimated shape (fallback policy, as in the simulation), LDA, QDA,
#     MCD-based Gaussian rule, multivariate-t rule (common and class-specific versions)
# Real data have no Bayes rule or true shape: only test errors are compared, no
# excess error. Dry run (B = 20):  options(pe.dry_run = TRUE); source("real_data_wdbc.R")
#                                  options(pe.dry_run = NULL)
# Needs: pe_discriminant_analysis_profile_optimized_6.R and robust_comparators.R in the
# working directory; packages dslabs, robustbase, MASS, parallel (ggplot2 for the figure).
# ============================================================================

suppressPackageStartupMessages({ library(MASS); library(parallel) })
source("pe_discriminant_analysis_profile_optimized_6.R")
source("robust_comparators.R")
if (!requireNamespace("dslabs", quietly = TRUE)) install.packages("dslabs")

DRY    <- isTRUE(getOption("pe.dry_run", FALSE))
B      <- if (DRY) 20L else 500L
NTR    <- c(25L, 50L, 100L)
NTE    <- 100L
SEED0  <- 20261106L
OUT    <- if (DRY) "real_data_wdbc_dry" else "real_data_wdbc"
CORES  <- max(1L, detectCores() - 1L)
dir.create(OUT, showWarnings = FALSE)

## ---- data --------------------------------------------------------------------------------
data("brca", package = "dslabs")
# The variables are standardised (z-scores over all patients, labels not used). All
# rules are affine-equivariant, so this changes no classification; it only avoids
# numerical problems from the very different scales of the measurements.
X0 <- scale(brca$x); y0 <- ifelse(brca$y == "B", 1L, 2L)
cat("patients:", nrow(X0), " benign:", sum(y0 == 1), " malignant:", sum(y0 == 2), "\n")
cn <- colnames(X0)
pick <- function(pat) { v <- cn[grep(pat, cn)]; stopifnot(length(v) >= 1); v[1] }
v_rad <- pick("^radius_mean$");     v_tex <- pick("^texture_mean$")
v_smo <- pick("^smoothness_mean$"); v_com <- pick("^compactness_mean$")
v_sym <- pick("^symmetry_mean$");   v_cav <- pick("^concavity_mean$")
v_cpt <- pick("^concave.*mean$");   v_fra <- pick("^fractal.*mean$")
VARS <- list(
  `2` = c(v_rad, v_tex),
  `5` = c(v_rad, v_tex, v_smo, v_com, v_sym),
  `8` = c(v_rad, v_tex, v_smo, v_com, v_sym, v_cav, v_cpt, v_fra))
cat("variables:\n"); print(VARS)

## ---- one random split -----------------------------------------------------------------------
raw_beta <- function(m, j = 1L) {
  b <- m$estimated_b
  if (is.null(b) || length(b) < j) NA_real_ else as.numeric(b[j])
}
split_once <- function(i, p, ntr) {
  set.seed(SEED0 + 1000L * p + 10L * ntr + i)
  vars <- VARS[[as.character(p)]]
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
  sep_acc  <- tryf(fit_pe_rule(Xtr, ytr, pooled = FALSE, est_beta = TRUE, boundary_policy = "accept"))
  sep_fb   <- tryf(fit_pe_rule(Xtr, ytr, pooled = FALSE, est_beta = TRUE, boundary_policy = "fallback"))
  yf  <- factor(ytr, levels = c(1, 2))
  lda_fit <- tryf(lda(Xtr, yf, prior = c(.5, .5)))
  qda_fit <- tryf(qda(Xtr, yf, prior = c(.5, .5)))
  err_class <- function(fit) if (is.null(fit)) NA_real_ else
    err_rate(as.integer(predict(fit, Xte)$class), yte)
  mcd <- tryf(fit_mcd(Xtr, ytr)); tt <- tryf(fit_t(Xtr, ytr))
  e_mcd <- function(S) if (is.null(mcd)) NA_real_ else
    tryCatch(err_rate(predict_gauss(mcd$m, mcd[[S]], Xte), yte), error = function(e) NA_real_)
  e_t <- function(k) if (is.null(tt)) NA_real_ else
    tryCatch(err_rate(predict_t(tt[[k]], Xte), yte), error = function(e) NA_real_)
  c(p = p, n = ntr, i = i,
    err_PE_pool = pe_err(pool_fb), err_PE_sep = pe_err(sep_fb),
    err_LDA = err_class(lda_fit), err_QDA = err_class(qda_fit),
    err_MCD_pool = e_mcd("S_pool"), err_MCD_sep = e_mcd("S_sep"),
    err_T_pool = e_t("pool"), err_T_sep = e_t("sep"),
    bhat_pool = if (is.null(pool_acc)) NA_real_ else raw_beta(pool_acc),
    bhat1_sep = if (is.null(sep_acc)) NA_real_ else raw_beta(sep_acc, 1L),
    bhat2_sep = if (is.null(sep_acc)) NA_real_ else raw_beta(sep_acc, 2L),
    fb_pool = if (is.null(pool_fb)) NA_real_ else as.numeric(isTRUE(pool_fb$fallback)),
    fb_sep  = if (is.null(sep_fb))  NA_real_ else as.numeric(isTRUE(sep_fb$fallback)),
    nu_pool = if (is.null(tt)) NA_real_ else tt$pool$nu[1])
}

## ---- run all splits ---------------------------------------------------------------------------
grid <- expand.grid(p = as.integer(names(VARS)), n = NTR)
t0 <- proc.time()[["elapsed"]]
sp <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  r <- mclapply(seq_len(B), split_once, p = grid$p[g], ntr = grid$n[g], mc.cores = CORES)
  do.call(rbind, r)
}))
sp <- as.data.frame(sp)
cat(sprintf("%d splits in %.1f minutes\n", nrow(sp), (proc.time()[["elapsed"]] - t0) / 60))
write.csv(sp, file.path(OUT, "splits.csv"), row.names = FALSE)

## ---- summaries ---------------------------------------------------------------------------------
mean_se <- function(v) { v <- v[is.finite(v)]; c(mean(v), sd(v) / sqrt(length(v))) }
paired <- function(a, b) { m <- mean_se(100 * (a - b)); c(diff = m[1], mcse = m[2]) }
rules <- c("PE_pool", "LDA", "MCD_pool", "T_pool", "PE_sep", "QDA", "MCD_sep", "T_sep")
tab <- do.call(rbind, lapply(split(sp, list(sp$p, sp$n)), function(s) {
  m <- sapply(rules, function(r) 100 * mean(s[[paste0("err_", r)]], na.rm = TRUE))
  data.frame(p = s$p[1], n = s$n[1], t(round(m, 2)),
             fail_PE_pool = mean(is.na(s$err_PE_pool)), fail_PE_sep = mean(is.na(s$err_PE_sep)))
}))
tab <- tab[order(tab$p, tab$n), ]
cat("\n=== Mean test error (%), 200 test patients per split, B =", B, "splits ===\n")
print(tab, row.names = FALSE)

cmp_defs <- list(
  "PE_pool - LDA"     = c("err_PE_pool", "err_LDA"),
  "PE_sep  - QDA"     = c("err_PE_sep",  "err_QDA"),
  "PE_pool - T_pool"  = c("err_PE_pool", "err_T_pool"),
  "PE_sep  - T_sep"   = c("err_PE_sep",  "err_T_sep"),
  "PE_pool - MCD_pool"= c("err_PE_pool", "err_MCD_pool"),
  "PE_sep  - MCD_sep" = c("err_PE_sep",  "err_MCD_sep"),
  "T_pool  - LDA"     = c("err_T_pool",  "err_LDA"),
  "T_sep   - QDA"     = c("err_T_sep",   "err_QDA"))
cmp <- do.call(rbind, lapply(split(sp, list(sp$p, sp$n)), function(s) {
  do.call(rbind, lapply(names(cmp_defs), function(nm) {
    z <- paired(s[[cmp_defs[[nm]][1]]], s[[cmp_defs[[nm]][2]]])
    data.frame(p = s$p[1], n = s$n[1], comparison = nm, diff_pp = z[["diff"]], mcse_pp = z[["mcse"]])
  }))
}))
cmp <- cmp[order(cmp$p, cmp$n, match(cmp$comparison, names(cmp_defs))), ]
cat("\n=== Paired differences in test error (percentage points; negative favours the first rule) ===\n")
print(cmp, row.names = FALSE)

shape <- do.call(rbind, lapply(split(sp, list(sp$p, sp$n)), function(s) {
  q <- function(v) quantile(v, c(.05, .5, .95), na.rm = TRUE)
  data.frame(p = s$p[1], n = s$n[1],
             bhat_pool_q05 = q(s$bhat_pool)[[1]], bhat_pool_med = q(s$bhat_pool)[[2]],
             bhat_pool_q95 = q(s$bhat_pool)[[3]],
             at_bound_pool = mean(s$bhat_pool <= 0.201 | s$bhat_pool >= 4.99, na.rm = TRUE),
             fallback_pool = mean(s$fb_pool, na.rm = TRUE),
             fallback_sep  = mean(s$fb_sep,  na.rm = TRUE),
             median_nu_pool = median(s$nu_pool, na.rm = TRUE))
}))
shape <- shape[order(shape$p, shape$n), ]
cat("\n=== Estimated shape (raw, before fallback), fallback rates, and chosen nu of the t rule ===\n")
print(shape, row.names = FALSE)

cat("\n=== PE fitted to ALL patients (descriptive): shape estimates ===\n")
for (p in names(VARS)) {
  f1 <- tryCatch(fit_pe_rule(X0[, VARS[[p]], drop = FALSE], y0, pooled = TRUE,  est_beta = TRUE,
                             boundary_policy = "accept"), error = function(e) NULL)
  f2 <- tryCatch(fit_pe_rule(X0[, VARS[[p]], drop = FALSE], y0, pooled = FALSE, est_beta = TRUE,
                             boundary_policy = "accept"), error = function(e) NULL)
  cat(sprintf("p = %s: pooled shape %.2f; class-specific shapes (benign, malignant) %.2f, %.2f\n",
              p, if (is.null(f1)) NA else raw_beta(f1),
              if (is.null(f2)) NA else raw_beta(f2, 1L), if (is.null(f2)) NA else raw_beta(f2, 2L)))
}
write.csv(tab,   file.path(OUT, "mean_errors.csv"),   row.names = FALSE)
write.csv(cmp,   file.path(OUT, "paired_differences.csv"), row.names = FALSE)
write.csv(shape, file.path(OUT, "shape_summary.csv"), row.names = FALSE)

cat("\nDone. Paste everything printed above, from 'patients:' to here.\n")
