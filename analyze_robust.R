# ============================================================================
# analyze_robust.R : summaries of the robust-comparator study
# Reads robust_results/robust_summary.csv (216 cells: 72 Block A, 144 Block C)
# and prints the tables needed for the Results. Nothing is written except
# robust_results/robust_tables.rds (all tables, for later use).
# ============================================================================
options(width = 160, digits = 4)
d <- read.csv("robust_results/robust_summary.csv", stringsAsFactors = FALSE)
stopifnot(nrow(d) == 216, all(d$repro_bayes_maxdiff == 0), all(d$repro_lda_maxdiff == 0))

TAU <- 0.5
d$group <- ifelse(d$source == "C", d$dist, paste0("A_beta", d$beta1))
d$grp_lab <- factor(d$group,
  levels = c("A_beta0.5", "A_beta1", "A_beta2", "t3", "cn", "st", "t5", "t10", "sn"))
d$n <- factor(d$n)
cat("cells per group:\n"); print(table(d$grp_lab, d$n))

out <- list()

## ---- 1. mean excess error (pp) by group and n -----------------------------------
mean_tab <- function(cols, labs) {
  a <- aggregate(d[cols], by = list(group = d$grp_lab, n = d$n), FUN = mean)
  names(a)[-(1:2)] <- labs
  a[order(a$group, a$n), ]
}
cat("\n=== 1a. Common covariance: mean excess error (pp) ===\n")
out$mean_pool <- mean_tab(
  c("excess_pp_PE1_pool_fallback", "excess_pp_LDA", "excess_pp_MCD_pool", "excess_pp_T_pool"),
  c("PE", "LDA", "MCD", "t"))
print(out$mean_pool, row.names = FALSE)
cat("\n=== 1b. Class-specific covariance: mean excess error (pp) ===\n")
out$mean_sep <- mean_tab(
  c("excess_pp_PE1_sep_fallback", "excess_pp_QDA", "excess_pp_MCD_sep", "excess_pp_T_sep"),
  c("PE", "QDA", "MCD", "t"))
print(out$mean_sep, row.names = FALSE)

## ---- 2. overall mean and worst case ---------------------------------------------------
cat("\n=== 2. Mean and worst-case (largest cell mean) excess error, by block ===\n")
ov <- function(sub, cols, labs) {
  t(sapply(cols, function(cn) c(mean = mean(sub[[cn]]), worst = max(sub[[cn]]))))
}
for (bl in c("A", "C")) {
  sub <- d[d$source == ifelse(bl == "A", "AB", "C"), ]
  cat("\nBlock", bl, "(", nrow(sub), "cells )\n")
  cat("common covariance\n")
  print(ov(sub, c("excess_pp_PE1_pool_fallback", "excess_pp_LDA",
                  "excess_pp_MCD_pool", "excess_pp_T_pool")))
  cat("class-specific covariance\n")
  print(ov(sub, c("excess_pp_PE1_sep_fallback", "excess_pp_QDA",
                  "excess_pp_MCD_sep", "excess_pp_T_sep")))
}

## ---- 3. decision counts (95% interval against tau = 0.5 pp) --------------------------
cls <- function(diff, mcse, tau = TAU) {
  lo <- diff - 1.96 * mcse; hi <- diff + 1.96 * mcse
  out <- ifelse(lo > -tau & hi < tau, "comparable",
         ifelse(hi < -tau, "first better",
         ifelse(lo > tau, "first worse", "inconclusive")))
  factor(out, levels = c("first better", "comparable", "first worse", "inconclusive"))
}
count_tab <- function(prefix) {
  cl <- cls(d[[paste0(prefix, "_diff_pp")]], d[[paste0(prefix, "_mcse_pp")]])
  tb <- table(d$grp_lab, cl)
  rbind(tb, All = colSums(tb))
}
cmp <- c(
  "PE_vs_MCD_pool", "PE_vs_MCD_sep", "PE_vs_T_pool", "PE_vs_T_sep",   # first = PE
  "MCD_vs_norm_pool", "MCD_vs_norm_sep", "T_vs_norm_pool", "T_vs_norm_sep")  # first = robust
cat("\n=== 3. Decision counts per comparison (first rule named first) ===\n")
cat("(PE_vs_X: PE minus X; X_vs_norm: X minus LDA/QDA; negative favours the first rule)\n")
out$counts <- lapply(setNames(cmp, cmp), count_tab)
for (nm in cmp) { cat("\n--", nm, "--\n"); print(out$counts[[nm]]) }

## ---- 4. mean paired differences (pp) by group and n -------------------------------------
cat("\n=== 4. Mean paired difference PE minus robust rule (pp), by group and n ===\n")
diff_tab <- function(prefix) {
  a <- aggregate(d[paste0(prefix, "_diff_pp")], by = list(group = d$grp_lab, n = d$n), FUN = mean)
  names(a)[3] <- prefix; a[order(a$group, a$n), ]
}
out$diffs <- Reduce(function(x, y) merge(x, y, by = c("group", "n"), sort = FALSE),
                    lapply(cmp[1:4], diff_tab))
print(out$diffs[order(out$diffs$group, out$diffs$n), ], row.names = FALSE)

## ---- 5. failures of the new rules ---------------------------------------------------------
cat("\n=== 5. Failure rates of the new rules (share of replicates without a prediction) ===\n")
fcols <- c("fail_MCD_pool", "fail_MCD_sep", "fail_T_pool", "fail_T_sep")
out$fail <- rbind(mean = colMeans(d[fcols]), max = sapply(d[fcols], max))
print(out$fail)

## ---- 6. degrees of freedom chosen by the t rule ----------------------------------------------
cat("\n=== 6. Median chosen nu (pooled; class-specific class 1 and 2), mean of cell medians ===\n")
a <- aggregate(d[c("median_nu_pool", "median_nu1_sep", "median_nu2_sep")],
               by = list(group = d$grp_lab, n = d$n), FUN = mean)
out$nu <- a[order(a$group, a$n), ]
print(out$nu, row.names = FALSE)

saveRDS(out, "robust_results/robust_tables.rds")
cat("\nDone. Paste everything above this line.\n")
