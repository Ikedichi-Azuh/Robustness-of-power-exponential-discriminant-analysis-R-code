# ============================================================================
# analyze_robust_B.R : robust comparators in Block B (unequal PE shapes), 72 cells.
# Reads robust_results_B/robust_summary.csv (written by run_robust_B.R).
# ============================================================================
options(width = 170, digits = 4)
d <- read.csv("robust_results_B/robust_summary.csv", stringsAsFactors = FALSE)
stopifnot(nrow(d) == 72, all(d$repro_bayes_maxdiff == 0), all(d$repro_lda_maxdiff == 0))
TAU <- 0.5
d$shapes <- factor(paste0("(", d$beta1, ", ", d$beta2, ")"),
                   levels = c("(0.5, 1)", "(1, 2)", "(0.5, 2)"))
d$n <- factor(d$n)
cat("cells per group:\n"); print(table(d$shapes, d$n))

out <- list()
mean_tab <- function(cols, labs) {
  a <- aggregate(d[cols], by = list(shapes = d$shapes, n = d$n), FUN = mean)
  names(a)[-(1:2)] <- labs; a[order(a$shapes, a$n), ]
}
cat("\n=== 1a. Common covariance: mean excess error (pp) ===\n")
out$mean_pool <- mean_tab(c("excess_pp_PE1_pool_fallback", "excess_pp_LDA",
                            "excess_pp_MCD_pool", "excess_pp_T_pool"), c("PE", "LDA", "MCD", "t"))
print(out$mean_pool, row.names = FALSE)
cat("\n=== 1b. Class-specific covariance: mean excess error (pp) ===\n")
out$mean_sep <- mean_tab(c("excess_pp_PE1_sep_fallback", "excess_pp_QDA",
                           "excess_pp_MCD_sep", "excess_pp_T_sep"), c("PE", "QDA", "MCD", "t"))
print(out$mean_sep, row.names = FALSE)

cat("\n=== 2. Mean and worst-case (largest cell mean) excess error over the 72 cells ===\n")
ov <- function(cols) t(sapply(cols, function(cn) c(mean = mean(d[[cn]]), worst = max(d[[cn]]))))
cat("common covariance\n")
print(ov(c("excess_pp_PE1_pool_fallback", "excess_pp_LDA", "excess_pp_MCD_pool", "excess_pp_T_pool")))
cat("class-specific covariance\n")
print(ov(c("excess_pp_PE1_sep_fallback", "excess_pp_QDA", "excess_pp_MCD_sep", "excess_pp_T_sep")))

cls <- function(diff, mcse, tau = TAU) {
  lo <- diff - 1.96 * mcse; hi <- diff + 1.96 * mcse
  o <- ifelse(lo > -tau & hi < tau, "comparable",
       ifelse(hi < -tau, "first better", ifelse(lo > tau, "first worse", "inconclusive")))
  factor(o, levels = c("first better", "comparable", "first worse", "inconclusive"))
}
cmp <- c("PE_vs_MCD_pool", "PE_vs_MCD_sep", "PE_vs_T_pool", "PE_vs_T_sep",
         "MCD_vs_norm_pool", "MCD_vs_norm_sep", "T_vs_norm_pool", "T_vs_norm_sep")
cat("\n=== 3. Decision counts (first rule named first; negative favours the first rule) ===\n")
out$counts <- lapply(setNames(cmp, cmp), function(p) {
  cl <- cls(d[[paste0(p, "_diff_pp")]], d[[paste0(p, "_mcse_pp")]])
  tb <- table(d$shapes, cl); rbind(tb, All = colSums(tb))
})
for (nm in cmp) { cat("\n--", nm, "--\n"); print(out$counts[[nm]]) }

cat("\n=== 4. Mean paired difference PE minus robust rule (pp), by shapes and n ===\n")
diff_tab <- function(p) {
  a <- aggregate(d[paste0(p, "_diff_pp")], by = list(shapes = d$shapes, n = d$n), FUN = mean)
  names(a)[3] <- p; a[order(a$shapes, a$n), ]
}
out$diffs <- Reduce(function(x, y) merge(x, y, by = c("shapes", "n"), sort = FALSE),
                    lapply(cmp[1:4], diff_tab))
print(out$diffs[order(out$diffs$shapes, out$diffs$n), ], row.names = FALSE)

cat("\n=== 4b. Largest single-cell differences PE minus t (class-specific) and PE minus MCD ===\n")
for (p in c("PE_vs_T_sep", "PE_vs_T_pool", "PE_vs_MCD_sep")) {
  v <- d[[paste0(p, "_diff_pp")]]
  cat(p, ": min", round(min(v), 2), " max", round(max(v), 2), "\n")
}

cat("\n=== 5. Failure rates of the new rules ===\n")
fcols <- c("fail_MCD_pool", "fail_MCD_sep", "fail_T_pool", "fail_T_sep")
print(rbind(mean = colMeans(d[fcols]), max = sapply(d[fcols], max)))

cat("\n=== 6. Median chosen nu (pooled; class 1 and class 2 of the class-specific t rule) ===\n")
a <- aggregate(d[c("median_nu_pool", "median_nu1_sep", "median_nu2_sep")],
               by = list(shapes = d$shapes, n = d$n), FUN = mean)
out$nu <- a[order(a$shapes, a$n), ]; print(out$nu, row.names = FALSE)

saveRDS(out, "robust_results_B/robust_tables_B.rds")
cat("\nDone. Paste everything above this line.\n")
