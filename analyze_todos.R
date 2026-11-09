# ============================================================================
# analyze_todos.R : the numbers needed to close the remaining \todo{} items in
# the Results of Blocks A and B. Reads grid_final/cell_summary.csv and the stored
# replicate files grid_final/cells/*B1000.rds. Prints four parts; paste everything.
#   1. sensitivity analyses: boundary estimates accepted; replicates with any
#      failed fit excluded (compared with the primary fallback analysis)
#   2. accuracy of the class-specific shape estimates (Blocks A and B)
#   3. largest Monte Carlo standard errors of the paired differences (adequacy of B)
#   4. cell-level count: class-specific against common covariance (unequal covariance)
# ============================================================================
options(width = 160, digits = 4)
TAU <- 0.5
d <- read.csv("grid_final/cell_summary.csv", stringsAsFactors = FALSE)
cat("cells in cell_summary.csv:", nrow(d), "\n\n")

## ---- 2. class-specific shape accuracy ------------------------------------------------
cat("=== 2. Shape estimation, class-specific fits (mean over cells) ===\n")
cols <- intersect(c("bias_pool", "rmse_pool", "bias1_sep", "bias2_sep", "rmse1_sep", "rmse2_sep"),
                  names(d))
cat("columns found:", paste(cols, collapse = ", "), "\n")
cat("all bias/rmse columns present:", paste(grep("bias|rmse", names(d), value = TRUE), collapse = ", "), "\n")
print(aggregate(d[cols], by = list(block = d$block, beta1 = d$beta1, beta2 = d$beta2, n = d$n),
                FUN = mean), row.names = FALSE)

## ---- 3. largest MCSE --------------------------------------------------------------------
cat("\n=== 3. Largest Monte Carlo standard error of the paired differences (pp), by block ===\n")
mc <- grep("mcse_pp", names(d), value = TRUE)
print(t(sapply(split(d[mc], d$block), function(x) sapply(x, max, na.rm = TRUE))))
cat("overall maximum:", max(unlist(d[mc]), na.rm = TRUE), "\n")

## ---- 1 and 4: replicate-level quantities ----------------------------------------------------
files <- list.files("grid_final/cells", pattern = "B1000\\.rds$", full.names = TRUE)
cat("\nreading", length(files), "cell files ...\n")
ex <- function(a, b) {
  v <- 100 * (a - b); v <- v[is.finite(v)]
  if (!length(v)) return(c(NA_real_, NA_real_))
  c(mean(v), sd(v) / sqrt(length(v)))
}
res <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f); z <- x$spec; r <- x$res
  fail <- r$fallback_pool_fallback == 1 | r$fallback_sep_fallback == 1
  keep <- !is.na(fail) & !fail
  e <- function(a, b, i = TRUE) ex(a[i], b[i])[1]
  cov_sep <- ex(r$PE1_sep_fallback, r$PE1_pool_fallback)
  data.frame(
    block = z$block, n = z$n, p = z$p, beta1 = z$beta1, beta2 = z$beta2,
    cov = z$cov_type, sep = z$separation,
    pool_prim = e(r$PE1_pool_fallback, r$LDA),
    pool_acc  = e(r$PE1_pool_accept,   r$LDA),
    pool_excl = e(r$PE1_pool_fallback, r$LDA, keep),
    sep_prim  = e(r$PE1_sep_fallback,  r$QDA),
    sep_acc   = e(r$PE1_sep_accept,    r$QDA),
    sep_excl  = e(r$PE1_sep_fallback,  r$QDA, keep),
    share_excl = mean(!keep),
    covdiff = cov_sep[1], covmcse = cov_sep[2])
}))

cat("\n=== 1. Sensitivity analyses: mean paired difference PE (estimated shape) minus normal rule (pp) ===\n")
cat("(pool = common covariance vs LDA; sep = class-specific vs QDA; prim = fallback policy;",
    "acc = boundary estimates accepted; excl = replicates with any fallback removed)\n")
sens <- aggregate(cbind(pool_prim, pool_acc, pool_excl, sep_prim, sep_acc, sep_excl, share_excl) ~
                    block + n, data = res, FUN = mean)
print(sens, row.names = FALSE)
cat("\nlargest |change| across cells (pp):\n")
print(c(pool_accept = max(abs(res$pool_acc - res$pool_prim), na.rm = TRUE),
        pool_excl   = max(abs(res$pool_excl - res$pool_prim), na.rm = TRUE),
        sep_accept  = max(abs(res$sep_acc - res$sep_prim), na.rm = TRUE),
        sep_excl    = max(abs(res$sep_excl - res$sep_prim), na.rm = TRUE)))
cat("cells where every replicate had a fallback (excl undefined):", sum(is.nan(res$pool_excl) | is.nan(res$sep_excl)), "\n")

cat("\n=== 4. Class-specific vs common covariance (approach 1), unequal covariance, Block A ===\n")
sub <- res[res$cov == "unequal" & res$block == "A", ]
lo <- sub$covdiff - 1.96 * sub$covmcse; hi <- sub$covdiff + 1.96 * sub$covmcse
sub$verdict <- ifelse(lo > -TAU & hi < TAU, "comparable",
               ifelse(hi < -TAU, "class-specific better",
               ifelse(lo > TAU, "common better", "inconclusive")))
print(table(n = sub$n, verdict = sub$verdict))
cat("\nDone. Paste everything above this line.\n")
