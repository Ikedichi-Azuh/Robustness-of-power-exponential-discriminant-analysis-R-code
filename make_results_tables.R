# ============================================================================
# make_results_tables.R : prints Tables 3-12 of the paper (Results, Blocks A, B, C)
# from the cell-level summaries, in the layout of the paper.
#
# Reads   grid_final/cell_summary.csv     (Blocks A and B, 288 cells)
#         grid_final_C/cell_summary.csv   (Block C, 144 cells)
# Writes  results_tables/table03.csv ... table12.csv and results_tables/tables.rds
#         (the tables are also printed to the console)
#
# Table 3  combinations with a fallback rate above 5%
# Table 4  Block A: bias, RMSE and boundary share of the pooled shape estimate
# Table 5  Block A, equal covariance: mean excess error of the common-covariance rules
# Table 6  decision map, approach 1 vs approach 3 (class-specific), Blocks A and B
# Table 7  Block A: common minus class-specific version of each approach
# Table 8  Block A: cost of estimating the shape (approach 1 minus approach 4)
# Table 9  Block C: pooled shape estimate
# Table 10 Block C: mean excess error by distribution
# Table 11 Block C: decision map by distribution
# Table 12 Block C: mean and worst-case excess error over the 144 cells
#
# Decision rule: 95% interval (diff +/- 1.96 MCSE) against tau = 0.5 pp.
# Group means are averages of cell means over the factors not shown.
# ============================================================================

IN_AB <- "grid_final/cell_summary.csv"
IN_C  <- "grid_final_C/cell_summary.csv"
OUT   <- "results_tables"
TAU   <- 0.5
dir.create(OUT, showWarnings = FALSE)

ab <- read.csv(IN_AB, stringsAsFactors = FALSE)
cc <- read.csv(IN_C,  stringsAsFactors = FALSE)
stopifnot(nrow(ab) == 288, nrow(cc) == 144)
A <- ab[ab$block == "A", ]
B <- ab[ab$block == "B", ]
stopifnot(nrow(A) == 216, nrow(B) == 72)

tables <- list()
show <- function(id, title, df, digits = 3) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
  print(df, row.names = FALSE, digits = digits)
  write.csv(df, file.path(OUT, sprintf("table%02d.csv", id)), row.names = FALSE)
  tables[[sprintf("table%02d", id)]] <<- df
}

# Class of one cell: better / comparable / worse / inconclusive (approach 1 vs comparator)
cls <- function(diff, mcse, tau = TAU) {
  lo <- diff - 1.96 * mcse; hi <- diff + 1.96 * mcse
  out <- ifelse(lo > -tau & hi < tau, "Comparable",
         ifelse(hi < -tau, "Better",
         ifelse(lo > tau, "Worse", "Inconclusive")))
  factor(out, levels = c("Better", "Comparable", "Worse", "Inconclusive"))
}

# ---- Table 3: fallback above 5% -------------------------------------------------
ab$bb <- paste0("(", ab$beta1, ", ", ab$beta2, ")")
t3 <- aggregate(cbind(Common = fallback_pool, Class = fallback_sep) ~
                  block + bb + n + p, data = ab, FUN = mean)
n_combo <- nrow(t3)
t3 <- t3[t3$Common > 0.05 | t3$Class > 0.05, ]
t3 <- t3[order(t3$block, t3$bb, t3$n, t3$p), ]
names(t3) <- c("Block", "Shapes", "n", "p", "Common", "Class-specific")
cat("\nCombinations of block, n, p and shapes:", n_combo,
    "; with fallback above 5%:", nrow(t3), "\n")
show(3, "Table 3. Combinations with a fallback rate above 5% (proportion of replicates)",
     t3)

# ---- Table 4: pooled shape estimate, Block A ----------------------------------
t4 <- aggregate(cbind(bias_pool, rmse_pool, boundary_pool) ~ beta1 + p + n,
                data = A[A$n %in% c(50, 500), ], FUN = mean)
w4 <- reshape(t4, idvar = c("beta1", "p"), timevar = "n", direction = "wide")
w4 <- w4[order(w4$beta1, w4$p), ]
names(w4) <- sub("beta1", "beta", names(w4))
show(4, "Table 4. Block A: bias, RMSE and boundary share of the pooled shape estimate",
     w4)

# ---- Table 5: common-covariance rules, Block A, equal covariance ----------------
E <- A[A$cov_type == "equal", ]
t5 <- aggregate(cbind(`PE (est.)` = excess_pp_PE1_pool_fallback,
                      LDA = excess_pp_LDA,
                      `PE (true)` = excess_pp_PE4_pool) ~ beta1 + n,
                data = E, FUN = mean)
names(t5)[1] <- "beta"
show(5, "Table 5. Block A, equal covariance: mean excess error (pp), common-covariance rules",
     t5)
cat("Check: max |PE2_pool - LDA| in these group means =",
    format(max(abs(tapply(E$excess_pp_PE2_pool, list(E$beta1, E$n), mean) -
                   tapply(E$excess_pp_LDA, list(E$beta1, E$n), mean))),
           digits = 3), "pp\n")

# ---- Table 6: decision map, approach 1 vs approach 3 (class-specific) -----------
d6 <- ab
d6$class <- cls(d6$RQ2_sep_diff_pp, d6$RQ2_sep_mcse_pp)
t6 <- as.data.frame.matrix(table(paste(d6$block, d6$cov_type), d6$class))
t6 <- cbind(Block = sub(" .*", "", rownames(t6)),
            Covariance = sub(".* ", "", rownames(t6)), t6)
names(t6)[3:6] <- c("Approach 1 better", "Comparable", "Approach 1 worse", "Inconclusive")
show(6, "Table 6. Decision map, approach 1 vs approach 3 (class-specific), number of cells",
     t6)
cl2 <- cls(ab$RQ1_sep_diff_pp, ab$RQ1_sep_mcse_pp)
cat("Approach 1 vs approach 2 (class-specific), all 288 cells:",
    paste(names(table(cl2)), table(cl2), collapse = "; "), "\n")
cat("Common covariance, 288 cells, approach 1 vs 2 / vs 3:",
    paste(table(cls(ab$RQ1_pool_diff_pp, ab$RQ1_pool_mcse_pp)), collapse = "/"), " ; ",
    paste(table(cls(ab$RQ2_pool_diff_pp, ab$RQ2_pool_mcse_pp)), collapse = "/"),
    "(order: better, comparable, worse, inconclusive)\n")

# ---- Table 7: common minus class-specific, Block A ---------------------------
t7 <- aggregate(cbind(`PE (est.)` = cov_PE1_diff_pp,
                      `PE (beta=1)` = cov_PE2_diff_pp,
                      `LDA vs QDA` = cov_normal_diff_pp) ~ cov_type + n,
                data = A, FUN = mean)
t7 <- t7[order(t7$cov_type, t7$n), ]
names(t7)[1] <- "Covariance"
show(7, "Table 7. Block A: mean paired difference (pp), common minus class-specific", t7)

# ---- Table 8: cost of estimating the shape, Block A --------------------------
t8 <- aggregate(cbind(Com = shape_est_pool_diff_pp, Class = shape_est_sep_diff_pp) ~
                  beta1 + n, data = A, FUN = mean)
w8 <- reshape(t8, idvar = "beta1", timevar = "n", direction = "wide")
names(w8)[1] <- "beta"
show(8, "Table 8. Block A: cost of estimating the shape (approach 1 minus approach 4), pp",
     w8)

# ---- Block C labels ---------------------------------------------------------
dist_lv  <- c("Student t, 3 d.f.", "Contaminated normal", "Skew-t, 5 d.f.",
              "Student t, 5 d.f.", "Student t, 10 d.f.", "Skew-normal")
dist_map <- c(t3 = dist_lv[1], cn = dist_lv[2], st = dist_lv[3],
              t5 = dist_lv[4], t10 = dist_lv[5], sn = dist_lv[6])
cc$Distribution <- factor(unname(dist_map[cc$dist]), levels = dist_lv)
stopifnot(!anyNA(cc$Distribution))

# ---- Table 9: Block C shape estimate -----------------------------------------
cc$at_bound_pct <- 100 * cc$lower_pool
t9 <- aggregate(cbind(Median = median_bhat_pool, P05 = q05_bhat_pool,
                      P95 = q95_bhat_pool, `At bound (%)` = at_bound_pct) ~
                  Distribution + n, data = cc, FUN = mean)
w9 <- reshape(t9, idvar = "Distribution", timevar = "n", direction = "wide")
show(9, "Table 9. Block C: pooled shape estimate (median, 5th and 95th percentile, % at lower bound 0.2)",
     w9, digits = 2)

# ---- Table 10: Block C mean excess error -------------------------------------
t10 <- aggregate(cbind(PE_common = excess_pp_PE1_pool_fallback, LDA = excess_pp_LDA,
                       PE_class = excess_pp_PE1_sep_fallback, QDA = excess_pp_QDA) ~
                   Distribution + n, data = cc, FUN = mean)
w10 <- reshape(t10, idvar = "Distribution", timevar = "n", direction = "wide")
show(10, "Table 10. Block C: mean excess error (pp), fallback policy", w10, digits = 3)

# ---- Table 11: Block C decision map by distribution --------------------------
cc$cls_pool <- cls(cc$RQ2_pool_diff_pp, cc$RQ2_pool_mcse_pp)
cc$cls_sep  <- cls(cc$RQ2_sep_diff_pp,  cc$RQ2_sep_mcse_pp)
tp <- as.data.frame.matrix(table(cc$Distribution, cc$cls_pool))
ts <- as.data.frame.matrix(table(cc$Distribution, cc$cls_sep))
names(tp) <- paste0("Common ", names(tp))
names(ts) <- paste0("Class ", names(ts))
t11 <- cbind(Distribution = rownames(tp), tp, ts)
show(11, "Table 11. Block C: number of the 24 cells per distribution (approach 1 vs LDA / QDA)",
     t11)
cat("Totals, common (vs LDA):", paste(colSums(tp), collapse = "/"),
    " class-specific (vs QDA):", paste(colSums(ts), collapse = "/"),
    "(better/comparable/worse/inconclusive)\n")
for (v in c("accept")) {
  cp <- cls(cc$RQ2_pool_accept_diff_pp, cc$RQ2_pool_accept_mcse_pp)
  cs <- cls(cc$RQ2_sep_accept_diff_pp,  cc$RQ2_sep_accept_mcse_pp)
  cat("Boundary estimates accepted: common", paste(table(cp), collapse = "/"),
      "; class-specific", paste(table(cs), collapse = "/"), "\n")
}

# ---- Table 12: mean and worst case over the 144 cells ------------------------
rl <- c("PE, estimated shape, common (fallback)"       = "excess_pp_PE1_pool_fallback",
        "PE, estimated shape, common (boundary accepted)" = "excess_pp_PE1_pool_accept",
        "LDA"                                          = "excess_pp_LDA",
        "PE, estimated shape, class-specific (fallback)" = "excess_pp_PE1_sep_fallback",
        "PE, estimated shape, class-specific (boundary accepted)" = "excess_pp_PE1_sep_accept",
        "QDA"                                          = "excess_pp_QDA")
t12 <- data.frame(Rule = names(rl),
                  Mean = sapply(rl, function(v) mean(cc[[v]])),
                  `Worst case` = sapply(rl, function(v) max(cc[[v]])),
                  check.names = FALSE, row.names = NULL)
show(12, "Table 12. Block C: mean and worst-case excess error (pp) over the 144 cells",
     t12, digits = 3)

saveRDS(tables, file.path(OUT, "tables.rds"))
cat("\nSaved", length(tables), "tables to", OUT, "/\n")
