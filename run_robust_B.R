# ============================================================================
# run_robust_B.R : robust comparators (MCD, multivariate-t) for the 72 Block B cells
# (unequal PE shapes; n in {50, 250}, p in {2, 5}), stored in grid_final.
# Separate from run_robust.R: writes only to robust_results_B/ and leaves
# robust_results/ untouched.
#
# For every Block B cell the script regenerates the SAME training and test sets (same
# seed, same draw order), fits the new rules, and compares them with the stored PE, LDA and
# QDA results. Each replicate is checked: the regenerated Bayes error and LDA error must
# equal the stored ones exactly.
#
# Dry run (5 replicates, 2 cells):
#   options(pe.dry_run = TRUE); source("run_robust_B.R"); options(pe.dry_run = NULL)
# Restartable: finished cells are skipped.
# Output: robust_results_B/ (per-cell RDS) and robust_results_B/robust_summary.csv
# Needs: pe_discriminant_analysis_profile_optimized_6.R, robust_comparators.R, grid_final/
# ============================================================================

suppressPackageStartupMessages(library(parallel))

SCRIPT_MAIN <- "pe_discriminant_analysis_profile_optimized_6.R"
SCRIPT_ROB  <- "robust_comparators.R"
DIR_AB      <- "grid_final"
OUT_DIR     <- "robust_results_B"
N_WORKERS   <- max(1L, detectCores() - 1L)
VAR2 <- 1.25; RHO2 <- 0.5
DRY_RUN <- isTRUE(getOption("pe.dry_run", FALSE))
if (DRY_RUN) { OUT_DIR <- "robust_results_B_dry"; cat("*** DRY RUN: 5 replicates, 2 cells ***\n") }

SCRIPT_MAIN <- normalizePath(SCRIPT_MAIN, mustWork = TRUE)
SCRIPT_ROB  <- normalizePath(SCRIPT_ROB,  mustWork = TRUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- normalizePath(OUT_DIR)
CONFIG <- list(main = unname(tools::md5sum(SCRIPT_MAIN)),
               rob = unname(tools::md5sum(SCRIPT_ROB)),
               VAR2 = VAR2, RHO2 = RHO2, dry = DRY_RUN)

## ---- task list: the stored Block B cells --------------------------------------------
fs <- list.files(file.path(DIR_AB, "cells"), pattern = "B1000\\.rds$", full.names = TRUE)
if (!length(fs)) stop("No stored cells found in ", file.path(DIR_AB, "cells"))
tasks <- do.call(rbind, lapply(fs, function(f) {
  z <- readRDS(f)$spec
  data.frame(file = f, block = z$block, n = z$n, p = z$p, cell_id = z$cell_id,
             stringsAsFactors = FALSE)
}))
tasks <- tasks[tasks$block == "B", ]
tasks$source <- "B"
if (!DRY_RUN && nrow(tasks) != 72) stop("expected 72 Block B cells, found ", nrow(tasks))
if (DRY_RUN) tasks <- head(tasks, 2)
tasks$out <- file.path(OUT_DIR, sprintf("B_cell_%04d.rds", tasks$cell_id))
stopifnot(!anyDuplicated(tasks$out))
cat("Block B cells to process:", nrow(tasks), "\n")

## ---- worker function ------------------------------------------------------------------
B_DRY <- 5L
run_task <- function(task) {
  if (file.exists(task$out)) {
    ok <- tryCatch({ x <- readRDS(task$out); identical(x$config, CONFIG) },
                   error = function(e) FALSE)
    if (ok) return(data.frame(out = task$out, status = "cached", seconds = NA_real_))
  }
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    x <- readRDS(task$file); z <- x$spec; stored <- x$res
    B <- if (DRY_RUN) B_DRY else nrow(stored)
    s <- envW$make_setup(z$p, z$beta1, z$beta2, z$cov_type, var2 = VAR2, rho2 = RHO2)
    set.seed(z$seed)
    res <- as.data.frame(t(replicate(B, envW$robust_rep(s, z$n, z$theta))))
    saveRDS(list(config = CONFIG, spec = z, source = "B", res = res,
                 stored_rows = seq_len(B)), task$out)
    data.frame(out = task$out, status = "done", seconds = proc.time()[["elapsed"]] - t0)
  }, error = function(e) {
    data.frame(out = task$out, status = paste("error:", conditionMessage(e)),
               seconds = NA_real_)
  })
}

cl <- makeCluster(N_WORKERS)
clusterExport(cl, c("SCRIPT_MAIN", "SCRIPT_ROB", "OUT_DIR", "CONFIG", "VAR2", "RHO2",
                    "DRY_RUN", "B_DRY", "run_task"), envir = environment())
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(MASS))
  envW <- new.env()
  sys.source(SCRIPT_MAIN, envir = envW)
  sys.source(SCRIPT_ROB,  envir = envW)
}))

t_start <- proc.time()[["elapsed"]]
run_log <- do.call(rbind, parLapplyLB(cl, lapply(seq_len(nrow(tasks)), function(i) tasks[i, ]),
                                  run_task))
stopCluster(cl)
cat("Status:", paste(names(table(run_log$status)), as.integer(table(run_log$status)), collapse = "; "), "\n")
cat(sprintf("Wall-clock time this session: %.1f minutes on %d workers\n",
            (proc.time()[["elapsed"]] - t_start) / 60, N_WORKERS))

## ---- summary: excess errors, paired differences, integrity check -----------------------
paired_pp <- function(a, b) {
  d <- 100 * (a - b); d <- d[is.finite(d)]
  if (!length(d)) return(c(diff_pp = NA_real_, mcse_pp = NA_real_))
  c(diff_pp = mean(d), mcse_pp = sd(d) / sqrt(length(d)))
}
rows <- lapply(seq_len(nrow(tasks)), function(i) {
  t <- tasks[i, ]
  if (!file.exists(t$out)) return(NULL)
  r <- readRDS(t$out); z <- r$spec; new <- r$res
  st <- readRDS(t$file)$res[r$stored_rows, ]
  ex <- function(v) 100 * mean(v - st$Bayes, na.rm = TRUE)
  cmp <- list(
    PE_vs_MCD_pool  = paired_pp(st$PE1_pool_fallback, new$MCD_pool),
    PE_vs_MCD_sep   = paired_pp(st$PE1_sep_fallback,  new$MCD_sep),
    PE_vs_T_pool    = paired_pp(st$PE1_pool_fallback, new$T_pool),
    PE_vs_T_sep     = paired_pp(st$PE1_sep_fallback,  new$T_sep),
    MCD_vs_norm_pool = paired_pp(new$MCD_pool, st$LDA),
    MCD_vs_norm_sep  = paired_pp(new$MCD_sep,  st$QDA),
    T_vs_norm_pool   = paired_pp(new$T_pool,   st$LDA),
    T_vs_norm_sep    = paired_pp(new$T_sep,    st$QDA)
  )
  cmpv <- unlist(lapply(names(cmp), function(nm) setNames(cmp[[nm]], paste0(nm, "_", names(cmp[[nm]])))))
  data.frame(
    source = "B", cell_id = z$cell_id, block = z$block, n = z$n, p = z$p,
    beta1 = z$beta1, beta2 = z$beta2, dist = NA_character_,
    cov_type = z$cov_type, separation = z$separation, target = z$target, B = nrow(new),
    excess_pp_PE1_pool_fallback = ex(st$PE1_pool_fallback),
    excess_pp_PE1_sep_fallback  = ex(st$PE1_sep_fallback),
    excess_pp_LDA = ex(st$LDA), excess_pp_QDA = ex(st$QDA),
    excess_pp_MCD_pool = ex(new$MCD_pool), excess_pp_MCD_sep = ex(new$MCD_sep),
    excess_pp_T_pool = ex(new$T_pool), excess_pp_T_sep = ex(new$T_sep),
    fail_MCD_pool = mean(is.na(new$MCD_pool)), fail_MCD_sep = mean(is.na(new$MCD_sep)),
    fail_T_pool = mean(is.na(new$T_pool)), fail_T_sep = mean(is.na(new$T_sep)),
    median_nu_pool = median(new$nu_pool, na.rm = TRUE),
    median_nu1_sep = median(new$nu1_sep, na.rm = TRUE),
    median_nu2_sep = median(new$nu2_sep, na.rm = TRUE),
    repro_bayes_maxdiff = max(abs(new$Bayes_check - st$Bayes)),
    repro_lda_maxdiff = max(abs(new$LDA_check - st$LDA), na.rm = TRUE),
    as.data.frame(as.list(cmpv))
  )
})
summ <- do.call(rbind, Filter(Negate(is.null), rows))
write.csv(summ, file.path(OUT_DIR, "robust_summary.csv"), row.names = FALSE)
cat("Summary rows:", nrow(summ), "of", nrow(tasks), "\n")
cat("Integrity check (must be 0): max Bayes diff =", format(max(summ$repro_bayes_maxdiff)),
    " max LDA diff =", format(max(summ$repro_lda_maxdiff)), "\n")
if (max(summ$repro_bayes_maxdiff) > 1e-12 || max(summ$repro_lda_maxdiff) > 1e-12) {
  warning("The regenerated data differ from the stored data: paired comparisons are NOT valid.")
}
cat("Done.\n")
