# ============================================================================
# run_grid3_restartable.R : restartable simulation of Blocks A and B (288 cells)
# Completed calibration targets and simulation cells are cached individually.
# A valid completed cell is skipped on restart; invalid/incomplete files rerun.
# Block C (distributional misspecification) is not included.
# ============================================================================

suppressPackageStartupMessages(library(parallel))

## ---- settings ---------------------------------------------------------------
SCRIPT    <- "pe_discriminant_analysis_profile_optimized_5.R"
OUT_DIR   <- "grid_final"
B_REPS    <- 1000
N_WORKERS <- max(1L, detectCores() - 1L)
M_CAL     <- 1e6
VAR2      <- 1.25
RHO2      <- 0.5
DRY_RUN   <- FALSE

if (DRY_RUN) {
  OUT_DIR <- "grid_dryrun"
  B_REPS <- 5
  M_CAL <- 1e5
}

SCRIPT <- normalizePath(SCRIPT, mustWork = TRUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "cells"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "calibrations"), recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- normalizePath(OUT_DIR)
SCRIPT_MD5 <- unname(tools::md5sum(SCRIPT))

CAL_CONFIG <- list(
  script_md5 = SCRIPT_MD5, M_CAL = M_CAL, VAR2 = VAR2, RHO2 = RHO2
)
CELL_CONFIG <- list(
  script_md5 = SCRIPT_MD5, B_REPS = B_REPS, M_CAL = M_CAL,
  VAR2 = VAR2, RHO2 = RHO2
)

## ---- atomic file writing ----------------------------------------------------
atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(
    pattern = paste0(".", basename(path), ".tmp_"),
    tmpdir = dirname(path)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not rename temporary file to: ", path)
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  tmp <- tempfile(
    pattern = paste0(".", basename(path), ".tmp_"),
    tmpdir = dirname(path)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  write.csv(x, tmp, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not rename temporary CSV to: ", path)
  invisible(path)
}

## ---- design -----------------------------------------------------------------
targets <- c(high = 0.05, moderate = 0.15, low = 0.30)
cols <- c("block", "n", "p", "beta1", "beta2", "cov_type", "target")

gA <- expand.grid(
  n = c(50, 100, 250, 500), p = c(2, 5, 10),
  beta1 = c(0.5, 1, 2), target = unname(targets),
  cov_type = c("equal", "unequal"), stringsAsFactors = FALSE
)
gA$beta2 <- gA$beta1
gA$block <- "A"

shape_pairs <- data.frame(beta1 = c(0.5, 1, 0.5), beta2 = c(1, 2, 2))
gB <- merge(
  expand.grid(
    n = c(50, 250), p = c(2, 5), target = unname(targets),
    cov_type = c("equal", "unequal"), stringsAsFactors = FALSE
  ),
  shape_pairs
)
gB$block <- "B"

grid <- rbind(gA[, cols], gB[, cols])
grid$separation <- names(targets)[match(grid$target, unname(targets))]
if (DRY_RUN) grid <- grid[c(1, 80, 160, 216, 217, 250, 288), ]
grid$key <- with(grid, paste(p, beta1, beta2, cov_type, target, sep = "|"))

calkeys <- grid[
  !duplicated(grid$key),
  c("key", "p", "beta1", "beta2", "cov_type", "target")
]
calkeys$cal_id <- seq_len(nrow(calkeys))

## ---- calibration cache and worker function ---------------------------------
calibration_file <- function(cal_id) {
  file.path(OUT_DIR, "calibrations", sprintf("cal_%04d.rds", cal_id))
}

calibration_cache_valid <- function(path, z) {
  if (!file.exists(path)) return(FALSE)
  tryCatch({
    x <- readRDS(path)
    is.list(x) && identical(x$config, CAL_CONFIG) &&
      is.data.frame(x$result) && nrow(x$result) == 1L &&
      identical(as.character(x$result$key[1]), as.character(z$key))
  }, error = function(e) FALSE)
}

calibrate_one <- function(z) {
  f <- calibration_file(z$cal_id)
  if (calibration_cache_valid(f, z)) {
    out <- readRDS(f)$result
    out$status <- "cached"
    return(out)
  }

  out <- tryCatch({
    s <- envW$make_setup(
      z$p, z$beta1, z$beta2, z$cov_type, var2 = VAR2, rho2 = RHO2
    )
    r <- envW$calibrate_theta(
      s, target = z$target, M = M_CAL, M_validate = M_CAL,
      seed = 20261030 + z$cal_id
    )
    data.frame(
      key = z$key, theta = r$theta, validated_bayes = r$bayes,
      validation_mcse = r$validation_mcse, note = r$note,
      status = "computed", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      key = z$key, theta = NA_real_, validated_bayes = NA_real_,
      validation_mcse = NA_real_, note = conditionMessage(e),
      status = "error", stringsAsFactors = FALSE
    )
  })

  atomic_save_rds(list(config = CAL_CONFIG, result = out), f)
  out
}

## ---- cell cache and worker function ----------------------------------------
cell_file <- function(cell_id) {
  file.path(OUT_DIR, "cells", sprintf("cell_%04d_B%d.rds", cell_id, B_REPS))
}

same_cell_spec <- function(saved, current) {
  fields <- c(
    "cell_id", "block", "n", "p", "beta1", "beta2", "cov_type",
    "target", "theta", "validated_bayes", "seed", "key"
  )
  if (!all(fields %in% names(saved)) || !all(fields %in% names(current))) {
    return(FALSE)
  }
  all(vapply(fields, function(nm) {
    a <- saved[[nm]]
    b <- current[[nm]]
    if (is.numeric(a) && is.numeric(b)) {
      isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 0))
    } else {
      identical(as.character(a), as.character(b))
    }
  }, logical(1)))
}

cell_cache_valid <- function(path, z) {
  if (!file.exists(path)) return(FALSE)
  tryCatch({
    x <- readRDS(path)
    required <- c(
      "Bayes",
      "PE1_pool_fallback", "PE1_pool_accept",
      "PE1_sep_fallback", "PE1_sep_accept",
      "PE2_pool", "PE2_sep", "PE4_pool", "PE4_sep", "LDA", "QDA",
      "converged_true_pool", "converged_true_sep",
      "converged_pool", "boundary_pool", "singular_pool", "stalled_pool",
      "fallback_pool_fallback",
      "converged_sep", "boundary_sep", "singular_sep", "stalled_sep",
      "fallback_sep_fallback",
      "bhat_pool", "bhat1_sep", "bhat2_sep"
    )
    is.list(x) && identical(x$config, CELL_CONFIG) &&
      same_cell_spec(x$spec, z) && is.data.frame(x$res) &&
      nrow(x$res) == B_REPS && all(required %in% names(x$res))
  }, error = function(e) FALSE)
}

run_one <- function(z) {
  f <- cell_file(z$cell_id)
  if (cell_cache_valid(f, z)) {
    return(data.frame(cell_id = z$cell_id, status = "cached", seconds = NA_real_))
  }
  if (!is.finite(z$theta)) {
    return(data.frame(cell_id = z$cell_id, status = "no_theta", seconds = NA_real_))
  }

  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    s <- envW$make_setup(
      z$p, z$beta1, z$beta2, z$cov_type, var2 = VAR2, rho2 = RHO2
    )
    res <- envW$run_cell(
      s, n = z$n, theta = z$theta, B = B_REPS, seed = z$seed
    )
    if (!is.data.frame(res) || nrow(res) != B_REPS) {
      stop("run_cell did not return exactly B_REPS rows.")
    }
    atomic_save_rds(list(spec = z, res = res, config = CELL_CONFIG), f)
    data.frame(
      cell_id = z$cell_id, status = "done",
      seconds = proc.time()[["elapsed"]] - t0
    )
  }, error = function(e) {
    data.frame(
      cell_id = z$cell_id, status = paste("error:", conditionMessage(e)),
      seconds = NA_real_
    )
  })
}

## ---- start workers ----------------------------------------------------------
cl <- makeCluster(N_WORKERS)
clusterExport(
  cl,
  c(
    "SCRIPT", "OUT_DIR", "B_REPS", "M_CAL", "VAR2", "RHO2",
    "CAL_CONFIG", "CELL_CONFIG", "atomic_save_rds",
    "calibration_file", "calibration_cache_valid", "calibrate_one",
    "cell_file", "same_cell_spec", "cell_cache_valid", "run_one"
  ),
  envir = environment()
)
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(MASS))
  envW <- new.env()
  sys.source(SCRIPT, envir = envW)
}))

## ---- 1. calibration lookup: cached separately per target --------------------
cal_start <- proc.time()[["elapsed"]]
cal_list <- parLapplyLB(
  cl,
  lapply(seq_len(nrow(calkeys)), function(i) calkeys[i, ]),
  calibrate_one
)
cal <- do.call(rbind, cal_list)
cal <- cal[match(calkeys$key, cal$key), ]

atomic_save_rds(
  list(config = CAL_CONFIG, cal = cal),
  file.path(OUT_DIR, "calibration.rds")
)
atomic_write_csv(cal, file.path(OUT_DIR, "calibration.csv"))

cat(
  "Calibration:", sum(is.finite(cal$theta)), "of", nrow(cal),
  "targets reachable\n"
)
if (any(!is.finite(cal$theta))) {
  print(cal[!is.finite(cal$theta), c("key", "note")])
}
cat(sprintf(
  "Calibration stage: %.1f minutes\n",
  (proc.time()[["elapsed"]] - cal_start) / 60
))

## ---- 2. simulation cells ----------------------------------------------------
grid <- merge(grid, cal, by = "key", all.x = TRUE, sort = FALSE)
grid$.original_order <- seq_len(nrow(grid))
grid <- grid[order(
  -(grid$n * grid$p), grid$block, grid$beta1, grid$beta2,
  grid$cov_type, grid$target, grid$.original_order
), ]
grid$.original_order <- NULL
grid$cell_id <- seq_len(nrow(grid))
grid$seed <- 20271030 + grid$cell_id

atomic_write_csv(grid, file.path(OUT_DIR, "design_with_seeds.csv"))

cat("Simulation cells:", nrow(grid), "\n")
cat("Replications per cell:", B_REPS, "\n")
cat("Workers:", N_WORKERS, "\n")
cat("Output directory:", OUT_DIR, "\n")

t_start <- proc.time()[["elapsed"]]
run_log <- do.call(rbind, parLapplyLB(
  cl, lapply(seq_len(nrow(grid)), function(i) grid[i, ]), run_one
))
atomic_write_csv(run_log, file.path(OUT_DIR, "run_log.csv"))

status_counts <- table(run_log$status)
cat(
  "Cells this session:",
  paste(names(status_counts), as.integer(status_counts), collapse = "; "),
  "\n"
)
cat(sprintf(
  "Wall-clock time this session: %.1f minutes on %d workers\n",
  (proc.time()[["elapsed"]] - t_start) / 60, N_WORKERS
))

## ---- 3. rebuild summary from all valid completed cells ----------------------
rules <- c(
  "PE1_pool_fallback", "PE1_pool_accept",
  "PE1_sep_fallback", "PE1_sep_accept",
  "PE2_pool", "PE2_sep", "PE4_pool", "PE4_sep", "LDA", "QDA"
)
comps <- list(
  RQ1_pool = c("PE1_pool_fallback", "PE2_pool"),
  RQ1_sep = c("PE1_sep_fallback", "PE2_sep"),
  RQ2_pool = c("PE1_pool_fallback", "LDA"),
  RQ2_sep = c("PE1_sep_fallback", "QDA"),
  accept_pool = c("PE1_pool_accept", "PE1_pool_fallback"),
  accept_sep = c("PE1_sep_accept", "PE1_sep_fallback"),
  cov_PE1 = c("PE1_pool_fallback", "PE1_sep_fallback"),
  cov_PE2 = c("PE2_pool", "PE2_sep"),
  cov_normal = c("LDA", "QDA"),
  check_pool = c("PE2_pool", "LDA"),
  check_sep = c("PE2_sep", "QDA"),
  # true-shape contrasts (excess_pp_PE4_* gives R_true - R_B)
  shape_est_pool = c("PE1_pool_fallback", "PE4_pool"),   # R_est - R_true
  shape_est_sep = c("PE1_sep_fallback", "PE4_sep"),
  fixed_true_pool = c("PE2_pool", "PE4_pool"),           # R_fixed - R_true
  fixed_true_sep = c("PE2_sep", "PE4_sep")
)

paired_pp <- function(a, b) {
  d <- 100 * (a - b)
  d <- d[is.finite(d)]
  if (!length(d)) return(c(diff_pp = NA_real_, mcse_pp = NA_real_))
  c(diff_pp = mean(d), mcse_pp = sd(d) / sqrt(length(d)))
}

summarise_cell_file <- function(f) {
  x <- readRDS(f)
  z <- x$spec
  r <- x$res

  ex <- 100 * colMeans(r[rules] - r$Bayes, na.rm = TRUE)
  names(ex) <- paste0("excess_pp_", rules)

  pc <- unlist(lapply(names(comps), function(nm) {
    v <- paired_pp(r[[comps[[nm]][1]]], r[[comps[[nm]][2]]])
    setNames(v, paste0(nm, "_", names(v)))
  }))

  fl <- c(
    nonconv_pool = 1 - mean(r$converged_pool),
    boundary_pool = mean(r$boundary_pool),
    singular_pool = mean(r$singular_pool),
    stalled_pool = mean(r$stalled_pool),
    fallback_pool = mean(r$fallback_pool_fallback),
    nonconv_sep = 1 - mean(r$converged_sep),
    boundary_sep = mean(r$boundary_sep),
    singular_sep = mean(r$singular_sep),
    stalled_sep = mean(r$stalled_sep),
    fallback_sep = mean(r$fallback_sep_fallback),
    nonconv_true_pool = if (all(is.na(r$converged_true_pool))) NA_real_ else
      1 - mean(r$converged_true_pool, na.rm = TRUE),
    nonconv_true_sep = 1 - mean(r$converged_true_sep, na.rm = TRUE)
  )

  same <- z$beta1 == z$beta2
  bs <- c(
    mean_bhat_pool = mean(r$bhat_pool),
    bias_pool = if (same) mean(r$bhat_pool) - z$beta1 else NA_real_,
    rmse_pool = if (same) sqrt(mean((r$bhat_pool - z$beta1)^2)) else NA_real_,
    bias1_sep = mean(r$bhat1_sep) - z$beta1,
    bias2_sep = mean(r$bhat2_sep) - z$beta2,
    rmse1_sep = sqrt(mean((r$bhat1_sep - z$beta1)^2)),
    rmse2_sep = sqrt(mean((r$bhat2_sep - z$beta2)^2)),
    median_bhat1_sep = median(r$bhat1_sep)
  )

  cbind(
    z[, c(
      "cell_id", "block", "n", "p", "beta1", "beta2", "cov_type",
      "separation", "target", "theta", "validated_bayes", "seed"
    )],
    B = nrow(r),
    as.data.frame(as.list(c(ex, pc, fl, bs)))
  )
}

valid_rows <- lapply(seq_len(nrow(grid)), function(i) {
  z <- grid[i, ]
  f <- cell_file(z$cell_id)
  if (cell_cache_valid(f, z)) {
    data.frame(cell_id = z$cell_id, file = f, stringsAsFactors = FALSE)
  } else {
    NULL
  }
})
valid_rows <- Filter(Negate(is.null), valid_rows)

if (length(valid_rows)) {
  valid_files <- do.call(rbind, valid_rows)
  summary_df <- do.call(
    rbind, lapply(valid_files$file, summarise_cell_file)
  )
  summary_df <- summary_df[order(summary_df$cell_id), ]
  atomic_write_csv(summary_df, file.path(OUT_DIR, "cell_summary.csv"))
  cat(
    "Summary contains", nrow(summary_df), "of", nrow(grid),
    "valid completed cells.\n"
  )
} else {
  cat("No valid completed simulation cells are available for the summary yet.\n")
}

n_incomplete <- nrow(grid) - length(valid_rows)
if (n_incomplete > 0L) {
  cat(
    n_incomplete,
    "cells remain incomplete. Rerun this same script with unchanged settings;",
    "valid completed cells will be skipped.\n"
  )
} else {
  cat("All simulation cells are complete.\n")
}

stopCluster(cl)
cat("Done.\n")
