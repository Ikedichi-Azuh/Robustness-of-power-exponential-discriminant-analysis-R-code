# ============================================================================
# robust_comparators.R : robust comparator rules for the PE discriminant study
# Sourced INTO the same environment as pe_discriminant_analysis_profile_optimized_6.R
# (it uses that file's gen_data(), logdet(), err_rate(), bayes_predict(), ...).
#
# Rules (equal priors, plug-in classifiers):
#   MCD_pool / MCD_sep : Gaussian discriminant rule with class-wise minimum
#       covariance determinant (MCD) location and scatter (robustbase::covMcd,
#       reweighted, with its consistency and small-sample corrections).
#       pooled scatter = (n1 S1 + n2 S2) / (n1 + n2).
#   T_pool / T_sep : multivariate-t discriminant rule. Class-wise location and
#       scatter come from MASS::cov.trob for each degrees of freedom nu on the
#       grid NU_GRID; nu is chosen by maximum likelihood on that grid (the
#       analogue of the PE profile likelihood for the shape).
#       T_sep: nu_k chosen separately for each class.
#       T_pool: one common nu and one common scatter (weighted mean of the
#       class scatters), class-specific locations; nu maximises the summed t
#       log-likelihood with that common scatter.
# ============================================================================

NU_GRID <- c(2, 3, 5, 10, 20, 50)

if (!requireNamespace("robustbase", quietly = TRUE)) {
  stop("Package 'robustbase' is required: install.packages('robustbase')")
}

gauss_score <- function(X, m, S) {
  R <- sweep(as.matrix(X), 2, m, "-")
  -0.5 * logdet(S) - 0.5 * rowSums((R %*% solve(S)) * R)
}

t_const <- function(p, nu) {
  lgamma((nu + p) / 2) - lgamma(nu / 2) - (p / 2) * log(nu * pi)
}

t_score <- function(X, m, S, nu) {
  X <- as.matrix(X); p <- ncol(X)
  R <- sweep(X, 2, m, "-")
  Q <- rowSums((R %*% solve(S)) * R)
  t_const(p, nu) - 0.5 * logdet(S) - 0.5 * (nu + p) * log1p(Q / nu)
}

t_loglik <- function(X, m, S, nu) sum(t_score(X, m, S, nu))

# ---- MCD rules -----------------------------------------------------------------
fit_mcd <- function(X, y) {
  g1 <- X[y == 1, , drop = FALSE]; g2 <- X[y == 2, , drop = FALSE]
  f1 <- suppressWarnings(robustbase::covMcd(g1))
  f2 <- suppressWarnings(robustbase::covMcd(g2))
  n1 <- nrow(g1); n2 <- nrow(g2)
  Sp <- (n1 * f1$cov + n2 * f2$cov) / (n1 + n2)
  list(m = list(f1$center, f2$center),
       S_sep = list(f1$cov, f2$cov), S_pool = list(Sp, Sp))
}

predict_gauss <- function(m, S, X) {
  ifelse(gauss_score(X, m[[1]], S[[1]]) > gauss_score(X, m[[2]], S[[2]]), 1L, 2L)
}

# ---- t rules ---------------------------------------------------------------------
fit_t_grid <- function(X, nus = NU_GRID) {
  lapply(nus, function(nu) {
    f <- suppressWarnings(MASS::cov.trob(X, nu = nu, maxit = 200, tol = 1e-6))
    list(nu = nu, m = f$center, S = f$cov, ll = t_loglik(X, f$center, f$cov, nu))
  })
}

fit_t <- function(X, y, nus = NU_GRID) {
  g1 <- X[y == 1, , drop = FALSE]; g2 <- X[y == 2, , drop = FALSE]
  n1 <- nrow(g1); n2 <- nrow(g2)
  fits1 <- fit_t_grid(g1, nus); fits2 <- fit_t_grid(g2, nus)

  # class-specific: best nu in each class
  b1 <- fits1[[which.max(vapply(fits1, function(f) f$ll, numeric(1)))]]
  b2 <- fits2[[which.max(vapply(fits2, function(f) f$ll, numeric(1)))]]

  # pooled: common nu and common scatter, class-specific locations
  ll_pool <- vapply(seq_along(nus), function(j) {
    Sp <- (n1 * fits1[[j]]$S + n2 * fits2[[j]]$S) / (n1 + n2)
    t_loglik(g1, fits1[[j]]$m, Sp, nus[j]) + t_loglik(g2, fits2[[j]]$m, Sp, nus[j])
  }, numeric(1))
  j <- which.max(ll_pool)
  Sp <- (n1 * fits1[[j]]$S + n2 * fits2[[j]]$S) / (n1 + n2)

  list(sep = list(m = list(b1$m, b2$m), S = list(b1$S, b2$S), nu = c(b1$nu, b2$nu)),
       pool = list(m = list(fits1[[j]]$m, fits2[[j]]$m), S = list(Sp, Sp),
                   nu = c(nus[j], nus[j])))
}

predict_t <- function(fit, X) {
  ifelse(t_score(X, fit$m[[1]], fit$S[[1]], fit$nu[1]) >
         t_score(X, fit$m[[2]], fit$S[[2]], fit$nu[2]), 1L, 2L)
}

# ---- one replicate ----------------------------------------------------------------
# Draws the training and test sets in exactly the same order as one_rep() in the
# main script, so that with the same seed they are the same data. The random
# state is saved after the draws and restored on exit, because the MCD fit uses
# random subsets and must not disturb the stream used by later replicates.
robust_rep <- function(s, n, theta, N = 10000) {
  tr <- gen_data(s, n, n, theta)
  te <- gen_data(s, N / 2, N / 2, theta)
  # The ORIGINAL one_rep() consumes random numbers after the data draws: MASS
  # predict.lda / predict.qda break exact ties with max.col(ties.method = "random").
  # To keep the random stream identical to the stored run, the LDA and QDA fits and
  # predictions are repeated here in the same order (lda fit, qda fit, lda
  # prediction, qda prediction). Only the LDA result is used (LDA_check).
  yf <- factor(tr$y, levels = c(1, 2))
  lda_fit <- tryCatch(MASS::lda(tr$X, yf, prior = c(0.5, 0.5)), error = function(e) NULL)
  qda_fit <- tryCatch(MASS::qda(tr$X, yf, prior = c(0.5, 0.5)), error = function(e) NULL)
  pred_lda <- if (is.null(lda_fit)) NULL else as.integer(predict(lda_fit, te$X)$class)
  if (!is.null(qda_fit)) invisible(predict(qda_fit, te$X)$class)

  # Save the stream state AFTER those draws. The MCD fit uses random subsets; the
  # state is restored on exit so later replicates are not disturbed.
  seed_state <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", seed_state, envir = .GlobalEnv), add = TRUE)

  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)

  mcd <- tryCatch(fit_mcd(tr$X, tr$y), error = function(e) NULL)
  tt  <- tryCatch(fit_t(tr$X, tr$y),   error = function(e) NULL)

  c(
    Bayes_check = err_rate(bayes_predict(s, te$X, theta), te$y),
    LDA_check   = if (is.null(pred_lda)) NA_real_ else err_rate(pred_lda, te$y),
    MCD_pool = if (is.null(mcd)) NA_real_ else
                 safe(err_rate(predict_gauss(mcd$m, mcd$S_pool, te$X), te$y)),
    MCD_sep  = if (is.null(mcd)) NA_real_ else
                 safe(err_rate(predict_gauss(mcd$m, mcd$S_sep, te$X), te$y)),
    T_pool   = if (is.null(tt)) NA_real_ else
                 safe(err_rate(predict_t(tt$pool, te$X), te$y)),
    T_sep    = if (is.null(tt)) NA_real_ else
                 safe(err_rate(predict_t(tt$sep, te$X), te$y)),
    nu_pool  = if (is.null(tt)) NA_real_ else tt$pool$nu[1],
    nu1_sep  = if (is.null(tt)) NA_real_ else tt$sep$nu[1],
    nu2_sep  = if (is.null(tt)) NA_real_ else tt$sep$nu[2]
  )
}

# ---- checks you should run before the big job ----------------------------------------
# (a) cov.trob returns the t SCATTER matrix, not the covariance;
# (b) covMcd is consistent for the covariance of normal data;
# (c) the saved-seed trick reproduces the stored data exactly.
self_test_robust <- function(stored_cell_file = NULL, B_check = 5) {
  set.seed(1)
  p <- 3
  Z <- matrix(rnorm(2e5 * p), ncol = p) / sqrt(rchisq(2e5, 5) / 5)   # t5, scatter = I
  f <- suppressWarnings(MASS::cov.trob(Z, nu = 5, maxit = 200, tol = 1e-8))
  cat("(a) cov.trob scatter for t5 data with scatter I (expect about 1 on the diagonal;",
      "5/3 = 1.67 would mean it returns the covariance):", round(diag(f$cov), 3), "\n")
  V <- cs_matrix(p, 0.5)
  X <- matrix(rnorm(2e4 * p), ncol = p) %*% chol(V)
  cat("(b) covMcd on normal data, max|cov - V| =",
      round(max(abs(suppressWarnings(robustbase::covMcd(X))$cov - V)), 4),
      "(expect < 0.05)\n")
  if (!is.null(stored_cell_file)) {
    x <- readRDS(stored_cell_file); z <- x$spec
    s <- if ("dist" %in% names(z)) {
      make_setup_C(z$p, z$dist, z$cov_type, var2 = 1.25, rho2 = 0.5)
    } else {
      make_setup(z$p, z$beta1, z$beta2, z$cov_type, var2 = 1.25, rho2 = 0.5)
    }
    set.seed(z$seed)
    r <- t(replicate(B_check, robust_rep(s, z$n, z$theta)))
    cat("(c) reproduction of stored data, max|Bayes diff| =",
        format(max(abs(r[, "Bayes_check"] - x$res$Bayes[1:B_check]))),
        " max|LDA diff| =",
        format(max(abs(r[, "LDA_check"] - x$res$LDA[1:B_check]))), "(expect 0)\n")
    print(round(r[, c("MCD_pool", "MCD_sep", "T_pool", "T_sep", "nu_pool")], 4))
  }
  invisible(TRUE)
}
