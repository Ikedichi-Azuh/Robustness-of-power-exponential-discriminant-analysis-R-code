# Power-exponential discriminant analysis simulation
# Profile-optimised pilot script: faster approximate profile evaluations followed
# by a full-tolerance final fit, plus the corrected diagnostics and self-tests.

suppressPackageStartupMessages(library(MASS))

## ---- 1. PE basics ---------------------------------------------------------
cfac <- function(p, b) {
  2^(1 / b) * gamma((p + 2) / (2 * b)) /
    (p * gamma(p / (2 * b)))
}

logC <- function(p, b) {
  log(p) + lgamma(p / 2) - (p / 2) * log(pi) -
    (1 + p / (2 * b)) * log(2) - lgamma(1 + p / (2 * b))
}

logdet <- function(S) {
  z <- determinant(S, logarithm = TRUE)
  if (z$sign <= 0) stop("Scale matrix is not positive definite.")
  as.numeric(z$modulus)
}

cs_matrix <- function(p, rho, var = 1) {
  var * ((1 - rho) * diag(p) + rho * matrix(1, p, p))
}

# X = mu + R A U, with A A' = Sigma and R^2 = (2W)^(1/b).
rpe <- function(n, mu, Sig, b) {
  p <- length(mu)
  A <- t(chol(Sig))
  U <- matrix(rnorm(n * p), nrow = n, ncol = p)
  U <- U / sqrt(rowSums(U^2))
  R <- (2 * rgamma(n, shape = p / (2 * b), rate = 1))^(1 / (2 * b))
  sweep((R * U) %*% t(A), 2, mu, "+")
}

# Log density; retain every normalising term because classes may differ in beta/Sigma.
lpe <- function(X, mu, Sig, b) {
  X <- as.matrix(X)
  d <- mahalanobis(X, center = mu, cov = Sig)
  logC(ncol(X), b) - 0.5 * logdet(Sig) - 0.5 * d^b
}

## ---- 2. Population set-up --------------------------------------------------
# Equal covariance: V1 = V2 = CS(rho = 0.5).
# Unequal covariance: V1 = CS(rho = 0.5),
# V2 = var2 * CS(rho = rho2). var2 and rho2 are ignored for cov_type="equal".
make_setup <- function(p, b1, b2 = b1,
                       cov_type = c("equal", "unequal"),
                       var2 = 2, rho2 = 0.2) {
  cov_type <- match.arg(cov_type)

  if (length(p) != 1L || !is.numeric(p) || !is.finite(p) ||
      p < 1 || p != as.integer(p)) {
    stop("p must be a single positive integer.")
  }
  if (length(var2) != 1L || !is.numeric(var2) ||
      !is.finite(var2) || var2 <= 0) {
    stop("var2 must be a single positive finite number.")
  }
  if (length(rho2) != 1L || !is.numeric(rho2) || !is.finite(rho2)) {
    stop("rho2 must be a single finite number.")
  }
  rho_lower <- if (p == 1) -1 else -1 / (p - 1)
  if (rho2 <= rho_lower || rho2 >= 1) {
    stop(sprintf("rho2 must lie in (%.6g, 1) for p = %d.", rho_lower, p))
  }

  V1 <- cs_matrix(p, 0.5)
  V2 <- if (cov_type == "equal") V1 else cs_matrix(p, rho2, var = var2)

  list(
    p = p, b1 = b1, b2 = b2, V1 = V1, V2 = V2,
    S1 = V1 / cfac(p, b1),
    S2 = V2 / cfac(p, b2),
    d = rep(1, p) / sqrt(p),
    cov_type = cov_type
  )
}

gen_data <- function(s, n1, n2, theta) {
  X1 <- rpe(n1, rep(0, s$p), s$S1, s$b1)
  X2 <- rpe(n2, theta * s$d, s$S2, s$b2)
  list(X = rbind(X1, X2), y = c(rep(1L, n1), rep(2L, n2)))
}

## ---- 3. True-parameter Bayes rule and calibration -------------------------
bayes_predict <- function(s, X, theta) {
  l1 <- lpe(X, rep(0, s$p), s$S1, s$b1)
  l2 <- lpe(X, theta * s$d, s$S2, s$b2)
  ifelse(l1 > l2, 1L, 2L)
}

bayes_error_mc <- function(s, theta, M = 1e6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Z1 <- rpe(M, rep(0, s$p), s$S1, s$b1)
  Z2 <- rpe(M, rep(0, s$p), s$S2, s$b2)
  X2 <- sweep(Z2, 2, theta * s$d, "+")
  e1 <- mean(bayes_predict(s, Z1, theta) != 1L)
  e2 <- mean(bayes_predict(s, X2, theta) != 2L)
  e <- 0.5 * (e1 + e2)
  # Approximate MCSE, treating the two class-specific Monte Carlo samples as independent.
  se <- sqrt(0.25 * (e1 * (1 - e1) + e2 * (1 - e2)) / M)
  c(error = e, mcse = se, class1_error = e1, class2_error = e2)
}

calibrate_theta <- function(s, target, M = 2e5, upper = 30, seed = 1,
                            M_validate = 1e6, max_upper = 1000) {
  if (!is.numeric(target) || length(target) != 1L ||
      !is.finite(target) || target <= 0 || target >= 0.5) {
    stop("target must be a single finite number strictly between 0 and 0.5.")
  }
  if (!is.numeric(upper) || length(upper) != 1L ||
      !is.finite(upper) || upper <= 0) {
    stop("upper must be a single positive finite number.")
  }
  if (!is.numeric(max_upper) || length(max_upper) != 1L ||
      !is.finite(max_upper) || max_upper < upper) {
    stop("max_upper must be finite and at least as large as upper.")
  }
  set.seed(seed)
  Z1 <- rpe(M, rep(0, s$p), s$S1, s$b1)
  Z2 <- rpe(M, rep(0, s$p), s$S2, s$b2)

  # Common random numbers reduce Monte Carlo noise during calibration.
  err <- function(theta) {
    X2 <- sweep(Z2, 2, theta * s$d, "+")
    e1 <- mean(bayes_predict(s, Z1, theta) != 1L)
    e2 <- mean(bayes_predict(s, X2, theta) != 2L)
    0.5 * (e1 + e2)
  }

  e_lo <- err(0)
  if (e_lo < target) {
    return(list(
      theta = NA_real_, bayes = NA_real_, calibration_bayes = e_lo,
      validation_mcse = NA_real_, lower = 0, upper = upper,
      note = sprintf(
        paste0("Target %.3f is below the estimated Bayes error at theta = 0 ",
               "(%.3f); the target is not bracketed at theta = 0 under the assumed decreasing-error calibration."),
        target, e_lo
      )
    ))
  }

  e_hi <- err(upper)
  while (e_hi > target && upper < max_upper) {
    upper <- min(2 * upper, max_upper)
    e_hi <- err(upper)
  }
  if (e_hi > target) {
    return(list(
      theta = NA_real_, bayes = NA_real_, calibration_bayes = e_hi,
      validation_mcse = NA_real_, lower = 0, upper = upper,
      note = sprintf(
        paste0("Target %.3f was not bracketed by theta in [0, %.1f]; ",
               "estimated errors were %.3f at 0 and %.3f at the upper bound. ",
               "Try a different covariance structure or inspect calibration noise."),
        target, upper, e_lo, e_hi
      )
    ))
  }

  root <- uniroot(function(t) err(t) - target, c(0, upper), tol = 1e-4)$root
  calibration_estimate <- err(root)

  # Independent validation sample: do not report the calibration sample's
  # near-target error as if it were independent confirmation.
  validation <- bayes_error_mc(
    s, root, M = M_validate,
    seed = if (is.null(seed)) NULL else seed + 1L
  )

  list(
    theta = root,
    calibration_lower = 0,
    calibration_upper = upper,
    bayes = unname(validation["error"]),
    calibration_bayes = calibration_estimate,
    validation_mcse = unname(validation["mcse"]),
    validation_class1_error = unname(validation["class1_error"]),
    validation_class2_error = unname(validation["class2_error"]),
    note = "ok"
  )
}

## ---- 4. PE maximum-likelihood fitting -------------------------------------
# Fit common Sigma and beta with group-specific locations.
# For K=1, this is a single-class fit; for K=2, it is a pooled fit.
#
# Convergence is determined from the final fixed-beta optimisation, not from
# failures encountered at unrelated beta values during profile evaluation.
fit_core <- function(X, g, b = 1, est_beta = FALSE,
                     bounds = c(0.2, 5), max_inner = 500,
                     tol = 1e-8, profile_tol = 1e-4,
                     profile_fit_tol = 1e-5) {
  if (!is.finite(profile_fit_tol) || profile_fit_tol <= 0) {
    stop("profile_fit_tol must be a positive finite number.")
  }
  X <- as.matrix(X)
  g <- as.integer(g)
  n <- nrow(X)
  p <- ncol(X)
  K <- max(g)

  if (length(g) != n || anyNA(g) || min(g) < 1 ||
      !all(seq_len(K) %in% unique(g))) {
    stop("g must be a complete group index from 1 to K.")
  }
  if (n <= p) {
    return(list(fail = TRUE, mus = matrix(NA_real_, K, p),
                S = diag(p), b = b, ll = NA_real_,
                converged = FALSE, boundary = FALSE, singular = TRUE))
  }

  # Cache group row indices once; avoid rebuilding them in every iteration.
  group_idx <- lapply(seq_len(K), function(k) which(g == k))
  initial_mus <- do.call(rbind, lapply(group_idx, function(idx) {
    colMeans(X[idx, , drop = FALSE])
  }))
  residual0 <- X - initial_mus[g, , drop = FALSE]
  initial_S <- crossprod(residual0) / n

  if (any(!is.finite(initial_S)) || rcond(initial_S) < 1e-12) {
    return(list(fail = TRUE, mus = initial_mus, S = initial_S, b = b,
                ll = NA_real_, converged = FALSE, boundary = FALSE,
                singular = TRUE))
  }

  # Evaluate likelihood and matrix quantities together. A single Cholesky
  # decomposition supplies both log|S| and S^{-1}; callers can reuse these.
  evaluate <- function(mus, S, beta) {
    if (!is.finite(beta) || beta <= 0 || any(!is.finite(S))) return(NULL)
    U <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(U)) return(NULL)
    ld <- 2 * sum(log(diag(U)))
    invS <- tryCatch(chol2inv(U), error = function(e) NULL)
    if (is.null(invS) || !is.finite(ld) || any(!is.finite(invS))) return(NULL)
    R <- X - mus[g, , drop = FALSE]
    d <- rowSums((R %*% invS) * R)
    if (any(!is.finite(d))) return(NULL)
    d <- pmax(d, 1e-12)
    ll <- n * logC(p, beta) - 0.5 * n * ld - 0.5 * sum(d^beta)
    if (!is.finite(ll)) return(NULL)
    list(ll = ll, invS = invS, logdet = ld)
  }

  # Profile evaluations only need a sufficiently accurate likelihood to locate
  # beta. Use a looser convergence tolerance during those repeated fits, then
  # refit the selected beta using the full `tol` below.
  inner <- function(mus_start, S_start, beta, inner_tol = tol) {
    mus <- mus_start
    S <- S_start
    ev <- evaluate(mus, S, beta)
    if (is.null(ev)) {
      return(list(mus = mus, S = S, ll = -Inf, converged = FALSE,
                  stalled = FALSE, iterations = 0L))
    }
    cur <- ev$ll
    converged <- FALSE
    stalled <- FALSE
    iter_used <- 0L

    for (it in seq_len(max_inner)) {
      iter_used <- it
      # Reuse the current iteration's inverse rather than decomposing S again.
      R <- X - mus[g, , drop = FALSE]
      d <- pmax(rowSums((R %*% ev$invS) * R), 1e-12)
      w <- d^(beta - 1)
      if (any(!is.finite(w)) || any(!is.finite(d))) {
        stalled <- TRUE
        break
      }

      mus_new <- do.call(rbind, lapply(group_idx, function(idx) {
        colSums(X[idx, , drop = FALSE] * w[idx]) / sum(w[idx])
      }))
      R_new <- X - mus_new[g, , drop = FALSE]
      S_new <- (beta / n) * crossprod(R_new * w, R_new)

      if (any(!is.finite(S_new)) || rcond(S_new) < 1e-14) {
        stalled <- TRUE
        break
      }

      # First trial step 1/beta (capped at 1). Along the scale direction the
      # undamped update has slope 1 - beta, so beta = 2 oscillates forever;
      # a step of 1/beta makes the slope 0 and converges quickly for beta > 1.
      step <- min(1, 1 / beta)
      accepted <- FALSE
      ev_try <- NULL
      while (step >= 1 / 4096) {
        mus_try <- mus + step * (mus_new - mus)
        S_try <- (1 - step) * S + step * S_new
        ev_try <- evaluate(mus_try, S_try, beta)
        if (!is.null(ev_try) && ev_try$ll >= cur - 1e-10) {
          accepted <- TRUE
          break
        }
        step <- step / 2
      }

      if (!accepted) {
        mu_change <- max(abs(mus_new - mus)) / (1 + max(abs(mus)))
        S_change <- max(abs(S_new - S)) / (1 + max(abs(S)))
        if (is.finite(mu_change) && is.finite(S_change) &&
            max(mu_change, S_change) < sqrt(inner_tol)) {
          converged <- TRUE
        } else {
          stalled <- TRUE
        }
        break
      }

      mu_change <- max(abs(mus_try - mus)) / (1 + max(abs(mus)))
      S_change <- max(abs(S_try - S)) / (1 + max(abs(S)))
      ll_change <- abs(ev_try$ll - cur) / (1 + abs(cur))

      mus <- mus_try
      S <- S_try
      ev <- ev_try
      cur <- ev$ll

      if (max(mu_change, S_change, ll_change) < inner_tol) {
        converged <- TRUE
        break
      }
    }

    list(mus = mus, S = S, ll = cur, converged = converged,
         stalled = stalled, iterations = iter_used)
  }

  fixed_fit <- tryCatch(inner(initial_mus, initial_S, b),
                        error = function(e) NULL)
  if (is.null(fixed_fit)) {
    return(list(fail = TRUE, mus = initial_mus, S = initial_S, b = b,
                ll = NA_real_, converged = FALSE, boundary = FALSE,
                singular = FALSE))
  }

  if (est_beta) {
    # Deterministic profile objective: every candidate starts from the same
    # initial values, avoiding evaluation-order-dependent warm starts.
    profile_objective <- function(beta_try) {
      f <- tryCatch(inner(initial_mus, initial_S, beta_try,
                          inner_tol = profile_fit_tol),
                    error = function(e) NULL)
      if (is.null(f) || !is.finite(f$ll)) return(.Machine$double.xmax / 100)
      -f$ll
    }
    opt <- tryCatch(
      optimize(profile_objective, interval = bounds, tol = profile_tol),
      error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$minimum)) {
      return(list(fail = TRUE, mus = fixed_fit$mus, S = fixed_fit$S,
                  b = b, ll = fixed_fit$ll, converged = FALSE,
                  boundary = FALSE, singular = FALSE))
    }
    b <- opt$minimum
    fixed_fit <- tryCatch(inner(initial_mus, initial_S, b),
                          error = function(e) NULL)
    if (is.null(fixed_fit)) {
      return(list(fail = TRUE, mus = initial_mus, S = initial_S, b = b,
                  ll = NA_real_, converged = FALSE,
                  boundary = FALSE, singular = FALSE))
    }
  }

  boundary <- est_beta &&
    (b <= bounds[1] + 1e-3 || b >= bounds[2] - 1e-3)
  rcond_S <- tryCatch(rcond(fixed_fit$S), error = function(e) NA_real_)
  singular <- !is.finite(rcond_S) || rcond_S < 1e-10
  converged <- isTRUE(fixed_fit$converged)
  fail <- !converged || boundary || singular

  list(
    mus = fixed_fit$mus, S = fixed_fit$S, b = b, ll = fixed_fit$ll,
    fail = fail, converged = converged, boundary = boundary,
    singular = singular, stalled = fixed_fit$stalled,
    iterations = fixed_fit$iterations
  )
}

# One PE rule. pooled = common Sigma/beta; otherwise class-specific Sigma/beta.
# On estimated-beta failure, use the corresponding beta=1 fit as primary
# analysis fallback, while retaining a failure flag for reporting.
fit_pe_rule <- function(
    X, y, pooled, est_beta,
    boundary_policy = c("accept", "fallback"),
    b_fixed = c(1, 1)
) {
  # b_fixed: shape(s) used when est_beta = FALSE. Default (1, 1) reproduces the
  # fixed-beta=1 fit. For a true-shape fit pass the data-generating shape(s);
  # a pooled fit uses b_fixed[1] only and is meaningful only if both are equal.
  stopifnot(length(b_fixed) == 2L, all(is.finite(b_fixed)), all(b_fixed > 0))
  boundary_policy <- match.arg(boundary_policy)
  
  diagnostics <- function(f) {
    list(
      converged = isTRUE(f$converged),
      boundary = isTRUE(f$boundary),
      singular = isTRUE(f$singular),
      stalled = isTRUE(f$stalled),
      iterations = if (is.null(f$iterations)) {
        NA_integer_
      } else {
        f$iterations
      },
      beta = f$b
    )
  }
  
  one <- function(estimate_beta) {
    if (pooled) {
      f <- fit_core(X, y, b = b_fixed[1], est_beta = estimate_beta)
      
      list(
        mus = f$mus,
        S = list(f$S, f$S),
        b = c(f$b, f$b),
        fail = isTRUE(f$fail),
        converged = isTRUE(f$converged),
        diagnostics = diagnostics(f)
      )
      
    } else {
      f1 <- fit_core(
        X[y == 1, , drop = FALSE],
        rep(1L, sum(y == 1)),
        b = b_fixed[1],
        est_beta = estimate_beta
      )
      
      f2 <- fit_core(
        X[y == 2, , drop = FALSE],
        rep(1L, sum(y == 2)),
        b = b_fixed[2],
        est_beta = estimate_beta
      )
      
      d1 <- diagnostics(f1)
      d2 <- diagnostics(f2)
      
      list(
        mus = rbind(f1$mus[1, ], f2$mus[1, ]),
        S = list(f1$S, f2$S),
        b = c(f1$b, f2$b),
        fail = isTRUE(f1$fail) || isTRUE(f2$fail),
        converged = isTRUE(f1$converged) &&
          isTRUE(f2$converged),
        diagnostics = list(
          converged = d1$converged && d2$converged,
          boundary = d1$boundary || d2$boundary,
          singular = d1$singular || d2$singular,
          stalled = d1$stalled || d2$stalled,
          iterations = d1$iterations + d2$iterations,
          beta = c(d1$beta, d2$beta)
        )
      )
    }
  }
  
  m <- one(est_beta)
  m$fallback <- FALSE
  m$estimated_diagnostics <- NULL
  m$estimated_b <- if (est_beta) m$b else NULL
  m$fallback_failed <- FALSE
  
  if (est_beta) {
    # Save the estimated fit before any fallback.
    estimated_fit <- m
    m$estimated_diagnostics <- estimated_fit$diagnostics
    m$estimated_b <- estimated_fit$b
    
    d <- estimated_fit$diagnostics
    
    # A boundary estimate is handled according to the selected policy.
    genuinely_failed <- !estimated_fit$converged ||
      d$singular ||
      d$stalled
    
    use_fallback <- genuinely_failed ||
      (isTRUE(d$boundary) && boundary_policy == "fallback")
    
    if (use_fallback) {
      fallback_fit <- one(FALSE)
      
      m$fallback <- TRUE
      m$fallback_failed <- isTRUE(fallback_fit$fail) ||
        !isTRUE(fallback_fit$converged) ||
        isTRUE(fallback_fit$diagnostics$singular) ||
        isTRUE(fallback_fit$diagnostics$stalled)
      
      if (!m$fallback_failed) {
        m$mus <- fallback_fit$mus
        m$S <- fallback_fit$S
        m$b <- fallback_fit$b
      } else {
        # Do not allow an unusable fallback to be treated as a valid fit.
        m$mus <- NULL
        m$S <- NULL
        m$b <- NULL
      }
      
      # Keep the estimated-fit diagnostics separate from fallback diagnostics.
      m$fallback_diagnostics <- fallback_fit$diagnostics
      
      # A fallback does not erase the fact that the estimated fit had a problem.
      m$fail <- TRUE
    }
    
    m$boundary_policy <- boundary_policy
    
  } else {
    m$fallback_diagnostics <- NULL
    m$boundary_policy <- NA_character_
  }
  
  m
}
predict_pe <- function(m, X) {
  l1 <- lpe(X, m$mus[1, ], m$S[[1]], m$b[1])
  l2 <- lpe(X, m$mus[2, ], m$S[[2]], m$b[2])
  ifelse(l1 > l2, 1L, 2L)
}

## ---- 5. Gaussian MLE classifier for exact implementation checks ------------
# This deliberately uses maximum-likelihood covariance estimates (divisor n_k),
# matching PE with beta=1 rather than MASS::qda's n_k - 1 convention.
fit_gaussian_mle <- function(X, y, pooled = FALSE) {
  X <- as.matrix(X)
  mu1 <- colMeans(X[y == 1, , drop = FALSE])
  mu2 <- colMeans(X[y == 2, , drop = FALSE])

  if (pooled) {
    R <- rbind(
      sweep(X[y == 1, , drop = FALSE], 2, mu1, "-"),
      sweep(X[y == 2, , drop = FALSE], 2, mu2, "-")
    )
    S <- crossprod(R) / nrow(X)
    Ss <- list(S, S)
  } else {
    R1 <- sweep(X[y == 1, , drop = FALSE], 2, mu1, "-")
    R2 <- sweep(X[y == 2, , drop = FALSE], 2, mu2, "-")
    Ss <- list(crossprod(R1) / nrow(R1), crossprod(R2) / nrow(R2))
  }
  list(mus = rbind(mu1, mu2), S = Ss)
}

predict_gaussian_mle <- function(fit, X) {
  X <- as.matrix(X)
  score <- vapply(1:2, function(k) {
    R <- sweep(X, 2, fit$mus[k, ], "-")
    -0.5 * logdet(fit$S[[k]]) -
      0.5 * rowSums((R %*% solve(fit$S[[k]])) * R)
  }, numeric(nrow(X)))
  if (is.null(dim(score))) score <- matrix(score, ncol = 2)
  ifelse(score[, 1] > score[, 2], 1L, 2L)
}

## ---- 6. One replicate and one scenario cell --------------------------------
err_rate <- function(pred, y) {
  0.5 * (mean(pred[y == 1] != 1L) + mean(pred[y == 2] != 2L))
}


one_rep <- function(s, n, theta, N = 10000) {
  
  tr <- gen_data(s, n, n, theta)
  te <- gen_data(s, N / 2, N / 2, theta)
  
  yf <- factor(tr$y, levels = c(1, 2))
  
  # Estimated-shape fits under both boundary policies.
  fit_est <- function(pooled, policy) {
    fit_pe_rule(
      tr$X, tr$y,
      pooled = pooled,
      est_beta = TRUE,
      boundary_policy = policy
    )
  }
  
  pool_accept <- fit_est(TRUE, "accept")
  pool_fallback <- fit_est(TRUE, "fallback")
  
  sep_accept <- fit_est(FALSE, "accept")
  sep_fallback <- fit_est(FALSE, "fallback")
  
  # Fixed-shape beta = 1 fits.
  m2p <- fit_pe_rule(
    tr$X, tr$y, pooled = TRUE, est_beta = FALSE
  )
  m2s <- fit_pe_rule(
    tr$X, tr$y, pooled = FALSE, est_beta = FALSE
  )
  
  # True-shape diagnostic fits: shape fixed at the data-generating value(s),
  # all other parameters estimated. The pooled fit exists only for equal true
  # shapes; it is NA otherwise. These are benchmarks, not Bayes rules.
  m4s <- fit_pe_rule(
    tr$X, tr$y, pooled = FALSE, est_beta = FALSE,
    b_fixed = c(s$b1, s$b2)
  )
  m4p <- if (isTRUE(all.equal(s$b1, s$b2))) {
    fit_pe_rule(
      tr$X, tr$y, pooled = TRUE, est_beta = FALSE,
      b_fixed = c(s$b1, s$b2)
    )
  } else {
    NULL
  }
  
  # Return NA if a fit cannot safely be used for prediction.
  safe_error <- function(m) {
    if (isTRUE(m$fallback_failed) ||
        is.null(m$mus) ||
        is.null(m$S) ||
        is.null(m$b)) {
      return(NA_real_)
    }
    
    tryCatch(
      err_rate(predict_pe(m, te$X), te$y),
      error = function(e) NA_real_
    )
  }
  
  lda_fit <- tryCatch(
    lda(tr$X, yf, prior = c(0.5, 0.5)),
    error = function(e) NULL
  )
  
  qda_fit <- tryCatch(
    qda(tr$X, yf, prior = c(0.5, 0.5)),
    error = function(e) NULL
  )
  
  pred_lda <- if (is.null(lda_fit)) {
    rep(NA_integer_, nrow(te$X))
  } else {
    as.integer(predict(lda_fit, te$X)$class)
  }
  
  pred_qda <- if (is.null(qda_fit)) {
    rep(NA_integer_, nrow(te$X))
  } else {
    as.integer(predict(qda_fit, te$X)$class)
  }
  
  # Estimated-fit diagnostics, independent of fallback parameters.
  get_diag <- function(m, field) {
    d <- m$estimated_diagnostics
    if (is.null(d)) d <- m$diagnostics
    
    value <- d[[field]]
    
    if (is.null(value) || length(value) != 1L) {
      return(NA_real_)
    }
    
    as.numeric(value)
  }
  
  # Raw beta estimates: never replace these with fallback beta = 1.
  raw_beta <- function(m, j = 1L) {
    b <- m$estimated_b
    if (is.null(b) || length(b) < j) return(NA_real_)
    as.numeric(b[j])
  }
  
  r <- c(
    Bayes = err_rate(
      bayes_predict(s, te$X, theta), te$y
    ),
    
    PE1_pool_accept = safe_error(pool_accept),
    PE1_pool_fallback = safe_error(pool_fallback),
    
    PE1_sep_accept = safe_error(sep_accept),
    PE1_sep_fallback = safe_error(sep_fallback),
    
    PE2_pool = safe_error(m2p),
    PE2_sep = safe_error(m2s),
    
    PE4_pool = if (is.null(m4p)) NA_real_ else safe_error(m4p),
    PE4_sep = safe_error(m4s),
    
    LDA = if (is.null(lda_fit)) {
      NA_real_
    } else {
      err_rate(pred_lda, te$y)
    },
    
    QDA = if (is.null(qda_fit)) {
      NA_real_
    } else {
      err_rate(pred_qda, te$y)
    }
  )
  
  info <- c(
    # Raw estimated shape parameters
    bhat_pool = raw_beta(pool_accept),
    bhat1_sep = raw_beta(sep_accept, 1L),
    bhat2_sep = raw_beta(sep_accept, 2L),
    
    # Pooled estimated-fit diagnostics
    converged_pool = get_diag(pool_accept, "converged"),
    boundary_pool = get_diag(pool_accept, "boundary"),
    singular_pool = get_diag(pool_accept, "singular"),
    stalled_pool = get_diag(pool_accept, "stalled"),
    fallback_pool_accept = as.numeric(pool_accept$fallback),
    fallback_pool_fallback = as.numeric(pool_fallback$fallback),
    fallback_failed_pool_accept = as.numeric(pool_accept$fallback_failed),
    fallback_failed_pool_fallback = as.numeric(pool_fallback$fallback_failed),
    
    # Separate estimated-fit diagnostics
    converged_sep = get_diag(sep_accept, "converged"),
    boundary_sep = get_diag(sep_accept, "boundary"),
    singular_sep = get_diag(sep_accept, "singular"),
    stalled_sep = get_diag(sep_accept, "stalled"),
    fallback_sep_accept = as.numeric(sep_accept$fallback),
    fallback_sep_fallback = as.numeric(sep_fallback$fallback),
    fallback_failed_sep_accept = as.numeric(sep_accept$fallback_failed),
    fallback_failed_sep_fallback = as.numeric(sep_fallback$fallback_failed),
    
    # True-shape fit convergence (fixed-beta fits: no boundary, no fallback)
    converged_true_pool = if (is.null(m4p)) NA_real_ else as.numeric(m4p$converged),
    converged_true_sep = as.numeric(m4s$converged)
  )
  
  c(r, info)
}


run_cell <- function(s, n, theta, B = 1000, N = 10000, seed = 2026) {
  if (!is.numeric(theta) || length(theta) != 1L || is.na(theta) ||
      !is.finite(theta) || theta < 0) {
    stop("theta must be a single non-negative finite numeric value.")
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != as.integer(n)) {
    stop("n must be a single positive integer.")
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != as.integer(B)) {
    stop("B must be a single positive integer.")
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
      N < 2 || N != as.integer(N) || N %% 2 != 0) {
    stop("N must be a positive even integer (the total test-set size).")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed)) {
    stop("seed must be a single finite numeric value.")
  }
  if (!is.list(s) || is.null(s$p) || !is.numeric(s$p) ||
      length(s$p) != 1L || !is.finite(s$p) ||
      s$p < 1 || s$p != as.integer(s$p)) {
    stop("s must be a valid setup object returned by make_setup().")
  }

  set.seed(seed)
  res <- t(replicate(B, one_rep(s, n, theta, N)))
  as.data.frame(res)
}

## ---- 7. Summaries ----------------------------------------------------------
paired <- function(a, b, tau = 0.005) {
  
  # Retain only replications where both methods have valid errors.
  valid <- is.finite(a) & is.finite(b)
  d <- a[valid] - b[valid]
  
  n_valid <- length(d)
  
  if (n_valid == 0L) {
    return(c(
      n = 0,
      mean = NA_real_,
      mcse = NA_real_,
      lo = NA_real_,
      hi = NA_real_,
      verdict = "insufficient data"
    ))
  }
  
  m <- mean(d)
  
  se <- if (n_valid > 1L) {
    sd(d) / sqrt(n_valid)
  } else {
    NA_real_
  }
  
  lo <- m - 1.96 * se
  hi <- m + 1.96 * se
  
  verdict <- if (!is.finite(se)) {
    "insufficient data"
  } else if (hi < -tau) {
    "approach 1 better"
  } else if (lo > tau) {
    "approach 1 worse"
  } else if (lo > -tau && hi < tau) {
    "comparable"
  } else {
    "inconclusive"
  }
  
  c(
    n = n_valid,
    mean = m,
    mcse = se,
    lo = lo,
    hi = hi,
    verdict = verdict
  )
}

summarise_cell <- function(res, tau = 0.005) {
  
  rules <- c(
    "PE1_pool_accept", "PE1_pool_fallback",
    "PE1_sep_accept", "PE1_sep_fallback",
    "PE2_pool", "PE2_sep", "PE4_pool", "PE4_sep", "LDA", "QDA"
  )
  
  excess <- sapply(rules, function(m) {
    e <- res[[m]] - res$Bayes
    
    c(
      mean = mean(e, na.rm = TRUE),
      mcse = sd(e, na.rm = TRUE) /
        sqrt(sum(is.finite(e)))
    )
  })
  
  comparison_pairs <- list(
    RQ1_pool_accept = c("PE1_pool_accept", "PE2_pool"),
    RQ1_pool_fallback = c("PE1_pool_fallback", "PE2_pool"),
    RQ1_sep_accept = c("PE1_sep_accept", "PE2_sep"),
    RQ1_sep_fallback = c("PE1_sep_fallback", "PE2_sep"),
    RQ2_pool_accept = c("PE1_pool_accept", "LDA"),
    RQ2_pool_fallback = c("PE1_pool_fallback", "LDA"),
    RQ2_sep_accept = c("PE1_sep_accept", "QDA"),
    RQ2_sep_fallback = c("PE1_sep_fallback", "QDA"),
    check_pool = c("PE2_pool", "LDA"),
    check_sep = c("PE2_sep", "QDA"),
    shape_est_pool = c("PE1_pool_fallback", "PE4_pool"),
    shape_est_sep = c("PE1_sep_fallback", "PE4_sep"),
    fixed_true_pool = c("PE2_pool", "PE4_pool"),
    fixed_true_sep = c("PE2_sep", "PE4_sep")
  )
  
  paired_list <- lapply(comparison_pairs, function(pair) {
    paired(
      res[[pair[1]]],
      res[[pair[2]]],
      tau = tau
    )
  })
  
  paired_numeric <- do.call(rbind, lapply(paired_list, function(x) {
    as.numeric(x[c("n", "mean", "mcse", "lo", "hi")])
  }))
  
  colnames(paired_numeric) <- c(
    "n", "mean", "mcse", "lo", "hi"
  )
  
  paired_verdict <- vapply(
    paired_list,
    function(x) as.character(x["verdict"]),
    character(1)
  )
  
  paired_results <- data.frame(
    paired_numeric,
    verdict = paired_verdict,
    row.names = names(comparison_pairs),
    check.names = FALSE
  )
  
  fail <- c(
    converged_pool = mean(res$converged_pool, na.rm = TRUE),
    boundary_pool = mean(res$boundary_pool, na.rm = TRUE),
    singular_pool = mean(res$singular_pool, na.rm = TRUE),
    stalled_pool = mean(res$stalled_pool, na.rm = TRUE),
    
    converged_sep = mean(res$converged_sep, na.rm = TRUE),
    boundary_sep = mean(res$boundary_sep, na.rm = TRUE),
    singular_sep = mean(res$singular_sep, na.rm = TRUE),
    stalled_sep = mean(res$stalled_sep, na.rm = TRUE),
    
    fallback_pool_accept = mean(res$fallback_pool_accept, na.rm = TRUE),
    fallback_pool_fallback = mean(res$fallback_pool_fallback, na.rm = TRUE),
    fallback_failed_pool_accept = mean(res$fallback_failed_pool_accept, na.rm = TRUE),
    fallback_failed_pool_fallback = mean(res$fallback_failed_pool_fallback, na.rm = TRUE),
    
    fallback_sep_accept = mean(res$fallback_sep_accept, na.rm = TRUE),
    fallback_sep_fallback = mean(res$fallback_sep_fallback, na.rm = TRUE),
    fallback_failed_sep_accept = mean(res$fallback_failed_sep_accept, na.rm = TRUE),
    fallback_failed_sep_fallback = mean(res$fallback_failed_sep_fallback, na.rm = TRUE)
  )
  
  list(
    excess = excess,
    paired = paired_results,
    fail = fail
  )
}

## ---- 8. Self-tests ---------------------------------------------------------
self_test <- function() {
  set.seed(1)

  # (a) Generator reproduces target covariance for several shape values.
  for (b in c(0.5, 1, 2)) {
    s <- make_setup(3, b, b, "equal")
    X <- rpe(4e5, rep(0, 3), s$S1, b)
    cat(sprintf(
      "beta=%.1f  c=%.4f  max|cov(X)-V| = %.4f\n",
      b, cfac(3, b), max(abs(cov(X) - s$V1))
    ))
  }

  # (b) Density integrates to 1 for p=2.
  for (b in c(0.5, 1, 2)) {
    Sig <- diag(2) / cfac(2, b)
    rad <- function(r) exp(lpe(cbind(r, 0), c(0, 0), Sig, b)) * 2 * pi * r
    cat(sprintf(
      "beta=%.1f  density integral (should be 1): %.5f\n",
      b, integrate(rad, 0, Inf)$value
    ))
  }

  # (c) Shape recovery: small diagnostic experiment, not a precision study.
  for (b in c(0.5, 1, 2)) {
    s <- make_setup(3, b, b, "equal")
    bh <- replicate(20, {
      d <- gen_data(s, 200, 200, 0)
      fit_core(d$X, d$y, est_beta = TRUE)$b
    })
    cat(sprintf("true beta=%.1f  mean est=%.3f\n", b, mean(bh, na.rm = TRUE)))
  }

  # (d) At beta=1, PE and matched Gaussian MLE decisions should agree.
  s <- make_setup(3, 1, 1, "unequal")
  d <- gen_data(s, 300, 300, 1.5)
  te <- gen_data(s, 5000, 5000, 1.5)

  m_sep <- fit_pe_rule(d$X, d$y, pooled = FALSE, est_beta = FALSE)
  g_sep <- fit_gaussian_mle(d$X, d$y, pooled = FALSE)
  m_pool <- fit_pe_rule(d$X, d$y, pooled = TRUE, est_beta = FALSE)
  g_pool <- fit_gaussian_mle(d$X, d$y, pooled = TRUE)

  cat("PE2_sep vs matched Gaussian MLE agreement:",
      mean(predict_pe(m_sep, te$X) == predict_gaussian_mle(g_sep, te$X)), "\n")
  cat("PE2_pool vs matched Gaussian MLE agreement:",
      mean(predict_pe(m_pool, te$X) == predict_gaussian_mle(g_pool, te$X)), "\n")

  invisible(TRUE)
}

## ---- 9. Demo ---------------------------------------------------------------
demo <- function(B = 50, N = 10000, M_cal = 2e5, M_validate = 1e6) {
  s <- make_setup(p = 2, b1 = 0.5, b2 = 0.5, cov_type = "equal")
  cal <- calibrate_theta(
    s, target = 0.15, M = M_cal, M_validate = M_validate
  )
  print(cal)

  if (!is.finite(cal$theta)) {
    stop("Calibration failed: target Bayes error was not bracketed.")
  }

  res <- run_cell(s, n = 100, theta = cal$theta, B = B, N = N)
  print(summarise_cell(res))
  invisible(list(calibration = cal, results = res))
}

# Suggested first checks:
# self_test()
# timing <- system.time(out <- demo(B = 50))
# print(timing)
# Larger pilot:
# out <- demo(B = 200)
