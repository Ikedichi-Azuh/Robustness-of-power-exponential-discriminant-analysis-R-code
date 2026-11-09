# ============================================================================
# make_calibration_tables.R : LaTeX tables of the calibrated location difference
# (theta) and the realised Bayes error, for Blocks A/B and Block C.
# Reads  grid_final/calibration.csv  and  grid_final_C/calibration.csv
# Writes calibration_tables.tex (an appendix section with two tables) and prints
# the largest and the median absolute deviation from target, to check the text.
# ============================================================================

IN_AB <- "grid_final/calibration.csv"
IN_C  <- "grid_final_C/calibration.csv"
OUT   <- "calibration_tables.tex"

fmt <- function(x, d) formatC(x, format = "f", digits = d)

read_cal <- function(file, kind) {
  cal <- read.csv(file, stringsAsFactors = FALSE)
  stopifnot(all(is.finite(cal$theta)), all(is.finite(cal$validated_bayes)))
  pr <- do.call(rbind, strsplit(cal$key, "|", fixed = TRUE))
  bayes <- cal$validated_bayes
  if (max(bayes) <= 1) bayes <- 100 * bayes            # work in per cent
  if (kind == "AB") {
    df <- data.frame(p = as.integer(pr[, 1]), b1 = as.numeric(pr[, 2]),
                     b2 = as.numeric(pr[, 3]), cov = pr[, 4],
                     target = as.numeric(pr[, 5]), stringsAsFactors = FALSE)
  } else {
    df <- data.frame(p = as.integer(pr[, 1]), dist = pr[, 2], cov = pr[, 3],
                     target = as.numeric(pr[, 4]), stringsAsFactors = FALSE)
  }
  df$theta <- cal$theta
  df$bayes <- bayes
  df
}

# one row per scenario, theta and Bayes error for the three targets side by side
wide <- function(df, idcols) {
  id <- do.call(paste, c(df[idcols], sep = "|"))
  tg <- sort(unique(df$target))
  stopifnot(length(tg) == 3)
  rows <- lapply(unique(id), function(k) {
    s <- df[id == k, ]
    s <- s[match(tg, s$target), ]
    stopifnot(nrow(s) == 3, !anyNA(s$theta))
    cbind(s[1, idcols, drop = FALSE],
          setNames(as.list(s$theta), paste0("theta_", 1:3)),
          setNames(as.list(s$bayes), paste0("bayes_", 1:3)))
  })
  do.call(rbind, rows)
}

cov_lab <- function(x) ifelse(x == "equal", "Equal", "Unequal")
num_cells <- function(r) paste(c(fmt(unlist(r[paste0("theta_", 1:3)]), 3),
                                 fmt(unlist(r[paste0("bayes_", 1:3)]), 2)),
                               collapse = " & ")
table_lines <- function(caption, label, first_col, body) c(
  "\\begin{table}[htbp]", "\\centering", "\\footnotesize",
  paste0("\\caption{", caption, "}"),
  paste0("\\label{", label, "}"),
  "\\setlength{\\tabcolsep}{3.5pt}",
  "\\begin{tabular}{llcrrrrrr}",
  "\\hline",
  paste0(first_col, " & Covariance & $p$ & \\multicolumn{3}{c}{$\\theta$ at target Bayes error} & ",
         "\\multicolumn{3}{c}{Realised Bayes error (\\%)} \\\\"),
  " & & & 5\\% & 15\\% & 30\\% & 5\\% & 15\\% & 30\\% \\\\ \\hline",
  body, "\\hline", "\\end{tabular}", "\\end{table}", "")

## ---- Blocks A and B -------------------------------------------------------------------
ab <- read_cal(IN_AB, "AB")
w  <- wide(ab, c("p", "b1", "b2", "cov"))
w  <- w[order(w$b1 != w$b2, w$b1, w$b2, w$cov != "equal", w$p), ]
shape_lab <- ifelse(w$b1 == w$b2,
                    paste0("Block A, $\\beta=", w$b1, "$"),
                    paste0("Block B, $(", w$b1, ",", w$b2, ")$"))
body_ab <- vapply(seq_len(nrow(w)), function(i) {
  paste0(shape_lab[i], " & ", cov_lab(w$cov[i]), " & ", w$p[i], " & ", num_cells(w[i, ]), " \\\\")
}, character(1))
lines_ab <- table_lines(
  paste0("Blocks~A and~B: calibrated location difference $\\theta$ and realised Bayes error ",
         "(validated on an independent sample of $10^6$ observations) for the three target Bayes errors. ",
         "Shape: common $\\beta$ (Block~A) or $(\\beta_1,\\beta_2)$ (Block~B). ",
         "Covariance: equal ($\\mathbf{V}_2=\\mathbf{V}_1$) or unequal ($\\mathbf{V}_2=1.25\\,\\mathbf{V}_1$)."),
  "tab:calibration_AB", "Scenario", body_ab)

## ---- Block C ---------------------------------------------------------------------------------
cc <- read_cal(IN_C, "C")
wc <- wide(cc, c("p", "dist", "cov"))
dist_lv  <- c("t3", "cn", "st", "t5", "t10", "sn")
dist_lab <- c(t3 = "Student $t$, 3 d.f.", cn = "Contaminated normal",
              st = "Skew-$t$, 5 d.f.", t5 = "Student $t$, 5 d.f.",
              t10 = "Student $t$, 10 d.f.", sn = "Skew-normal")
wc <- wc[order(match(wc$dist, dist_lv), wc$cov != "equal", wc$p), ]
body_c <- vapply(seq_len(nrow(wc)), function(i) {
  paste0(dist_lab[[wc$dist[i]]], " & ", cov_lab(wc$cov[i]), " & ", wc$p[i], " & ",
         num_cells(wc[i, ]), " \\\\")
}, character(1))
lines_c <- table_lines(
  paste0("Block~C: calibrated location difference $\\theta$ and realised Bayes error ",
         "(validated on an independent sample of $10^6$ observations) for the three target Bayes errors. ",
         "Covariance: equal ($\\mathbf{V}_2=\\mathbf{V}_1$) or unequal ($\\mathbf{V}_2=1.25\\,\\mathbf{V}_1$)."),
  "tab:calibration_C", "Distribution", body_c)

writeLines(c("\\section{Calibration of the class separation}", "\\label{app:calibration}", "",
             lines_ab, lines_c), OUT)
cat("written:", normalizePath(OUT), "\n\n")

dev_ab <- abs(ab$bayes - 100 * ab$target)
dev_c  <- abs(cc$bayes - 100 * cc$target)
cat("Blocks A/B:", nrow(ab), "targets; max |deviation| =", fmt(max(dev_ab), 3),
    "pp; median =", fmt(median(dev_ab), 3), "pp\n")
cat("Block C   :", nrow(cc), "targets; max |deviation| =", fmt(max(dev_c), 3),
    "pp; median =", fmt(median(dev_c), 3), "pp\n")
cat("(The text says: A/B within 0.09 pp, median 0.02; C within 0.12 pp, median 0.02.)\n")
