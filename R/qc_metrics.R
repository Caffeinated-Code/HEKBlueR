safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

safe_mad <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::mad(x, constant = 1.4826)
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

cv_percent <- function(x) {
  m <- safe_mean(x)
  if (is.na(m) || abs(m) < .Machine$double.eps) return(NA_real_)
  100 * safe_sd(x) / abs(m)
}

z_prime <- function(pos, neg) {
  denom <- abs(safe_mean(pos) - safe_mean(neg))
  if (is.na(denom) || denom == 0) return(NA_real_)
  1 - (3 * (safe_sd(pos) + safe_sd(neg)) / denom)
}

robust_z_prime <- function(pos, neg) {
  denom <- abs(safe_median(pos) - safe_median(neg))
  if (is.na(denom) || denom == 0) return(NA_real_)
  1 - (3 * (safe_mad(pos) + safe_mad(neg)) / denom)
}

ssmd <- function(pos, neg) {
  denom <- sqrt(safe_sd(pos)^2 + safe_sd(neg)^2)
  if (is.na(denom) || denom == 0) return(NA_real_)
  (safe_mean(pos) - safe_mean(neg)) / denom
}

qc_status <- function(value, pass, warn = NULL, direction = c("gte", "lte")) {
  direction <- match.arg(direction)
  if (is.na(value)) return("FAIL")
  if (direction == "gte") {
    if (value >= pass) return("PASS")
    if (!is.null(warn) && value >= warn) return("WARN")
    return("FAIL")
  }
  if (value <= pass) return("PASS")
  if (!is.null(warn) && value <= warn) return("WARN")
  "FAIL"
}

edge_effect_score <- function(df, value_col = "blank_corrected_od") {
  edge <- df$row %in% c("A", "H") | df$col %in% c("01", "12")
  inner <- !edge
  if (!any(edge, na.rm = TRUE) || !any(inner, na.rm = TRUE)) return(NA_real_)
  abs(safe_median(df[[value_col]][edge]) - safe_median(df[[value_col]][inner]))
}

row_bias_score <- function(df, value_col = "blank_corrected_od") {
  vals <- tapply(df[[value_col]], df$row, safe_median)
  if (length(vals) < 2) return(NA_real_)
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2) return(NA_real_)
  max(vals) - min(vals)
}

col_bias_score <- function(df, value_col = "blank_corrected_od") {
  vals <- tapply(df[[value_col]], df$col, safe_median)
  if (length(vals) < 2) return(NA_real_)
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2) return(NA_real_)
  max(vals) - min(vals)
}
