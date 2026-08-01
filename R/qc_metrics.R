safe_mean <- function(x) mean(x, na.rm = TRUE)
safe_sd <- function(x) stats::sd(x, na.rm = TRUE)
safe_median <- function(x) stats::median(x, na.rm = TRUE)
safe_mad <- function(x) stats::mad(x, constant = 1.4826, na.rm = TRUE)

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
  max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE)
}

col_bias_score <- function(df, value_col = "blank_corrected_od") {
  vals <- tapply(df[[value_col]], df$col, safe_median)
  if (length(vals) < 2) return(NA_real_)
  max(vals, na.rm = TRUE) - min(vals, na.rm = TRUE)
}

