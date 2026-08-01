source_if_needed <- function(path) {
  if (file.exists(path)) source(path)
}

clean_well_data <- function(raw_data) {
  df <- raw_data
  df$row <- substr(df$well, 1, 1)
  df$col <- substr(df$well, 2, 3)
  blanks <- aggregate(raw_od ~ plate_id, df[df$control_type == "blank", ], safe_median)
  names(blanks)[2] <- "blank_od"
  df <- merge(df, blanks, by = "plate_id", all.x = TRUE)
  df$blank_corrected_od <- df$raw_od - df$blank_od
  df$saturated_flag <- !is.na(df$raw_od) & df$raw_od >= 3
  df$negative_corrected_flag <- !is.na(df$blank_corrected_od) & df$blank_corrected_od < -0.02
  df$missing_flag <- is.na(df$raw_od)
  by_plate_med <- ave(df$blank_corrected_od, df$plate_id, FUN = function(v) safe_median(v))
  by_plate_mad <- ave(df$blank_corrected_od, df$plate_id, FUN = function(v) max(safe_mad(v), 0.001))
  df$robust_z <- (df$blank_corrected_od - by_plate_med) / by_plate_mad
  df$outlier_flag <- !is.na(df$robust_z) & abs(df$robust_z) > 4 & !df$control_type %in% c("positive_control", "agonist_challenge_control")
  df$cleaning_action <- ifelse(df$missing_flag, "MISSING_RAW_OD",
    ifelse(df$saturated_flag, "SATURATED_REVIEW",
      ifelse(df$negative_corrected_flag, "NEGATIVE_AFTER_BLANK_REVIEW",
        ifelse(df$outlier_flag, "ROBUST_OUTLIER_REVIEW", "KEEP")
      )
    )
  )
  df
}

input_eda <- function(raw_data, plate_map) {
  required_raw <- c("plate_id", "well", "assay_mode", "sample_id", "peptide_id", "control_type", "concentration_uM", "raw_od")
  missing_raw <- setdiff(required_raw, names(raw_data))
  plate_counts <- length(unique(raw_data$plate_id))
  well_counts <- nrow(raw_data)
  duplicate_wells <- sum(duplicated(paste(raw_data$plate_id, raw_data$well)))
  missing_od <- sum(is.na(raw_data$raw_od))
  saturated <- sum(!is.na(raw_data$raw_od) & raw_data$raw_od >= 3)
  controls <- paste(sort(unique(raw_data$control_type)), collapse = ";")
  data.frame(
    metric = c("plates", "wells", "missing_required_columns", "duplicate_plate_wells", "missing_raw_od", "saturated_od", "control_types"),
    value = c(plate_counts, well_counts, paste(missing_raw, collapse = ";"), duplicate_wells, missing_od, saturated, controls),
    status = c(
      ifelse(plate_counts >= 1, "PASS", "FAIL"),
      ifelse(well_counts > 0, "PASS", "FAIL"),
      ifelse(length(missing_raw) == 0, "PASS", "FAIL"),
      ifelse(duplicate_wells == 0, "PASS", "FAIL"),
      ifelse(missing_od == 0, "PASS", "WARN"),
      ifelse(saturated == 0, "PASS", "WARN"),
      ifelse(nchar(controls) > 0, "PASS", "FAIL")
    ),
    interpretation = c(
      "Number of plate IDs in the upload.",
      "Number of well-level rows in the upload.",
      "Columns needed for standard HEKBlueR analysis.",
      "A plate and well should appear once.",
      "Missing raw OD values cannot be interpreted without review.",
      "Very high OD values may exceed the useful reader range.",
      "Control labels found in the uploaded file."
    ),
    stringsAsFactors = FALSE
  )
}

metadata_completeness <- function(metadata) {
  required <- metadata[metadata$required %in% c(TRUE, "TRUE", "true", "1"), , drop = FALSE]
  optional <- metadata[!metadata$field %in% required$field, , drop = FALSE]
  required_complete <- sum(!is.na(required$value) & required$value != "")
  optional_complete <- sum(!is.na(optional$value) & optional$value != "")
  total <- nrow(metadata)
  complete <- sum(!is.na(metadata$value) & metadata$value != "")
  data.frame(
    metric = "metadata_completeness_percent",
    value = round(100 * complete / max(total, 1), 1),
    required_complete = required_complete,
    required_total = nrow(required),
    optional_complete = optional_complete,
    optional_total = nrow(optional),
    status = ifelse(required_complete == nrow(required) && complete / max(total, 1) >= 0.8, "PASS", ifelse(complete / max(total, 1) >= 0.5, "WARN", "FAIL")),
    stringsAsFactors = FALSE
  )
}

validate_design <- function(plate_map, metadata) {
  controls <- split(plate_map, plate_map$plate_id)
  out <- lapply(names(controls), function(pid) {
    x <- controls[[pid]]
    assay_mode <- unique(x$assay_mode)[1]
    required_controls <- switch(
      assay_mode,
      agonist = c("blank", "negative_control", "positive_control"),
      antagonist = c("blank", "negative_control", "agonist_challenge_control", "known_antagonist_control"),
      counter = c("blank", "negative_control", "positive_control", "viability_counter", "no_cell_interference", "unrelated_reporter", "null_cell_reporter"),
      c("blank", "negative_control")
    )
    present <- unique(x$control_type)
    missing_controls <- setdiff(required_controls, present)
    control_counts <- table(x$control_type)
    min_control_wells <- suppressWarnings(min(control_counts[required_controls], na.rm = TRUE))
    if (!is.finite(min_control_wells)) min_control_wells <- 0
    test_rows <- x[x$control_type == "test_sample", ]
    reps <- if (nrow(test_rows)) {
      aggregate(well ~ plate_id + sample_id + peptide_id + concentration_uM, test_rows, length)
    } else {
      data.frame()
    }
    min_tech <- if (nrow(reps)) min(reps$well, na.rm = TRUE) else NA_integer_
    dose_counts <- if (nrow(test_rows)) {
      aggregate(concentration_uM ~ peptide_id, test_rows, function(v) length(unique(v[!is.na(v)])))
    } else {
      data.frame()
    }
    min_doses <- if (nrow(dose_counts)) min(dose_counts$concentration_uM) else NA_integer_
    ipc_present <- "inter_plate_calibrator" %in% present
    status <- "PASS"
    reasons <- character()
    if (length(missing_controls)) {
      status <- "FAIL"
      reasons <- c(reasons, paste("missing controls:", paste(missing_controls, collapse = ",")))
    }
    if (!is.na(min_tech) && min_tech < 3) {
      status <- ifelse(status == "FAIL", "FAIL", "WARN")
      reasons <- c(reasons, "technical replicates below 3")
    }
    if (!is.na(min_doses) && assay_mode %in% c("agonist", "antagonist") && min_doses < 8) {
      status <- ifelse(min_doses < 5, "FAIL", ifelse(status == "FAIL", "FAIL", "WARN"))
      reasons <- c(reasons, "dose-response has fewer than 8 concentrations")
    }
    if (min_control_wells < 4) {
      status <- "FAIL"
      reasons <- c(reasons, "too few required control wells")
    } else if (min_control_wells < 8 && status != "FAIL") {
      status <- "WARN"
      reasons <- c(reasons, "control wells below preferred count of 8")
    }
    if (!ipc_present && assay_mode != "antagonist" && status != "FAIL") {
      status <- "WARN"
      reasons <- c(reasons, "inter-plate calibrator missing")
    }
    data.frame(
      plate_id = pid,
      assay_mode = assay_mode,
      design_status = status,
      missing_controls = paste(missing_controls, collapse = ";"),
      min_control_wells = min_control_wells,
      min_technical_replicates = min_tech,
      min_dose_points = min_doses,
      inter_plate_calibrator_present = ipc_present,
      notes = paste(reasons, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  design <- do.call(rbind, out)
  meta <- metadata_completeness(metadata)
  design$metadata_status <- meta$status
  design$metadata_completeness_percent <- meta$value
  design
}

calculate_plate_qc <- function(cleaned) {
  plates <- split(cleaned, cleaned$plate_id)
  out <- lapply(names(plates), function(pid) {
    x <- plates[[pid]]
    pos <- x$blank_corrected_od[x$control_type %in% c("positive_control", "agonist_challenge_control")]
    neg <- x$blank_corrected_od[x$control_type == "negative_control"]
    blank <- x$raw_od[x$control_type == "blank"]
    zp <- z_prime(pos, neg)
    rzp <- robust_z_prime(pos, neg)
    sb <- safe_mean(pos) / max(safe_mean(neg), .Machine$double.eps)
    signal_window <- safe_mean(pos) - safe_mean(neg)
    edge <- edge_effect_score(x)
    row_bias <- row_bias_score(x)
    col_bias <- col_bias_score(x)
    status <- "PASS"
    notes <- character()
    if (is.na(zp) || zp < 0.3) {
      status <- "FAIL"
      notes <- c(notes, "Z-prime below 0.3")
    } else if (zp < 0.5) {
      status <- "WARN"
      notes <- c(notes, "Z-prime below 0.5")
    }
    if (!is.na(cv_percent(pos)) && cv_percent(pos) > 20) {
      status <- ifelse(status == "FAIL", "FAIL", "WARN")
      notes <- c(notes, "positive control CV above 20%")
    }
    if (!is.na(cv_percent(neg)) && cv_percent(neg) > 20) {
      status <- ifelse(status == "FAIL", "FAIL", "WARN")
      notes <- c(notes, "negative control CV above 20%")
    }
    if (!is.na(edge) && edge > 0.15) {
      status <- ifelse(status == "FAIL", "FAIL", "WARN")
      notes <- c(notes, "edge effect detected")
    }
    data.frame(
      plate_id = pid,
      assay_mode = unique(x$assay_mode)[1],
      plate_qc_status = status,
      z_prime = round(zp, 3),
      robust_z_prime = round(rzp, 3),
      ssmd = round(ssmd(pos, neg), 3),
      signal_to_background = round(sb, 3),
      signal_window = round(signal_window, 3),
      positive_control_cv = round(cv_percent(pos), 2),
      negative_control_cv = round(cv_percent(neg), 2),
      blank_cv = round(cv_percent(blank), 2),
      edge_effect = round(edge, 3),
      row_bias = round(row_bias, 3),
      column_bias = round(col_bias, 3),
      saturated_wells = sum(x$saturated_flag, na.rm = TRUE),
      missing_wells = sum(x$missing_flag, na.rm = TRUE),
      notes = paste(notes, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

reference_control_qc <- function(cleaned) {
  control_types <- c("negative_control", "positive_control", "agonist_challenge_control", "known_antagonist_control", "inter_plate_calibrator")
  x <- cleaned[cleaned$control_type %in% control_types, ]
  if (!nrow(x)) return(data.frame())
  agg <- aggregate(blank_corrected_od ~ plate_id + control_type, x, function(v) c(mean = safe_mean(v), median = safe_median(v), cv = cv_percent(v), n = sum(!is.na(v))))
  out <- cbind(agg[c("plate_id", "control_type")], as.data.frame(agg$blank_corrected_od))
  names(out)[3:6] <- c("mean_od", "median_od", "cv_percent", "n_wells")
  out$status <- ifelse(out$n_wells < 4, "FAIL", ifelse(out$cv_percent > 20, "WARN", "PASS"))
  out$interpretation <- ifelse(out$status == "PASS", "Control is stable enough for normalization.",
    ifelse(out$status == "WARN", "Control variability is high. Review plate handling and pipetting.", "Too few control wells for reliable normalization.")
  )
  out
}

interplate_calibration <- function(cleaned) {
  calibrator <- cleaned[cleaned$control_type == "inter_plate_calibrator", ]
  method <- "inter_plate_calibrator"
  if (!nrow(calibrator)) {
    calibrator <- cleaned[cleaned$control_type %in% c("positive_control", "agonist_challenge_control"), ]
    method <- "shared_positive_control"
  }
  if (!nrow(calibrator)) {
    return(list(
      factors = data.frame(plate_id = unique(cleaned$plate_id), calibration_method = "none", calibration_factor = 1, calibration_status = "WARN", notes = "No shared calibrator or positive control found."),
      calibrated = transform(cleaned, calibrated_od = blank_corrected_od)
    ))
  }
  med <- aggregate(blank_corrected_od ~ plate_id, calibrator, safe_median)
  global <- safe_median(med$blank_corrected_od)
  med$calibration_method <- method
  med$calibration_factor <- global - med$blank_corrected_od
  med$calibrator_median_od <- med$blank_corrected_od
  med$calibration_status <- ifelse(abs(med$calibration_factor) <= 0.15, "PASS", "WARN")
  med$notes <- ifelse(med$calibration_status == "PASS", "Plate is aligned with shared reference.", "Plate shows calibration drift above 0.15 OD.")
  med$blank_corrected_od <- NULL
  calibrated <- merge(cleaned, med[c("plate_id", "calibration_factor")], by = "plate_id", all.x = TRUE)
  calibrated$calibration_factor[is.na(calibrated$calibration_factor)] <- 0
  calibrated$calibrated_od <- calibrated$blank_corrected_od + calibrated$calibration_factor
  list(factors = med, calibrated = calibrated)
}

normalize_responses <- function(cleaned) {
  plates <- split(cleaned, cleaned$plate_id)
  out <- lapply(plates, function(x) {
    value_col <- if ("calibrated_od" %in% names(x)) "calibrated_od" else "blank_corrected_od"
    neg <- safe_median(x[[value_col]][x$control_type == "negative_control"])
    pos <- safe_median(x[[value_col]][x$control_type %in% c("positive_control", "agonist_challenge_control")])
    agonist <- safe_median(x[[value_col]][x$control_type == "agonist_challenge_control"])
    antagonist <- safe_median(x[[value_col]][x$control_type == "known_antagonist_control"])
    denom_activation <- pos - neg
    denom_inhibition <- agonist - antagonist
    x$percent_activation <- 100 * (x[[value_col]] - neg) / denom_activation
    x$percent_inhibition <- 100 * (agonist - x[[value_col]]) / denom_inhibition
    x
  })
  do.call(rbind, out)
}

summarize_primary <- function(normalized) {
  test <- normalized[normalized$control_type == "test_sample", ]
  if (!nrow(test)) return(data.frame())
  keys <- unique(test[c("plate_id", "assay_mode", "peptide_id", "concentration_uM")])
  result <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    x <- test[
      test$plate_id == key$plate_id &
        test$assay_mode == key$assay_mode &
        test$peptide_id == key$peptide_id &
        test$concentration_uM == key$concentration_uM,
    ]
    data.frame(
      key,
      activation_mean = safe_mean(x$percent_activation),
      activation_sd = safe_sd(x$percent_activation),
      activation_cv = cv_percent(x$percent_activation),
      inhibition_mean = safe_mean(x$percent_inhibition),
      inhibition_sd = safe_sd(x$percent_inhibition),
      inhibition_cv = cv_percent(x$percent_inhibition),
      od_mean = safe_mean(x$blank_corrected_od),
      od_sd = safe_sd(x$blank_corrected_od),
      od_cv = cv_percent(x$blank_corrected_od),
      stringsAsFactors = FALSE
    )
  }))
  result$primary_hit <- with(result, ifelse(
    assay_mode == "agonist" & activation_mean >= 50 & activation_cv <= 25,
    "AGONIST_HIT",
    ifelse(assay_mode == "antagonist" & inhibition_mean >= 50 & inhibition_cv <= 25, "ANTAGONIST_HIT", "NO_HIT")
  ))
  result$primary_status <- ifelse(result$primary_hit == "NO_HIT", "PASS", "WARN")
  result
}

fit_one_curve <- function(df, response_col) {
  df <- df[is.finite(df$concentration_uM) & is.finite(df[[response_col]]) & df$concentration_uM > 0, ]
  if (length(unique(df$concentration_uM)) < 5) return(NULL)
  y <- df[[response_col]]
  x <- df$concentration_uM
  start <- list(bottom = min(y), top = max(y), ec50 = stats::median(x), hill = 1)
  fit <- try(stats::nls(
    y ~ bottom + (top - bottom) / (1 + (ec50 / x)^hill),
    start = start,
    control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
  ), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  pred <- stats::predict(fit)
  list(coef = coef(fit), rmse = sqrt(mean((y - pred)^2, na.rm = TRUE)), residual_max = max(abs(y - pred), na.rm = TRUE))
}

fit_dose_response <- function(primary) {
  if (!nrow(primary)) return(list(results = data.frame(), qc = data.frame()))
  keys <- unique(primary[c("plate_id", "assay_mode", "peptide_id")])
  rows <- list()
  qc <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    x <- primary[primary$plate_id == key$plate_id & primary$peptide_id == key$peptide_id, ]
    response_col <- if (key$assay_mode == "antagonist") "inhibition_mean" else "activation_mean"
    fit <- fit_one_curve(x, response_col)
    response <- x[[response_col]]
    n_dose_points <- length(unique(x$concentration_uM[is.finite(x$concentration_uM)]))
    response_min <- min(response, na.rm = TRUE)
    response_max <- max(response, na.rm = TRUE)
    dynamic_range <- response_max - response_min
    ordered_x <- x[order(x$concentration_uM), ]
    ordered_response <- ordered_x[[response_col]]
    monotonic_violations <- if (key$assay_mode == "antagonist") {
      sum(diff(ordered_response) < -10, na.rm = TRUE)
    } else {
      sum(diff(ordered_response) < -10, na.rm = TRUE)
    }
    max_rep_cv <- max(x[[paste0(ifelse(key$assay_mode == "antagonist", "inhibition", "activation"), "_cv")]], na.rm = TRUE)
    if (!is.finite(max_rep_cv)) max_rep_cv <- NA_real_
    if (is.null(fit)) {
      rows[[i]] <- data.frame(key, n_dose_points = n_dose_points, bottom = NA, top = NA, ec50_ic50 = NA, hill = NA, response_min = round(response_min, 4), response_max = round(response_max, 4), dynamic_range = round(dynamic_range, 4), rmse = NA, max_residual = NA, curve_status = "FIT_FAILED")
      qc[[i]] <- data.frame(key, curve_qc_status = "FAIL", n_dose_points = n_dose_points, max_replicate_cv = round(max_rep_cv, 4), monotonic_violations = monotonic_violations, dynamic_range = round(dynamic_range, 4), ec50_in_range = FALSE, top_plateau_observed = FALSE, bottom_plateau_observed = FALSE, curve_flags = "FIT_FAILED")
      next
    }
    coef <- fit$coef
    dose_min <- min(x$concentration_uM, na.rm = TRUE)
    dose_max <- max(x$concentration_uM, na.rm = TRUE)
    flags <- character()
    ec50_in_range <- coef["ec50"] >= dose_min && coef["ec50"] <= dose_max
    top_plateau <- abs(response_max - coef["top"]) <= 20
    bottom_plateau <- abs(response_min - coef["bottom"]) <= 20
    if (n_dose_points < 8) flags <- c(flags, "LOW_DOSE_COUNT")
    if (!ec50_in_range) flags <- c(flags, "EC50_OUT_OF_RANGE")
    if (abs(coef["top"] - coef["bottom"]) < 30) flags <- c(flags, "WEAK_RESPONSE")
    if (!top_plateau) flags <- c(flags, "NO_TOP_PLATEAU")
    if (!bottom_plateau) flags <- c(flags, "NO_BOTTOM_PLATEAU")
    if (abs(coef["hill"]) > 4 || abs(coef["hill"]) < 0.2) flags <- c(flags, "BAD_HILL_SLOPE")
    if (monotonic_violations > 1) flags <- c(flags, "NON_MONOTONIC")
    high_noise <- !is.na(max_rep_cv) && max_rep_cv > 30
    if (high_noise) flags <- c(flags, "HIGH_REPLICATE_NOISE")
    status <- ifelse(length(flags) == 0, "PASS", ifelse(any(flags %in% c("FIT_FAILED", "LOW_DOSE_COUNT", "EC50_OUT_OF_RANGE", "WEAK_RESPONSE")), "WARN", "WARN"))
    rows[[i]] <- data.frame(
      key,
      n_dose_points = n_dose_points,
      bottom = round(coef["bottom"], 4),
      top = round(coef["top"], 4),
      ec50_ic50 = round(coef["ec50"], 4),
      hill = round(coef["hill"], 4),
      response_min = round(response_min, 4),
      response_max = round(response_max, 4),
      dynamic_range = round(dynamic_range, 4),
      rmse = round(fit$rmse, 4),
      max_residual = round(fit$residual_max, 4),
      curve_status = ifelse(status == "PASS", "GOOD_CURVE", "REVIEW_CURVE"),
      stringsAsFactors = FALSE
    )
    qc[[i]] <- data.frame(key, curve_qc_status = status, n_dose_points = n_dose_points, max_replicate_cv = round(max_rep_cv, 4), monotonic_violations = monotonic_violations, dynamic_range = round(dynamic_range, 4), ec50_in_range = ec50_in_range, top_plateau_observed = top_plateau, bottom_plateau_observed = bottom_plateau, curve_flags = ifelse(length(flags), paste(flags, collapse = ";"), "GOOD_CURVE"))
  }
  list(results = do.call(rbind, rows), qc = do.call(rbind, qc))
}

counter_assay_qc <- function(normalized) {
  x <- normalized[normalized$assay_mode == "counter" & !is.na(normalized$peptide_id), ]
  if (!nrow(x)) return(data.frame())
  peptides <- unique(x$peptide_id)
  out <- lapply(peptides, function(p) {
    z <- x[x$peptide_id == p, ]
    viability <- safe_mean(z$raw_od[z$control_type == "viability_counter"])
    no_cell <- safe_mean(z$raw_od[z$control_type == "no_cell_interference"])
    unrelated <- safe_mean(z$raw_od[z$control_type == "unrelated_reporter"])
    null_cell <- safe_mean(z$raw_od[z$control_type == "null_cell_reporter"])
    flags <- character()
    if (!is.na(viability) && viability < 0.7) flags <- c(flags, "CYTOTOXICITY_CONFOUNDED")
    if (!is.na(no_cell) && no_cell > 0.25) flags <- c(flags, "ASSAY_INTERFERENCE")
    if (!is.na(unrelated) && unrelated > 0.4) flags <- c(flags, "UNRELATED_REPORTER_ACTIVITY")
    if (!is.na(null_cell) && null_cell > 0.35) flags <- c(flags, "NULL_CELL_ACTIVITY")
    data.frame(
      peptide_id = p,
      counter_qc_status = ifelse(length(flags), "WARN", "PASS"),
      viability_signal = round(viability, 3),
      no_cell_signal = round(no_cell, 3),
      unrelated_reporter_signal = round(unrelated, 3),
      null_cell_signal = round(null_cell, 3),
      artifact_flags = ifelse(length(flags), paste(flags, collapse = ";"), "NONE"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

make_final_qc <- function(design_qc, plate_qc, dose_qc, counter_qc, primary) {
  peptides <- unique(primary$peptide_id)
  out <- lapply(peptides, function(p) {
    p_primary <- primary[primary$peptide_id == p, ]
    true_hits <- unique(p_primary$primary_hit[p_primary$primary_hit != "NO_HIT"])
    hit <- ifelse(length(true_hits), paste(true_hits, collapse = ";"), "NO_HIT")
    p_curve <- dose_qc[dose_qc$peptide_id == p, , drop = FALSE]
    p_counter <- counter_qc[counter_qc$peptide_id == p, , drop = FALSE]
    artifact <- if (nrow(p_counter)) p_counter$artifact_flags[1] else "NOT_TESTED"
    final <- "PASS_ADVANCE"
    if (hit == "NO_HIT") final <- "NO_ADVANCE"
    if (nrow(p_curve) && any(p_curve$curve_qc_status == "FAIL")) final <- "REVIEW_CURVE_NOISE"
    if (!artifact %in% c("NONE", "NOT_TESTED")) final <- "LIKELY_ARTIFACT"
    data.frame(
      final_status = ifelse(final %in% c("PASS_ADVANCE", "NO_ADVANCE"), "PASS", ifelse(final == "LIKELY_ARTIFACT", "FAIL", "WARN")),
      peptide_id = p,
      hit_call = hit,
      curve_qc = ifelse(nrow(p_curve), paste(unique(p_curve$curve_qc_status), collapse = ";"), "NOT_RUN"),
      counter_assay_qc = ifelse(nrow(p_counter), p_counter$counter_qc_status[1], "NOT_RUN"),
      artifact_flags = artifact,
      final_action = final,
      stringsAsFactors = FALSE
    )
  })
  plate_summary <- data.frame(
    design_overall = ifelse(any(design_qc$design_status == "FAIL"), "FAIL", ifelse(any(design_qc$design_status == "WARN"), "WARN", "PASS")),
    plate_overall = ifelse(any(plate_qc$plate_qc_status == "FAIL"), "FAIL", ifelse(any(plate_qc$plate_qc_status == "WARN"), "WARN", "PASS")),
    stringsAsFactors = FALSE
  )
  final <- do.call(rbind, out)
  cbind(plate_summary[rep(1, nrow(final)), ], final)
}

run_hekblue_analysis <- function(raw_data, plate_map, metadata, output_dir = NULL) {
  cleaned <- clean_well_data(raw_data)
  eda <- input_eda(raw_data, plate_map)
  design <- validate_design(plate_map, metadata)
  plate <- calculate_plate_qc(cleaned)
  ref_qc <- reference_control_qc(cleaned)
  calibration <- interplate_calibration(cleaned)
  normalized <- normalize_responses(calibration$calibrated)
  primary <- summarize_primary(normalized)
  dose <- fit_dose_response(primary)
  counter <- counter_assay_qc(normalized)
  final <- make_final_qc(design, plate, dose$qc, counter, primary)
  exclusions <- cleaned[cleaned$saturated_flag | cleaned$negative_corrected_flag | cleaned$missing_flag, ]
  results <- list(
    metadata_completeness = metadata_completeness(metadata),
    input_eda = eda,
    cleaned_well_data = cleaned,
    interplate_calibration = calibration$factors,
    reference_control_qc = ref_qc,
    normalized_results = normalized,
    design_qc = design,
    plate_qc = plate,
    primary_results = primary,
    dose_response_results = dose$results,
    dose_response_qc = dose$qc,
    counter_assay_qc = counter,
    hit_calls = primary[primary$primary_hit != "NO_HIT", ],
    exclusions = exclusions,
    final_qc_table = final
  )
  if (!is.null(output_dir)) export_hekblue_results(results, raw_data, plate_map, metadata, output_dir)
  results
}
