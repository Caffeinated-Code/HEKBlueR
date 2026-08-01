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
  df
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

normalize_responses <- function(cleaned) {
  plates <- split(cleaned, cleaned$plate_id)
  out <- lapply(plates, function(x) {
    neg <- safe_median(x$blank_corrected_od[x$control_type == "negative_control"])
    pos <- safe_median(x$blank_corrected_od[x$control_type %in% c("positive_control", "agonist_challenge_control")])
    agonist <- safe_median(x$blank_corrected_od[x$control_type == "agonist_challenge_control"])
    antagonist <- safe_median(x$blank_corrected_od[x$control_type == "known_antagonist_control"])
    denom_activation <- pos - neg
    denom_inhibition <- agonist - antagonist
    x$percent_activation <- 100 * (x$blank_corrected_od - neg) / denom_activation
    x$percent_inhibition <- 100 * (agonist - x$blank_corrected_od) / denom_inhibition
    x
  })
  do.call(rbind, out)
}

summarize_primary <- function(normalized) {
  test <- normalized[normalized$control_type == "test_sample", ]
  if (!nrow(test)) return(data.frame())
  agg <- aggregate(
    cbind(percent_activation, percent_inhibition, blank_corrected_od) ~ plate_id + assay_mode + peptide_id + concentration_uM,
    test,
    function(v) c(mean = safe_mean(v), sd = safe_sd(v), cv = cv_percent(v))
  )
  flatten <- function(mat, prefix) {
    out <- as.data.frame(mat)
    names(out) <- paste0(prefix, c("_mean", "_sd", "_cv"))
    out
  }
  result <- cbind(
    agg[c("plate_id", "assay_mode", "peptide_id", "concentration_uM")],
    flatten(agg$percent_activation, "activation"),
    flatten(agg$percent_inhibition, "inhibition"),
    flatten(agg$blank_corrected_od, "od")
  )
  result$primary_hit <- with(result, ifelse(
    assay_mode == "agonist" & activation_mean >= 50 & activation_cv <= 25,
    "AGONIST_HIT",
    ifelse(assay_mode == "antagonist" & inhibition_mean >= 50 & inhibition_cv <= 25, "ANTAGONIST_HIT", "NO_HIT")
  ))
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
  coef(fit)
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
    if (is.null(fit)) {
      rows[[i]] <- data.frame(key, bottom = NA, top = NA, ec50_ic50 = NA, hill = NA, curve_status = "FIT_FAILED")
      qc[[i]] <- data.frame(key, curve_qc_status = "FAIL", curve_flags = "FIT_FAILED")
      next
    }
    dose_min <- min(x$concentration_uM, na.rm = TRUE)
    dose_max <- max(x$concentration_uM, na.rm = TRUE)
    flags <- character()
    if (fit["ec50"] < dose_min || fit["ec50"] > dose_max) flags <- c(flags, "EC50_OUT_OF_RANGE")
    if (abs(fit["top"] - fit["bottom"]) < 30) flags <- c(flags, "WEAK_RESPONSE")
    if (abs(fit["hill"]) > 4 || abs(fit["hill"]) < 0.2) flags <- c(flags, "BAD_HILL_SLOPE")
    high_noise <- any(x[[paste0(ifelse(key$assay_mode == "antagonist", "inhibition", "activation"), "_cv")]] > 30, na.rm = TRUE)
    if (high_noise) flags <- c(flags, "HIGH_REPLICATE_NOISE")
    status <- ifelse(length(flags) == 0, "PASS", ifelse(any(flags %in% c("EC50_OUT_OF_RANGE", "WEAK_RESPONSE")), "WARN", "WARN"))
    rows[[i]] <- data.frame(
      key,
      bottom = round(fit["bottom"], 3),
      top = round(fit["top"], 3),
      ec50_ic50 = round(fit["ec50"], 4),
      hill = round(fit["hill"], 3),
      curve_status = ifelse(status == "PASS", "GOOD_CURVE", "REVIEW_CURVE"),
      stringsAsFactors = FALSE
    )
    qc[[i]] <- data.frame(key, curve_qc_status = status, curve_flags = ifelse(length(flags), paste(flags, collapse = ";"), "GOOD_CURVE"))
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
  design <- validate_design(plate_map, metadata)
  plate <- calculate_plate_qc(cleaned)
  normalized <- normalize_responses(cleaned)
  primary <- summarize_primary(normalized)
  dose <- fit_dose_response(primary)
  counter <- counter_assay_qc(normalized)
  final <- make_final_qc(design, plate, dose$qc, counter, primary)
  exclusions <- cleaned[cleaned$saturated_flag | cleaned$negative_corrected_flag | cleaned$missing_flag, ]
  results <- list(
    metadata_completeness = metadata_completeness(metadata),
    cleaned_well_data = cleaned,
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
