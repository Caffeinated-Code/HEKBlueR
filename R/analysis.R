source_if_needed <- function(path) {
  if (file.exists(path)) source(path)
}

enrich_raw_with_plate_map <- function(raw_data, plate_map) {
  join_cols <- intersect(c("plate_id", "well"), intersect(names(raw_data), names(plate_map)))
  if (!all(c("plate_id", "well") %in% join_cols)) return(raw_data)
  extra_cols <- setdiff(
    intersect(c("target_id", "biological_replicate", "technical_replicate", "expected_activity", "compound_role", "cpd_role"), names(plate_map)),
    names(raw_data)
  )
  if (!length(extra_cols)) return(raw_data)
  merge(raw_data, plate_map[c(join_cols, extra_cols)], by = join_cols, all.x = TRUE, sort = FALSE)
}

object_md5 <- function(x) {
  tmp <- tempfile(fileext = ".rds")
  saveRDS(x, tmp, version = 2)
  on.exit(unlink(tmp), add = TRUE)
  unname(tools::md5sum(tmp))
}

metadata_value <- function(metadata, field, default = "") {
  if (!all(c("field", "value") %in% names(metadata))) return(default)
  value <- metadata$value[match(field, metadata$field)]
  ifelse(is.na(value) || value == "", default, as.character(value))
}

safe_id_part <- function(x, fallback = "NA", max_chars = 28) {
  x <- toupper(gsub("[^A-Za-z0-9]+", "-", as.character(x)))
  x <- gsub("^-|-$", "", x)
  x <- ifelse(is.na(x) || x == "", fallback, x)
  substr(x, 1, max_chars)
}

hash_to_6_digits <- function(hash) {
  sprintf("%06d", as.integer(strtoi(substr(hash, 1, 7), base = 16L)) %% 1000000)
}

thresholds_equal <- function(x, y, tolerance = 1e-9) {
  if (is.list(x) || is.list(y)) {
    if (!is.list(x) || !is.list(y)) return(FALSE)
    if (!identical(sort(names(x)), sort(names(y)))) return(FALSE)
    return(all(vapply(sort(names(x)), function(nm) thresholds_equal(x[[nm]], y[[nm]], tolerance), logical(1))))
  }
  if (length(x) != length(y)) return(FALSE)
  if (is.numeric(x) || is.numeric(y)) {
    x_num <- as.numeric(x)
    y_num <- as.numeric(y)
    if (!identical(is.na(x_num), is.na(y_num))) return(FALSE)
    return(all(abs(x_num[!is.na(x_num)] - y_num[!is.na(y_num)]) <= tolerance))
  }
  identical(as.character(x), as.character(y))
}

thresholds_changed_from_default <- function(thresholds = default_qc_thresholds()) {
  !thresholds_equal(thresholds, default_qc_thresholds())
}

assay_manifest <- function(raw_data, plate_map, metadata, thresholds = default_qc_thresholds(), threshold_change_note = "") {
  raw_sig <- object_md5(raw_data)
  plate_sig <- object_md5(plate_map)
  meta_sig <- object_md5(metadata)
  threshold_sig <- object_md5(thresholds)
  strategy <- "HEKBlueR plate QC, normalization, primary summary, dose-response, counter-assay review"
  combined <- object_md5(list(raw = raw_sig, plate_map = plate_sig, metadata = meta_sig, thresholds = threshold_sig, threshold_change_note = trimws(threshold_change_note), strategy = strategy))
  project <- metadata_value(metadata, "project")
  target <- metadata_value(metadata, "target_id", "UNKNOWN_TARGET")
  assay_date <- metadata_value(metadata, "assay_date", format(Sys.Date(), "%Y-%m-%d"))
  scientist <- metadata_value(metadata, "scientist", "UNKNOWN_PERSON")
  assay_types <- if ("assay_mode" %in% names(raw_data)) paste(sort(unique(raw_data$assay_mode)), collapse = "+") else "unknown"
  assay_type_id <- if ("assay_mode" %in% names(raw_data)) {
    modes <- sort(unique(raw_data$assay_mode))
    paste(vapply(modes, function(mode) switch(mode, agonist = "AGO", antagonist = "ANT", counter = "CTR", unknown = "UNK", safe_id_part(mode, "ASSAY", 6)), character(1)), collapse = "+")
  } else {
    "UNK"
  }
  peptides <- if ("peptide_id" %in% names(raw_data)) sort(unique(raw_data$peptide_id[!is.na(raw_data$peptide_id) & raw_data$peptide_id != ""])) else character()
  peptide_label <- if (length(peptides)) paste(head(peptides, 3), collapse = "+") else "NO_PEPTIDE"
  if (length(peptides) > 3) peptide_label <- paste0(peptide_label, "+N", length(peptides))
  change_code <- hash_to_6_digits(combined)
  assay_id <- paste(
    safe_id_part("HEKBLUER", max_chars = 8),
    safe_id_part(target, "TARGET"),
    safe_id_part(assay_type_id, "ASSAY", 18),
    safe_id_part(peptide_label, "PEPTIDE", 32),
    safe_id_part(assay_date, "DATE", 10),
    safe_id_part(scientist, "PERSON", 18),
    change_code,
    sep = "-"
  )
  data.frame(
    assay_identifier = assay_id,
    change_code = change_code,
    input_signature = combined,
    raw_data_signature = raw_sig,
    plate_map_signature = plate_sig,
    metadata_signature = meta_sig,
    threshold_signature = threshold_sig,
    threshold_changed_from_default = thresholds_changed_from_default(thresholds),
    threshold_change_note = trimws(threshold_change_note),
    project = project,
    target_id = target,
    assay_type = assay_types,
    peptide_ids = paste(peptides, collapse = ";"),
    scientist = scientist,
    assay_date = assay_date,
    analysis_strategy = strategy,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    stringsAsFactors = FALSE
  )
}

run_documentation_table <- function(manifest, metadata, thresholds = default_qc_thresholds(), threshold_change_note = "") {
  threshold_rows <- qc_threshold_table(thresholds)
  threshold_doc <- data.frame(
    assay_identifier = manifest$assay_identifier[1],
    section = "qc_threshold",
    parameter = threshold_rows$metric,
    value = threshold_rows$threshold,
    source = threshold_rows$qc_layer,
    notes = threshold_rows$interpretation,
    stringsAsFactors = FALSE
  )
  metadata_doc <- data.frame(
    assay_identifier = manifest$assay_identifier[1],
    section = "metadata",
    parameter = metadata$field,
    value = as.character(metadata$value),
    source = ifelse(metadata$required %in% c(TRUE, "TRUE", "true", "1"), "required", "optional"),
    notes = "",
    stringsAsFactors = FALSE
  )
  manifest_doc <- data.frame(
    assay_identifier = manifest$assay_identifier[1],
    section = "assay_manifest",
    parameter = names(manifest),
    value = as.character(unlist(manifest[1, ], use.names = FALSE)),
    source = "computed",
    notes = vapply(names(manifest), function(name) switch(name,
      assay_identifier = "Readable deterministic run identifier.",
      change_code = "Six-digit code that changes when inputs, metadata, thresholds, or analysis strategy change.",
      input_signature = "Combined signature of raw data, plate map, metadata, thresholds, and analysis strategy.",
      raw_data_signature = "Raw data signature.",
      plate_map_signature = "Plate map signature.",
      metadata_signature = "Metadata signature.",
      threshold_signature = "QC threshold signature.",
      threshold_changed_from_default = "TRUE when active QC thresholds differ from the defaults.",
      threshold_change_note = "Required rationale when active QC thresholds differ from the defaults.",
      project = "Project metadata.",
      target_id = "Target metadata.",
      assay_type = "Assay modules found in the raw data.",
      peptide_ids = "Peptide or compound IDs found in the raw data.",
      scientist = "Scientist or analyst supplied in metadata.",
      assay_date = "Assay date supplied in metadata.",
      analysis_strategy = "Analysis workflow included in the input signature.",
      created_at = "Run timestamp.",
      "Computed assay manifest field."
    ), character(1)),
    stringsAsFactors = FALSE
  )
  threshold_note_doc <- data.frame(
    assay_identifier = manifest$assay_identifier[1],
    section = "analysis_parameter",
    parameter = "threshold_change_note",
    value = ifelse(trimws(threshold_change_note) == "", manifest$threshold_change_note[1], trimws(threshold_change_note)),
    source = "scientist-entered",
    notes = "Required when QC thresholds differ from the defaults.",
    stringsAsFactors = FALSE
  )
  config_doc <- data.frame(
    assay_identifier = manifest$assay_identifier[1],
    section = "analysis_parameter",
    parameter = c("assay_type", "normalization", "curve_model", "table_precision_export", "table_precision_app"),
    value = c("HEK-Blue SEAP reporter", "blank correction plus control scaling", "four-parameter logistic with base R nls", "4 decimal places for analysis exports", "2 decimal places where decimals are useful"),
    source = "HEKBlueR",
    notes = c(
      "Reporter assay readout.",
      "Applied before percent activation and inhibition.",
      "Used for dose-response summaries when enough dose points are available.",
      "Keeps downstream analysis and database loads precise.",
      "Keeps visual review readable."
    ),
    stringsAsFactors = FALSE
  )
  rbind(manifest_doc, metadata_doc, threshold_doc, threshold_note_doc, config_doc)
}

fill_target_id <- function(df, metadata) {
  metadata_target <- if (all(c("field", "value") %in% names(metadata))) {
    metadata$value[match("target_id", metadata$field)]
  } else {
    NA_character_
  }
  if (!"target_id" %in% names(df)) df$target_id <- NA_character_
  df$target_id[is.na(df$target_id) | df$target_id == ""] <- ifelse(
    is.na(metadata_target) || metadata_target == "",
    "UNKNOWN_TARGET",
    metadata_target
  )
  df
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

raw_data_summary <- function(raw_data) {
  assay_modes <- if ("assay_mode" %in% names(raw_data)) paste(sort(unique(raw_data$assay_mode)), collapse = ";") else ""
  data.frame(
    metric = c(
      "rows",
      "plates",
      "assay_modes",
      "test_samples",
      "control_types",
      "concentration_points",
      "technical_replicate_values",
      "biological_replicate_values",
      "raw_od_min",
      "raw_od_median",
      "raw_od_max"
    ),
    value = c(
      nrow(raw_data),
      if ("plate_id" %in% names(raw_data)) length(unique(raw_data$plate_id)) else NA_integer_,
      assay_modes,
      if ("peptide_id" %in% names(raw_data)) length(unique(raw_data$peptide_id[!is.na(raw_data$peptide_id) & raw_data$peptide_id != ""])) else NA_integer_,
      if ("control_type" %in% names(raw_data)) length(unique(raw_data$control_type)) else NA_integer_,
      if ("concentration_uM" %in% names(raw_data)) length(unique(raw_data$concentration_uM[!is.na(raw_data$concentration_uM)])) else NA_integer_,
      if ("technical_replicate" %in% names(raw_data)) paste(sort(unique(raw_data$technical_replicate[!is.na(raw_data$technical_replicate)])), collapse = ";") else "",
      if ("biological_replicate" %in% names(raw_data)) paste(sort(unique(raw_data$biological_replicate[!is.na(raw_data$biological_replicate)])), collapse = ";") else "",
      if ("raw_od" %in% names(raw_data)) round(min(raw_data$raw_od, na.rm = TRUE), 4) else NA_real_,
      if ("raw_od" %in% names(raw_data)) round(safe_median(raw_data$raw_od), 4) else NA_real_,
      if ("raw_od" %in% names(raw_data)) round(max(raw_data$raw_od, na.rm = TRUE), 4) else NA_real_
    ),
    interpretation = c(
      "Total uploaded well-level records.",
      "Number of unique plates in the upload.",
      "Assay modules detected in the upload.",
      "Number of unique tested peptide or compound IDs.",
      "Number of well roles detected.",
      "Number of non-missing dose levels.",
      "Technical replicate labels detected.",
      "Biological replicate labels detected.",
      "Lowest raw optical density.",
      "Median raw optical density.",
      "Highest raw optical density."
    ),
    stringsAsFactors = FALSE
  )
}

metadata_summary <- function(metadata) {
  required <- metadata[metadata$required %in% c(TRUE, "TRUE", "true", "1"), , drop = FALSE]
  optional <- metadata[!(metadata$required %in% c(TRUE, "TRUE", "true", "1")), , drop = FALSE]
  data.frame(
    metric = c("metadata_fields", "required_fields", "optional_fields", "blank_values", "required_blank_values", "key_fields_present"),
    value = c(
      nrow(metadata),
      nrow(required),
      nrow(optional),
      sum(is.na(metadata$value) | metadata$value == ""),
      sum(is.na(required$value) | required$value == ""),
      paste(intersect(c("scientist", "project", "assay_date", "target_id", "cell_line", "protocol_version"), metadata$field[!is.na(metadata$value) & metadata$value != ""]), collapse = ";")
    ),
    interpretation = c(
      "Total metadata rows supplied.",
      "Fields treated as required for reproducibility.",
      "Fields useful for audit, troubleshooting, and search.",
      "Metadata rows without values.",
      "Required metadata rows without values.",
      "High-value documentation fields found."
    ),
    stringsAsFactors = FALSE
  )
}

metadata_completeness <- function(metadata, thresholds = default_qc_thresholds()) {
  required <- metadata[metadata$required %in% c(TRUE, "TRUE", "true", "1"), , drop = FALSE]
  optional <- metadata[!metadata$field %in% required$field, , drop = FALSE]
  required_complete <- sum(!is.na(required$value) & required$value != "")
  optional_complete <- sum(!is.na(optional$value) & optional$value != "")
  total <- nrow(metadata)
  complete <- sum(!is.na(metadata$value) & metadata$value != "")
  pass_cutoff <- thresholds$metadata$completeness_pass_percent
  warn_cutoff <- thresholds$metadata$completeness_warn_percent
  total_percent <- 100 * complete / max(total, 1)
  data.frame(
    metric = "metadata_completeness_percent",
    value = round(100 * complete / max(total, 1), 1),
    required_percent = round(100 * required_complete / max(nrow(required), 1), 1),
    optional_percent = round(100 * optional_complete / max(nrow(optional), 1), 1),
    required_complete = required_complete,
    required_total = nrow(required),
    optional_complete = optional_complete,
    optional_total = nrow(optional),
    pass_threshold_percent = pass_cutoff,
    warn_threshold_percent = warn_cutoff,
    status = ifelse(required_complete == nrow(required) && total_percent >= pass_cutoff, "PASS", ifelse(total_percent >= warn_cutoff, "WARN", "FAIL")),
    stringsAsFactors = FALSE
  )
}

validate_design <- function(plate_map, metadata, thresholds = default_qc_thresholds()) {
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
    flagged_metrics <- character()
    missing_controls_status <- "PASS"
    technical_replicates_status <- "PASS"
    dose_points_status <- "PASS"
    control_wells_status <- "PASS"
    inter_plate_calibrator_status <- "PASS"
    if (length(missing_controls)) {
      status <- "FAIL"
      missing_controls_status <- "FAIL"
      flagged_metrics <- c(flagged_metrics, "missing_controls")
      reasons <- c(reasons, paste("missing controls:", paste(missing_controls, collapse = ",")))
    }
    if (!is.na(min_tech) && min_tech < thresholds$design$technical_replicates_pass) {
      status <- ifelse(status == "FAIL", "FAIL", "WARN")
      technical_replicates_status <- "WARN"
      flagged_metrics <- c(flagged_metrics, "min_technical_replicates")
      reasons <- c(reasons, paste("technical replicates below", thresholds$design$technical_replicates_pass))
    }
    dose_response_design <- !is.na(min_doses) && min_doses > 1 && assay_mode %in% c("agonist", "antagonist")
    if (dose_response_design && min_doses < thresholds$design$dose_points_pass) {
      status <- ifelse(min_doses < thresholds$design$dose_points_warn, "FAIL", ifelse(status == "FAIL", "FAIL", "WARN"))
      dose_points_status <- ifelse(min_doses < thresholds$design$dose_points_warn, "FAIL", "WARN")
      flagged_metrics <- c(flagged_metrics, "min_dose_points")
      reasons <- c(reasons, paste("dose-response has fewer than", thresholds$design$dose_points_pass, "concentrations"))
    } else if (!is.na(min_doses) && min_doses == 1 && assay_mode %in% c("agonist", "antagonist")) {
      reasons <- c(reasons, "single-dose primary design detected; dose-response QC is not required")
    }
    if (min_control_wells < thresholds$design$control_wells_fail) {
      status <- "FAIL"
      control_wells_status <- "FAIL"
      flagged_metrics <- c(flagged_metrics, "min_control_wells")
      reasons <- c(reasons, "too few required control wells")
    } else if (min_control_wells < thresholds$design$control_wells_pass && status != "FAIL") {
      status <- "WARN"
      control_wells_status <- "WARN"
      flagged_metrics <- c(flagged_metrics, "min_control_wells")
      reasons <- c(reasons, paste("control wells below preferred count of", thresholds$design$control_wells_pass))
    }
    if (!ipc_present && assay_mode != "antagonist" && status != "FAIL") {
      status <- "WARN"
      inter_plate_calibrator_status <- "WARN"
      flagged_metrics <- c(flagged_metrics, "inter_plate_calibrator_present")
      reasons <- c(reasons, "inter-plate calibrator missing")
    }
    data.frame(
      plate_id = pid,
      assay_mode = assay_mode,
      design_status = status,
      flagged_metrics = ifelse(length(flagged_metrics), paste(unique(flagged_metrics), collapse = ";"), "0"),
      missing_controls = ifelse(length(missing_controls), paste(missing_controls, collapse = ";"), "0"),
      missing_controls_status = missing_controls_status,
      min_control_wells = min_control_wells,
      control_wells_status = control_wells_status,
      min_technical_replicates = min_tech,
      technical_replicates_status = technical_replicates_status,
      min_dose_points = min_doses,
      dose_points_status = dose_points_status,
      technical_replicate_threshold = thresholds$design$technical_replicates_pass,
      dose_points_pass_threshold = thresholds$design$dose_points_pass,
      control_wells_pass_threshold = thresholds$design$control_wells_pass,
      inter_plate_calibrator_present = ipc_present,
      inter_plate_calibrator_status = inter_plate_calibrator_status,
      notes = paste(reasons, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  design <- do.call(rbind, out)
  meta <- metadata_completeness(metadata, thresholds)
  design$metadata_status <- meta$status
  if (meta$status != "PASS") {
    design$flagged_metrics <- ifelse(design$flagged_metrics == "0", "metadata_completeness_percent", paste(design$flagged_metrics, "metadata_completeness_percent", sep = ";"))
  }
  design$metadata_completeness_percent <- meta$value
  design
}

calculate_plate_qc <- function(cleaned, thresholds = default_qc_thresholds()) {
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
    if (is.na(zp) || zp < thresholds$plate$z_prime_warn) {
      status <- "FAIL"
      notes <- c(notes, paste("Z-prime below", thresholds$plate$z_prime_warn))
    } else if (zp < thresholds$plate$z_prime_pass) {
      status <- "WARN"
      notes <- c(notes, paste("Z-prime below", thresholds$plate$z_prime_pass))
    }
    if (!is.na(cv_percent(pos)) && cv_percent(pos) > thresholds$plate$control_cv_warn_percent) {
      status <- ifelse(status == "FAIL", "FAIL", "WARN")
      notes <- c(notes, paste("positive control CV above", paste0(thresholds$plate$control_cv_warn_percent, "%")))
    }
    if (!is.na(cv_percent(neg)) && cv_percent(neg) > thresholds$plate$control_cv_warn_percent) {
      status <- ifelse(status == "FAIL", "FAIL", "WARN")
      notes <- c(notes, paste("negative control CV above", paste0(thresholds$plate$control_cv_warn_percent, "%")))
    }
    if (!is.na(edge) && edge > thresholds$plate$edge_effect_warn) {
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
      z_prime_pass_threshold = thresholds$plate$z_prime_pass,
      z_prime_warn_threshold = thresholds$plate$z_prime_warn,
      control_cv_warn_threshold = thresholds$plate$control_cv_warn_percent,
      edge_effect_warn_threshold = thresholds$plate$edge_effect_warn,
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

intraplate_variability_qc <- function(cleaned, thresholds = default_qc_thresholds()) {
  plates <- split(cleaned, cleaned$plate_id)
  out <- lapply(names(plates), function(pid) {
    x <- plates[[pid]]
    control <- x[x$control_type != "test_sample" & x$control_type != "empty", ]
    control_cv <- cv_percent(control$blank_corrected_od)
    edge <- edge_effect_score(x)
    row_bias <- row_bias_score(x)
    col_bias <- col_bias_score(x)
    outlier_rate <- 100 * mean(x$outlier_flag | x$saturated_flag | x$negative_corrected_flag | x$missing_flag, na.rm = TRUE)
    status <- "PASS"
    notes <- character()
    if (!is.na(control_cv) && control_cv > thresholds$plate$intraplate_cv_warn_percent) {
      status <- "WARN"
      notes <- c(notes, "within-plate control CV is high")
    }
    if (any(c(edge, row_bias, col_bias) > thresholds$plate$spatial_bias_warn, na.rm = TRUE)) {
      status <- "WARN"
      notes <- c(notes, "spatial bias exceeds threshold")
    }
    if (!is.na(outlier_rate) && outlier_rate > thresholds$plate$outlier_rate_warn_percent) {
      status <- "WARN"
      notes <- c(notes, "flagged-well rate exceeds threshold")
    }
    data.frame(
      plate_id = pid,
      assay_mode = unique(x$assay_mode)[1],
      intraplate_status = status,
      control_cv_percent = round(control_cv, 4),
      edge_effect = round(edge, 4),
      row_bias = round(row_bias, 4),
      column_bias = round(col_bias, 4),
      flagged_well_rate_percent = round(outlier_rate, 4),
      intraplate_cv_warn_threshold = thresholds$plate$intraplate_cv_warn_percent,
      spatial_bias_warn_threshold = thresholds$plate$spatial_bias_warn,
      outlier_rate_warn_threshold = thresholds$plate$outlier_rate_warn_percent,
      notes = paste(notes, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

reference_control_qc <- function(cleaned, thresholds = default_qc_thresholds()) {
  control_types <- c("negative_control", "positive_control", "agonist_challenge_control", "known_antagonist_control", "inter_plate_calibrator")
  x <- cleaned[cleaned$control_type %in% control_types, ]
  if (!nrow(x)) return(data.frame())
  agg <- aggregate(blank_corrected_od ~ plate_id + control_type, x, function(v) c(mean = safe_mean(v), median = safe_median(v), cv = cv_percent(v), n = sum(!is.na(v))))
  out <- cbind(agg[c("plate_id", "control_type")], as.data.frame(agg$blank_corrected_od))
  names(out)[3:6] <- c("mean_od", "median_od", "cv_percent", "n_wells")
  out$status <- ifelse(out$n_wells < thresholds$design$control_wells_fail, "FAIL", ifelse(out$cv_percent > thresholds$plate$control_cv_warn_percent, "WARN", "PASS"))
  out$min_wells_threshold <- thresholds$design$control_wells_fail
  out$cv_warn_threshold <- thresholds$plate$control_cv_warn_percent
  out$interpretation <- ifelse(out$status == "PASS", "Control is stable enough for normalization.",
    ifelse(out$status == "WARN", "Control variability is high. Review plate handling and pipetting.", "Too few control wells for reliable normalization.")
  )
  out
}

interplate_calibration <- function(cleaned, thresholds = default_qc_thresholds()) {
  calibrator <- cleaned[cleaned$control_type == "inter_plate_calibrator", ]
  method <- "inter_plate_calibrator"
  if (!nrow(calibrator)) {
    calibrator <- cleaned[cleaned$control_type %in% c("positive_control", "agonist_challenge_control"), ]
    method <- "shared_positive_control"
  }
  if (!nrow(calibrator)) {
    return(list(
      factors = data.frame(plate_id = unique(cleaned$plate_id), calibration_method = "none", calibration_factor = 1, calibration_drift_warn_threshold = thresholds$plate$calibration_drift_warn_od, calibration_status = "WARN", notes = "No shared calibrator or positive control found."),
      calibrated = transform(cleaned, calibrated_od = blank_corrected_od)
    ))
  }
  med <- aggregate(blank_corrected_od ~ plate_id, calibrator, safe_median)
  global <- safe_median(med$blank_corrected_od)
  med$calibration_method <- method
  med$calibration_factor <- global - med$blank_corrected_od
  med$calibrator_median_od <- med$blank_corrected_od
  med$calibration_drift_warn_threshold <- thresholds$plate$calibration_drift_warn_od
  med$calibration_status <- ifelse(abs(med$calibration_factor) <= thresholds$plate$calibration_drift_warn_od, "PASS", "WARN")
  med$notes <- ifelse(med$calibration_status == "PASS", "Plate is aligned with shared reference.", paste("Plate shows calibration drift above", thresholds$plate$calibration_drift_warn_od, "OD."))
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

summarize_primary <- function(normalized, thresholds = default_qc_thresholds()) {
  test <- normalized[normalized$control_type == "test_sample", ]
  if (!nrow(test)) return(data.frame())
  if (!"target_id" %in% names(test)) test$target_id <- NA_character_
  if (!"expected_activity" %in% names(test)) {
    test$expected_activity <- ifelse("compound_role" %in% names(test), test$compound_role, ifelse("cpd_role" %in% names(test), test$cpd_role, "unknown"))
  }
  test$expected_activity[is.na(test$expected_activity) | test$expected_activity == ""] <- "unknown"
  keys <- unique(test[c("plate_id", "assay_mode", "target_id", "peptide_id", "expected_activity", "concentration_uM")])
  result <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    x <- test[
      test$plate_id == key$plate_id &
        test$assay_mode == key$assay_mode &
        test$target_id == key$target_id &
        test$peptide_id == key$peptide_id &
        test$expected_activity == key$expected_activity &
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
  active_cv <- ifelse(result$assay_mode == "antagonist", result$inhibition_cv, result$activation_cv)
  active_mean <- ifelse(result$assay_mode == "antagonist", result$inhibition_mean, result$activation_mean)
  interpretable_cv <- is.finite(active_cv) & abs(active_mean) >= 10
  result$primary_hit <- with(result, ifelse(
    assay_mode == "agonist" & activation_mean >= thresholds$primary$agonist_hit_percent & activation_cv <= thresholds$primary$replicate_cv_warn_percent,
    "AGONIST_HIT",
    ifelse(assay_mode == "antagonist" & inhibition_mean >= thresholds$primary$antagonist_hit_percent & inhibition_cv <= thresholds$primary$replicate_cv_warn_percent, "ANTAGONIST_HIT", "NO_HIT")
  ))
  result$observed_direction <- with(result, ifelse(
    activation_mean >= thresholds$primary$agonist_hit_percent & inhibition_mean < thresholds$primary$antagonist_hit_percent,
    "agonist_like",
    ifelse(inhibition_mean >= thresholds$primary$antagonist_hit_percent & activation_mean < thresholds$primary$agonist_hit_percent,
      "antagonist_like",
      ifelse(activation_mean >= thresholds$primary$agonist_hit_percent & inhibition_mean >= thresholds$primary$antagonist_hit_percent, "mixed", "no_clear_effect")
    )
  ))
  expected <- tolower(result$expected_activity)
  result$direction_status <- ifelse(
    expected %in% c("unknown", "na", ""),
    "WARN",
    ifelse((expected %in% c("agonist", "agonist_like") & result$observed_direction == "agonist_like") |
      (expected %in% c("antagonist", "antagonist_like") & result$observed_direction == "antagonist_like"),
    "PASS", "WARN")
  )
  result$direction_note <- ifelse(
    result$direction_status == "PASS",
    "Observed response matches expected compound role.",
    ifelse(expected %in% c("unknown", "na", ""), "Expected compound role is unknown.", "Observed response does not match expected compound role.")
  )
  result$primary_status <- ifelse(
    interpretable_cv & active_cv > thresholds$primary$replicate_cv_fail_percent,
    "FAIL",
    ifelse(result$primary_hit != "NO_HIT" | (interpretable_cv & active_cv > thresholds$primary$replicate_cv_warn_percent), "WARN", "PASS")
  )
  result$primary_rule <- paste0(
    "hit >= ", ifelse(result$assay_mode == "antagonist", thresholds$primary$antagonist_hit_percent, thresholds$primary$agonist_hit_percent),
    "% and replicate CV <= ", thresholds$primary$replicate_cv_warn_percent,
    "%; FAIL if CV > ", thresholds$primary$replicate_cv_fail_percent, "%"
  )
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

fit_dose_response <- function(primary, thresholds = default_qc_thresholds()) {
  if (!nrow(primary)) return(list(results = data.frame(), qc = data.frame()))
  if (!"target_id" %in% names(primary)) primary$target_id <- NA_character_
  keys <- unique(primary[c("plate_id", "assay_mode", "target_id", "peptide_id")])
  rows <- list()
  qc <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    x <- primary[primary$plate_id == key$plate_id & primary$assay_mode == key$assay_mode & primary$target_id == key$target_id & primary$peptide_id == key$peptide_id, ]
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
      qc[[i]] <- data.frame(key, curve_qc_status = "FAIL", n_dose_points = n_dose_points, max_replicate_cv = round(max_rep_cv, 4), monotonic_violations = monotonic_violations, dynamic_range = round(dynamic_range, 4), ec50_in_range = FALSE, top_plateau_observed = FALSE, bottom_plateau_observed = FALSE, dose_points_pass_threshold = thresholds$curve$dose_points_pass, replicate_cv_warn_threshold = thresholds$curve$replicate_cv_warn_percent, dynamic_range_warn_threshold = thresholds$curve$dynamic_range_warn_percent, curve_flags = "FIT_FAILED")
      next
    }
    coef <- fit$coef
    dose_min <- min(x$concentration_uM, na.rm = TRUE)
    dose_max <- max(x$concentration_uM, na.rm = TRUE)
    flags <- character()
    ec50_in_range <- coef["ec50"] >= dose_min && coef["ec50"] <= dose_max
    top_plateau <- abs(response_max - coef["top"]) <= thresholds$curve$plateau_tolerance_percent
    bottom_plateau <- abs(response_min - coef["bottom"]) <= thresholds$curve$plateau_tolerance_percent
    if (n_dose_points < thresholds$curve$dose_points_pass) flags <- c(flags, "LOW_DOSE_COUNT")
    if (!ec50_in_range) flags <- c(flags, "EC50_OUT_OF_RANGE")
    if (abs(coef["top"] - coef["bottom"]) < thresholds$curve$dynamic_range_warn_percent) flags <- c(flags, "WEAK_RESPONSE")
    if (!top_plateau) flags <- c(flags, "NO_TOP_PLATEAU")
    if (!bottom_plateau) flags <- c(flags, "NO_BOTTOM_PLATEAU")
    if (abs(coef["hill"]) > thresholds$curve$hill_max || abs(coef["hill"]) < thresholds$curve$hill_min) flags <- c(flags, "BAD_HILL_SLOPE")
    if (monotonic_violations > thresholds$curve$monotonic_violations_warn) flags <- c(flags, "NON_MONOTONIC")
    high_noise <- !is.na(max_rep_cv) && max_rep_cv > thresholds$curve$replicate_cv_warn_percent
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
    qc[[i]] <- data.frame(key, curve_qc_status = status, n_dose_points = n_dose_points, max_replicate_cv = round(max_rep_cv, 4), monotonic_violations = monotonic_violations, dynamic_range = round(dynamic_range, 4), ec50_in_range = ec50_in_range, top_plateau_observed = top_plateau, bottom_plateau_observed = bottom_plateau, dose_points_pass_threshold = thresholds$curve$dose_points_pass, replicate_cv_warn_threshold = thresholds$curve$replicate_cv_warn_percent, dynamic_range_warn_threshold = thresholds$curve$dynamic_range_warn_percent, curve_flags = ifelse(length(flags), paste(flags, collapse = ";"), "GOOD_CURVE"))
  }
  list(results = do.call(rbind, rows), qc = do.call(rbind, qc))
}

counter_assay_qc <- function(normalized, thresholds = default_qc_thresholds()) {
  x <- normalized[normalized$assay_mode == "counter" & !is.na(normalized$peptide_id), ]
  if (!nrow(x)) return(data.frame())
  if (!"target_id" %in% names(x)) x$target_id <- NA_character_
  keys <- unique(x[c("target_id", "peptide_id")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    z <- x[x$target_id == key$target_id & x$peptide_id == key$peptide_id, ]
    viability <- safe_mean(z$raw_od[z$control_type == "viability_counter"])
    no_cell <- safe_mean(z$raw_od[z$control_type == "no_cell_interference"])
    unrelated <- safe_mean(z$raw_od[z$control_type == "unrelated_reporter"])
    null_cell <- safe_mean(z$raw_od[z$control_type == "null_cell_reporter"])
    flags <- character()
    if (!is.na(viability) && viability < thresholds$counter$viability_min_od) flags <- c(flags, "CYTOTOXICITY_CONFOUNDED")
    if (!is.na(no_cell) && no_cell > thresholds$counter$no_cell_max_od) flags <- c(flags, "ASSAY_INTERFERENCE")
    if (!is.na(unrelated) && unrelated > thresholds$counter$unrelated_reporter_max_od) flags <- c(flags, "UNRELATED_REPORTER_ACTIVITY")
    if (!is.na(null_cell) && null_cell > thresholds$counter$null_cell_max_od) flags <- c(flags, "NULL_CELL_ACTIVITY")
    data.frame(
      key,
      counter_qc_status = ifelse(length(flags), "WARN", "PASS"),
      viability_signal = round(viability, 3),
      no_cell_signal = round(no_cell, 3),
      unrelated_reporter_signal = round(unrelated, 3),
      null_cell_signal = round(null_cell, 3),
      viability_min_threshold = thresholds$counter$viability_min_od,
      no_cell_max_threshold = thresholds$counter$no_cell_max_od,
      unrelated_reporter_max_threshold = thresholds$counter$unrelated_reporter_max_od,
      null_cell_max_threshold = thresholds$counter$null_cell_max_od,
      artifact_flags = ifelse(length(flags), paste(flags, collapse = ";"), "NONE"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

sample_qc_summary <- function(primary, dose_qc, counter_qc) {
  if (!nrow(primary)) return(data.frame())
  groups <- unique(primary[c("target_id", "peptide_id", "assay_mode")])
  out <- lapply(seq_len(nrow(groups)), function(i) {
    key <- groups[i, ]
    p_primary <- primary[primary$target_id == key$target_id & primary$peptide_id == key$peptide_id & primary$assay_mode == key$assay_mode, , drop = FALSE]
    p_curve <- dose_qc[dose_qc$target_id == key$target_id & dose_qc$peptide_id == key$peptide_id & dose_qc$assay_mode == key$assay_mode, , drop = FALSE]
    p_counter <- counter_qc[counter_qc$target_id == key$target_id & counter_qc$peptide_id == key$peptide_id, , drop = FALSE]
    sample_status <- "PASS"
    reasons <- character()
    if (any(p_primary$primary_status == "FAIL", na.rm = TRUE)) {
      sample_status <- "FAIL"
      reasons <- c(reasons, "primary replicate QC failed")
    } else if (any(p_primary$primary_status == "WARN", na.rm = TRUE)) {
      sample_status <- "WARN"
      reasons <- c(reasons, "primary hit or replicate QC warning")
    }
    if (nrow(p_curve) && any(p_curve$curve_qc_status == "FAIL", na.rm = TRUE)) {
      sample_status <- "FAIL"
      reasons <- c(reasons, "curve QC failed")
    } else if (nrow(p_curve) && any(p_curve$curve_qc_status == "WARN", na.rm = TRUE) && sample_status != "FAIL") {
      sample_status <- "WARN"
      reasons <- c(reasons, "curve QC warning")
    }
    if (nrow(p_counter) && any(!p_counter$artifact_flags %in% c("NONE", NA), na.rm = TRUE)) {
      sample_status <- "FAIL"
      reasons <- c(reasons, "counter-assay artifact flag")
    }
    active_cv <- ifelse(p_primary$assay_mode == "antagonist", p_primary$inhibition_cv, p_primary$activation_cv)
    data.frame(
      key,
      sample_status = sample_status,
      n_doses = length(unique(p_primary$concentration_uM)),
      max_replicate_cv = round(max(active_cv, na.rm = TRUE), 4),
      hit_calls = paste(unique(p_primary$primary_hit[p_primary$primary_hit != "NO_HIT"]), collapse = ";"),
      curve_status = ifelse(nrow(p_curve), paste(unique(p_curve$curve_qc_status), collapse = ";"), "NOT_RUN"),
      artifact_flags = ifelse(nrow(p_counter), paste(unique(p_counter$artifact_flags), collapse = ";"), "NOT_TESTED"),
      review_notes = paste(reasons, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, out)
  result$hit_calls[result$hit_calls == ""] <- "NO_HIT"
  row.names(result) <- NULL
  result
}

make_final_qc <- function(design_qc, plate_qc, dose_qc, counter_qc, primary) {
  if (!nrow(primary)) {
    base <- data.frame(
      design_overall = ifelse(nrow(design_qc) && any(design_qc$design_status == "FAIL"), "FAIL", ifelse(nrow(design_qc) && any(design_qc$design_status == "WARN"), "WARN", "PASS")),
      plate_overall = ifelse(nrow(plate_qc) && any(plate_qc$plate_qc_status == "FAIL"), "FAIL", ifelse(nrow(plate_qc) && any(plate_qc$plate_qc_status == "WARN"), "WARN", "PASS")),
      stringsAsFactors = FALSE
    )
    if (nrow(counter_qc)) {
      out <- data.frame(
        final_status = ifelse(counter_qc$counter_qc_status == "PASS", "PASS", "WARN"),
        target_id = counter_qc$target_id,
        peptide_id = counter_qc$peptide_id,
        hit_call = "NO_PRIMARY_DATA",
        curve_qc = "NOT_RUN",
        counter_assay_qc = counter_qc$counter_qc_status,
        artifact_flags = counter_qc$artifact_flags,
        final_action = ifelse(counter_qc$artifact_flags == "NONE", "COUNTER_ONLY_REVIEW", "COUNTER_ONLY_ARTIFACT_REVIEW"),
        stringsAsFactors = FALSE
      )
      return(cbind(base[rep(1, nrow(out)), ], out))
    }
    return(cbind(base, data.frame(final_status = "WARN", target_id = NA_character_, peptide_id = NA_character_, hit_call = "NO_PRIMARY_DATA", curve_qc = "NOT_RUN", counter_assay_qc = "NOT_RUN", artifact_flags = "NOT_TESTED", final_action = "PARTIAL_DATA_REVIEW", stringsAsFactors = FALSE)))
  }
  peptides <- unique(primary[c("target_id", "peptide_id")])
  out <- lapply(seq_len(nrow(peptides)), function(i) {
    key <- peptides[i, ]
    p_primary <- primary[primary$target_id == key$target_id & primary$peptide_id == key$peptide_id, ]
    true_hits <- unique(p_primary$primary_hit[p_primary$primary_hit != "NO_HIT"])
    hit <- ifelse(length(true_hits), paste(true_hits, collapse = ";"), "NO_HIT")
    p_curve <- dose_qc[dose_qc$target_id == key$target_id & dose_qc$peptide_id == key$peptide_id, , drop = FALSE]
    p_counter <- counter_qc[counter_qc$target_id == key$target_id & counter_qc$peptide_id == key$peptide_id, , drop = FALSE]
    artifact <- if (nrow(p_counter)) p_counter$artifact_flags[1] else "NOT_TESTED"
    final <- "PASS_ADVANCE"
    if (hit == "NO_HIT") final <- "NO_ADVANCE"
    if (nrow(p_curve) && any(p_curve$curve_qc_status == "FAIL")) final <- "REVIEW_CURVE_NOISE"
    if (!artifact %in% c("NONE", "NOT_TESTED")) final <- "LIKELY_ARTIFACT"
    data.frame(
      final_status = ifelse(final %in% c("PASS_ADVANCE", "NO_ADVANCE"), "PASS", ifelse(final == "LIKELY_ARTIFACT", "FAIL", "WARN")),
      key,
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
  row.names(final) <- NULL
  cbind(plate_summary[rep(1, nrow(final)), ], final)
}

run_hekblue_analysis <- function(raw_data, plate_map, metadata, output_dir = NULL, thresholds = default_qc_thresholds(), threshold_change_note = "") {
  manifest <- assay_manifest(raw_data, plate_map, metadata, thresholds, threshold_change_note)
  raw_data <- enrich_raw_with_plate_map(raw_data, plate_map)
  raw_data <- fill_target_id(raw_data, metadata)
  plate_map <- fill_target_id(plate_map, metadata)
  cleaned <- clean_well_data(raw_data)
  eda <- input_eda(raw_data, plate_map)
  design <- validate_design(plate_map, metadata, thresholds)
  plate <- calculate_plate_qc(cleaned, thresholds)
  intraplate <- intraplate_variability_qc(cleaned, thresholds)
  ref_qc <- reference_control_qc(cleaned, thresholds)
  calibration <- interplate_calibration(cleaned, thresholds)
  normalized <- normalize_responses(calibration$calibrated)
  primary <- summarize_primary(normalized, thresholds)
  dose <- fit_dose_response(primary, thresholds)
  counter <- counter_assay_qc(normalized, thresholds)
  sample_qc <- sample_qc_summary(primary, dose$qc, counter)
  final <- make_final_qc(design, plate, dose$qc, counter, primary)
  exclusions <- cleaned[cleaned$saturated_flag | cleaned$negative_corrected_flag | cleaned$missing_flag, ]
  results <- list(
    qc_thresholds = qc_threshold_table(thresholds),
    assay_manifest = manifest,
    run_documentation = run_documentation_table(manifest, metadata, thresholds, threshold_change_note),
    metadata_completeness = metadata_completeness(metadata, thresholds),
    input_eda = eda,
    raw_data_summary = raw_data_summary(raw_data),
    metadata_summary = metadata_summary(metadata),
    cleaned_well_data = cleaned,
    interplate_calibration = calibration$factors,
    reference_control_qc = ref_qc,
    normalized_results = normalized,
    design_qc = design,
    plate_qc = plate,
    intraplate_variability_qc = intraplate,
    primary_results = primary,
    sample_qc_table = sample_qc,
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
