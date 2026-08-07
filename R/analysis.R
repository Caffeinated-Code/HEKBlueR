source_if_needed <- function(path) {
  if (file.exists(path)) source(path)
}

normalize_column_names <- function(df) {
  names(df) <- tolower(gsub("[^A-Za-z0-9]+", "_", trimws(names(df))))
  names(df) <- gsub("^_|_$", "", names(df))
  df
}

first_present <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else NA_character_
}

rename_first_present <- function(df, canonical, candidates) {
  if (canonical %in% names(df)) return(df)
  hit <- first_present(df, candidates)
  if (!is.na(hit)) names(df)[names(df) == hit] <- canonical
  df
}

normalize_well <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-H0-9]", "", x)
  row <- substr(x, 1, 1)
  col <- suppressWarnings(as.integer(gsub("^[A-H]", "", x)))
  ifelse(row %in% LETTERS[1:8] & !is.na(col) & col >= 1 & col <= 12, paste0(row, sprintf("%02d", col)), NA_character_)
}

standardize_control_type <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[^a-z0-9]+", "_", y)
  y <- gsub("^_|_$", "", y)
  synonyms <- c(
    vehicle = "negative_control",
    veh = "negative_control",
    mock = "negative_control",
    untreated = "negative_control",
    untreated_control = "negative_control",
    neg = "negative_control",
    negative = "negative_control",
    blank_control = "blank",
    media = "blank",
    medium = "blank",
    pos = "positive_control",
    positive = "positive_control",
    positive_agonist = "positive_control",
    agonist_challenge = "agonist_challenge_control",
    antagonist_control = "known_antagonist_control",
    known_antagonist = "known_antagonist_control",
    ipc = "inter_plate_calibrator",
    calibrator = "inter_plate_calibrator",
    test = "test_sample",
    treatment = "test_sample",
    compound = "test_sample",
    sample = "test_sample",
    viability = "viability_counter",
    cytotoxicity = "viability_counter",
    no_cell = "no_cell_interference",
    unrelated = "unrelated_reporter",
    null_cell = "null_cell_reporter",
    empty_well = "empty"
  )
  out <- unname(synonyms[y])
  out[is.na(out)] <- y[is.na(out)]
  out[is.na(out) | out == ""] <- "unknown"
  out
}

standardize_assay_mode <- function(x, control_type = NULL) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[^a-z0-9]+", "_", y)
  y[y %in% c("activation", "agonism", "agonist_screen")] <- "agonist"
  y[y %in% c("inhibition", "antagonism", "antagonist_screen")] <- "antagonist"
  y[y %in% c("counter_assay", "counterassay", "artifact", "viability")] <- "counter"
  y[is.na(y) | y == ""] <- "unknown"
  if (!is.null(control_type)) {
    y[control_type %in% c("viability_counter", "no_cell_interference", "unrelated_reporter", "null_cell_reporter")] <- "counter"
  }
  y
}

standardize_assay_stage <- function(x, assay_mode = NULL) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[^a-z0-9]+", "_", y)
  y[y %in% c("primary_screen", "screen", "fixed_dose")] <- "primary"
  y[y %in% c("secondary_confirmation", "confirmation", "dose_response", "curve")] <- "secondary"
  y[y %in% c("counter_assay", "counterassay", "artifact")] <- "counter"
  y[is.na(y) | y == ""] <- "unknown"
  if (!is.null(assay_mode)) y[assay_mode == "counter"] <- "counter"
  y
}

standardize_metadata <- function(metadata) {
  if (is.null(metadata) || !is.data.frame(metadata) || !nrow(metadata)) {
    return(data.frame(field = character(), value = character(), required = logical(), stringsAsFactors = FALSE))
  }
  metadata <- normalize_column_names(metadata)
  metadata <- rename_first_present(metadata, "field", c("metadata_field", "parameter", "key", "name"))
  metadata <- rename_first_present(metadata, "value", c("metadata_value", "val"))
  if (!"field" %in% names(metadata) || !"value" %in% names(metadata)) {
    if (nrow(metadata) == 1) {
      metadata <- data.frame(field = names(metadata), value = as.character(unlist(metadata[1, ], use.names = FALSE)), stringsAsFactors = FALSE)
    } else {
      stop("Metadata must contain field and value columns, or be a one-row wide metadata table.", call. = FALSE)
    }
  }
  if (!"required" %in% names(metadata)) metadata$required <- FALSE
  metadata$field <- tolower(gsub("[^A-Za-z0-9]+", "_", trimws(as.character(metadata$field))))
  metadata$field <- gsub("^_|_$", "", metadata$field)
  metadata$value <- trimws(as.character(metadata$value))
  metadata$required <- metadata$required %in% c(TRUE, "TRUE", "true", "1", "yes", "YES", "required", "Required")
  metadata[, c("field", "value", "required"), drop = FALSE]
}

standardize_well_table <- function(df, table_name = "raw data") {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) stop(table_name, " is empty.", call. = FALSE)
  df <- normalize_column_names(df)
  df <- rename_first_present(df, "plate_id", c("plate", "plate_name", "barcode", "plate_barcode"))
  df <- rename_first_present(df, "well", c("well_id", "well_position", "position"))
  df <- rename_first_present(df, "raw_od", c("od", "od655", "od_655", "absorbance", "absorbance_655", "value", "signal"))
  df <- rename_first_present(df, "sample_id", c("sample", "sample_name", "compound_id", "compound", "cpd_id", "treatment_id"))
  df <- rename_first_present(df, "peptide_id", c("peptide", "compound", "compound_id", "cpd_id", "treatment"))
  df <- rename_first_present(df, "control_type", c("well_role", "role", "sample_type", "type"))
  df <- rename_first_present(df, "concentration_uM", c("concentration", "conc", "dose", "dose_um", "concentration_um", "concentration_umol", "concentration_umol_l"))
  df <- rename_first_present(df, "assay_mode", c("mode", "direction", "assay_direction"))
  df <- rename_first_present(df, "assay_stage", c("stage", "assay_module", "module", "screen_stage"))
  missing_core <- setdiff(c("plate_id", "well"), names(df))
  if (length(missing_core)) stop(table_name, " is missing required column(s): ", paste(missing_core, collapse = ", "), call. = FALSE)
  df$plate_id <- trimws(as.character(df$plate_id))
  df$well <- normalize_well(df$well)
  if (any(is.na(df$well))) stop(table_name, " contains invalid well IDs. Use A01-H12 or A1-H12 style wells.", call. = FALSE)
  if (!"control_type" %in% names(df)) df$control_type <- "test_sample"
  df$control_type <- standardize_control_type(df$control_type)
  if (!"assay_mode" %in% names(df)) df$assay_mode <- "unknown"
  df$assay_mode <- standardize_assay_mode(df$assay_mode, df$control_type)
  if (!"assay_stage" %in% names(df)) df$assay_stage <- "unknown"
  df$assay_stage <- standardize_assay_stage(df$assay_stage, df$assay_mode)
  if (!"sample_id" %in% names(df)) df$sample_id <- ifelse(df$control_type == "test_sample", df$well, df$control_type)
  if (!"peptide_id" %in% names(df)) df$peptide_id <- NA_character_
  if (!"expected_activity" %in% names(df)) df$expected_activity <- NA_character_
  if (!"concentration_uM" %in% names(df)) df$concentration_uM <- NA_real_
  if (!"technical_replicate" %in% names(df)) df$technical_replicate <- NA_integer_
  if (!"biological_replicate" %in% names(df)) df$biological_replicate <- 1L
  df$row <- substr(df$well, 1, 1)
  df$col <- substr(df$well, 2, 3)
  df$concentration_uM <- suppressWarnings(as.numeric(df$concentration_uM))
  df$technical_replicate <- suppressWarnings(as.integer(df$technical_replicate))
  df$biological_replicate <- suppressWarnings(as.integer(df$biological_replicate))
  text_cols <- intersect(c("plate_id", "sample_id", "peptide_id", "expected_activity", "target_id"), names(df))
  for (col in text_cols) df[[col]] <- trimws(as.character(df[[col]]))
  if ("raw_od" %in% names(df)) df$raw_od <- suppressWarnings(as.numeric(df$raw_od))
  df <- df[!is.na(df$plate_id) & df$plate_id != "" & !is.na(df$well), , drop = FALSE]
  df
}

prepare_hekblue_inputs <- function(raw_data, plate_map = NULL, metadata = NULL) {
  metadata <- standardize_metadata(metadata)
  raw_data <- standardize_well_table(raw_data, "raw data")
  if (!"raw_od" %in% names(raw_data)) stop("raw data is missing required column: raw_od.", call. = FALSE)
  if (all(is.na(raw_data$raw_od))) stop("raw data has no numeric raw_od values.", call. = FALSE)
  if (any(duplicated(paste(raw_data$plate_id, raw_data$well)))) stop("raw data contains duplicate plate_id + well rows.", call. = FALSE)
  if (is.null(plate_map) || !is.data.frame(plate_map) || !nrow(plate_map)) {
    plate_map <- raw_data
  } else {
    plate_map <- standardize_well_table(plate_map, "plate map")
  }
  if (any(duplicated(paste(plate_map$plate_id, plate_map$well)))) stop("plate map contains duplicate plate_id + well rows.", call. = FALSE)
  missing_map <- setdiff(paste(raw_data$plate_id, raw_data$well), paste(plate_map$plate_id, plate_map$well))
  if (length(missing_map)) stop("plate map is missing ", length(missing_map), " raw data well(s).", call. = FALSE)
  list(raw_data = raw_data, plate_map = plate_map, metadata = metadata)
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
  blank_rows <- df[df$control_type == "blank", , drop = FALSE]
  if (nrow(blank_rows)) {
    blanks <- aggregate(raw_od ~ plate_id, blank_rows, safe_median)
    names(blanks)[2] <- "blank_od"
  } else {
    blanks <- data.frame(plate_id = unique(df$plate_id), blank_od = 0, stringsAsFactors = FALSE)
  }
  df <- merge(df, blanks, by = "plate_id", all.x = TRUE)
  df$blank_od[is.na(df$blank_od)] <- 0
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
  raw_od_values <- if ("raw_od" %in% names(raw_data)) raw_data$raw_od else numeric()
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
      round(safe_min(raw_od_values), 4),
      round(safe_median(raw_od_values), 4),
      round(safe_max(raw_od_values), 4)
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
  if (is.null(plate_map) || !nrow(plate_map)) return(data.frame())
  controls <- split(plate_map, plate_map$plate_id)
  out <- lapply(names(controls), function(pid) {
    x <- controls[[pid]]
    assay_mode <- unique(x$assay_mode)[1]
    assay_stage <- if ("assay_stage" %in% names(x)) unique(x$assay_stage)[1] else "unknown"
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
    dose_response_design <- assay_stage %in% c("secondary", "unknown") && !is.na(min_doses) && min_doses > 1 && assay_mode %in% c("agonist", "antagonist")
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
    if (!ipc_present && assay_mode != "antagonist" && assay_stage != "counter" && status != "FAIL") {
      status <- "WARN"
      inter_plate_calibrator_status <- "WARN"
      flagged_metrics <- c(flagged_metrics, "inter_plate_calibrator_present")
      reasons <- c(reasons, "inter-plate calibrator missing")
    }
    data.frame(
      plate_id = pid,
      assay_stage = assay_stage,
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
  if (is.null(cleaned) || !nrow(cleaned)) return(data.frame())
  plates <- split(cleaned, cleaned$plate_id)
  out <- lapply(names(plates), function(pid) {
    x <- plates[[pid]]
    pos <- x$blank_corrected_od[x$control_type %in% c("positive_control", "agonist_challenge_control")]
    neg <- x$blank_corrected_od[x$control_type == "negative_control"]
    blank <- x$raw_od[x$control_type == "blank"]
    zp <- z_prime(pos, neg)
    rzp <- robust_z_prime(pos, neg)
    neg_mean <- safe_mean(neg)
    sb <- ifelse(is.na(neg_mean) || abs(neg_mean) < .Machine$double.eps, NA_real_, safe_mean(pos) / abs(neg_mean))
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
  if (is.null(cleaned) || !nrow(cleaned)) return(data.frame())
  plates <- split(cleaned, cleaned$plate_id)
  out <- lapply(names(plates), function(pid) {
    x <- plates[[pid]]
    control <- x[x$control_type != "test_sample" & x$control_type != "empty", ]
    control_cv <- cv_percent(control$blank_corrected_od)
    edge <- edge_effect_score(x)
    row_bias <- row_bias_score(x)
    col_bias <- col_bias_score(x)
    flagged <- x$outlier_flag | x$saturated_flag | x$negative_corrected_flag | x$missing_flag
    outlier_rate <- if (length(flagged)) 100 * mean(flagged, na.rm = TRUE) else NA_real_
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

score_component <- function(value, good, warn, direction = c("lte", "gte"), missing_grade = "FAIL") {
  direction <- match.arg(direction)
  if (is.na(value)) {
    return(list(points = ifelse(missing_grade == "WARN", 0.5, 0), grade = missing_grade))
  }
  if (direction == "lte") {
    if (value <= good) return(list(points = 1, grade = "PASS"))
    if (value <= warn) return(list(points = 0.5, grade = "WARN"))
    return(list(points = 0, grade = "FAIL"))
  }
  if (value >= good) return(list(points = 1, grade = "PASS"))
  if (value >= warn) return(list(points = 0.5, grade = "WARN"))
  list(points = 0, grade = "FAIL")
}

reference_stability_row <- function(n_wells, cv_percent_value, calibration_drift_abs, spatial_bias, outlier_rate, thresholds = default_qc_thresholds()) {
  weights <- c(well_count = 25, control_cv = 25, calibration_drift = 25, spatial_bias = 15, outlier_rate = 10)
  well <- score_component(n_wells, thresholds$design$control_wells_pass, thresholds$design$control_wells_fail, direction = "gte")
  cv <- score_component(cv_percent_value, thresholds$plate$control_cv_warn_percent, thresholds$plate$control_cv_fail_percent, direction = "lte")
  drift <- score_component(calibration_drift_abs, thresholds$plate$calibration_drift_warn_od, thresholds$plate$calibration_drift_fail_od, direction = "lte")
  spatial <- score_component(spatial_bias, thresholds$plate$spatial_bias_warn, thresholds$plate$spatial_bias_fail, direction = "lte", missing_grade = "WARN")
  outlier <- score_component(outlier_rate, thresholds$plate$outlier_rate_warn_percent, thresholds$plate$outlier_rate_fail_percent, direction = "lte", missing_grade = "WARN")
  components <- list(well_count = well, control_cv = cv, calibration_drift = drift, spatial_bias = spatial, outlier_rate = outlier)
  score <- sum(vapply(names(components), function(name) weights[[name]] * components[[name]]$points, numeric(1)))
  hard_fail <- well$grade == "FAIL" || cv$grade == "FAIL" || drift$grade == "FAIL"
  status <- ifelse(
    hard_fail || score < thresholds$plate$reference_stability_warn_score,
    "FAIL",
    ifelse(score >= thresholds$plate$reference_stability_pass_score && !any(vapply(components, function(x) x$grade == "WARN", logical(1))), "PASS", "WARN")
  )
  flags <- names(components)[vapply(components, function(x) x$grade != "PASS", logical(1))]
  list(
    score = round(score, 1),
    status = status,
    flags = ifelse(length(flags), paste(flags, collapse = ";"), "0"),
    well_count_grade = well$grade,
    control_cv_grade = cv$grade,
    calibration_drift_grade = drift$grade,
    spatial_bias_grade = spatial$grade,
    outlier_rate_grade = outlier$grade
  )
}

reference_control_qc <- function(cleaned, thresholds = default_qc_thresholds()) {
  control_types <- c("negative_control", "positive_control", "agonist_challenge_control", "known_antagonist_control", "inter_plate_calibrator")
  x <- cleaned[cleaned$control_type %in% control_types, ]
  if (!nrow(x)) return(data.frame())
  agg <- aggregate(blank_corrected_od ~ plate_id + control_type, x, function(v) c(mean = safe_mean(v), median = safe_median(v), cv = cv_percent(v), n = sum(!is.na(v))))
  out <- cbind(agg[c("plate_id", "control_type")], as.data.frame(agg$blank_corrected_od))
  names(out)[3:6] <- c("mean_od", "median_od", "cv_percent", "n_wells")
  out$global_control_median_od <- ave(out$median_od, out$control_type, FUN = safe_median)
  out$calibration_drift_od <- out$global_control_median_od - out$median_od
  out$absolute_calibration_drift_od <- abs(out$calibration_drift_od)
  plate_spatial <- intraplate_variability_qc(cleaned, thresholds)
  out$spatial_bias <- NA_real_
  out$outlier_rate_percent <- NA_real_
  if (nrow(plate_spatial)) {
    out$spatial_bias <- pmax(
      plate_spatial$edge_effect[match(out$plate_id, plate_spatial$plate_id)],
      plate_spatial$row_bias[match(out$plate_id, plate_spatial$plate_id)],
      plate_spatial$column_bias[match(out$plate_id, plate_spatial$plate_id)],
      na.rm = TRUE
    )
    out$spatial_bias[!is.finite(out$spatial_bias)] <- NA_real_
    out$outlier_rate_percent <- plate_spatial$flagged_well_rate_percent[match(out$plate_id, plate_spatial$plate_id)]
  }
  stability <- lapply(seq_len(nrow(out)), function(i) {
    reference_stability_row(out$n_wells[i], out$cv_percent[i], out$absolute_calibration_drift_od[i], out$spatial_bias[i], out$outlier_rate_percent[i], thresholds)
  })
  out$reference_stability_score <- vapply(stability, `[[`, numeric(1), "score")
  out$reference_stability_status <- vapply(stability, `[[`, character(1), "status")
  out$stability_flags <- vapply(stability, `[[`, character(1), "flags")
  out$well_count_grade <- vapply(stability, `[[`, character(1), "well_count_grade")
  out$control_cv_grade <- vapply(stability, `[[`, character(1), "control_cv_grade")
  out$calibration_drift_grade <- vapply(stability, `[[`, character(1), "calibration_drift_grade")
  out$spatial_bias_grade <- vapply(stability, `[[`, character(1), "spatial_bias_grade")
  out$outlier_rate_grade <- vapply(stability, `[[`, character(1), "outlier_rate_grade")
  out$status <- out$reference_stability_status
  out$min_wells_threshold <- thresholds$design$control_wells_fail
  out$cv_warn_threshold <- thresholds$plate$control_cv_warn_percent
  out$cv_fail_threshold <- thresholds$plate$control_cv_fail_percent
  out$calibration_drift_warn_threshold <- thresholds$plate$calibration_drift_warn_od
  out$calibration_drift_fail_threshold <- thresholds$plate$calibration_drift_fail_od
  out$reference_stability_pass_threshold <- thresholds$plate$reference_stability_pass_score
  out$reference_stability_warn_threshold <- thresholds$plate$reference_stability_warn_score
  out$interpretation <- ifelse(out$status == "PASS", "Control is stable enough for unflagged normalization.",
    ifelse(out$status == "WARN", "Control is used, but normalized values should carry a stability warning.", "Control is dropped from calibration, percent-response, and fold-change calculations.")
  )
  out
}

interplate_calibration <- function(cleaned, reference_qc = reference_control_qc(cleaned, thresholds), thresholds = default_qc_thresholds()) {
  all_plates <- data.frame(plate_id = unique(cleaned$plate_id), stringsAsFactors = FALSE)
  if (is.null(reference_qc) || !nrow(reference_qc)) {
    factors <- transform(all_plates, calibration_method = "none", control_type = NA_character_, calibrator_median_od = NA_real_, calibration_factor = 0, reference_stability_score = NA_real_, reference_stability_status = "FAIL", calibration_drift_warn_threshold = thresholds$plate$calibration_drift_warn_od, calibration_status = "WARN", notes = "No reference control QC was available.")
    return(list(factors = factors, calibrated = transform(cleaned, calibrated_od = blank_corrected_od)))
  }

  eligible_qc <- reference_qc[reference_qc$reference_stability_status != "FAIL", , drop = FALSE]
  priority <- c("inter_plate_calibrator", "positive_control", "agonist_challenge_control")
  chosen_type <- NA_character_
  for (candidate in priority) {
    candidate_plates <- unique(eligible_qc$plate_id[eligible_qc$control_type == candidate])
    if (length(candidate_plates) >= min(2, length(unique(cleaned$plate_id)))) {
      chosen_type <- candidate
      break
    }
  }

  if (is.na(chosen_type)) {
    factors <- transform(all_plates, calibration_method = "none", control_type = NA_character_, calibrator_median_od = NA_real_, calibration_factor = 0, reference_stability_score = NA_real_, reference_stability_status = "FAIL", calibration_drift_warn_threshold = thresholds$plate$calibration_drift_warn_od, calibration_status = "WARN", notes = "No eligible shared reference control passed stability screening.")
    return(list(factors = factors, calibrated = transform(cleaned, calibrated_od = blank_corrected_od)))
  }

  stable_keys <- paste(eligible_qc$plate_id, eligible_qc$control_type)
  calibrator <- cleaned[cleaned$control_type == chosen_type & paste(cleaned$plate_id, cleaned$control_type) %in% stable_keys, , drop = FALSE]
  med <- aggregate(blank_corrected_od ~ plate_id + control_type, calibrator, safe_median)
  global <- safe_median(med$blank_corrected_od)
  med$calibration_method <- ifelse(chosen_type == "inter_plate_calibrator", "inter_plate_calibrator", "shared_positive_control")
  med$calibration_factor <- global - med$blank_corrected_od
  med$calibrator_median_od <- med$blank_corrected_od
  key <- paste(med$plate_id, med$control_type)
  qc_key <- paste(eligible_qc$plate_id, eligible_qc$control_type)
  med$reference_stability_score <- eligible_qc$reference_stability_score[match(key, qc_key)]
  med$reference_stability_status <- eligible_qc$reference_stability_status[match(key, qc_key)]
  med$calibration_drift_warn_threshold <- thresholds$plate$calibration_drift_warn_od
  med$calibration_status <- ifelse(abs(med$calibration_factor) <= thresholds$plate$calibration_drift_warn_od & med$reference_stability_status == "PASS", "PASS", "WARN")
  med$notes <- ifelse(med$calibration_status == "PASS", "Plate is aligned with a stable shared reference.", "Calibration reference is used with warning; review drift and stability score.")
  med$blank_corrected_od <- NULL

  factors <- merge(all_plates, med, by = "plate_id", all.x = TRUE, sort = FALSE)
  factors$calibration_method[is.na(factors$calibration_method)] <- "none_dropped_or_missing"
  factors$calibration_factor[is.na(factors$calibration_factor)] <- 0
  factors$reference_stability_status[is.na(factors$reference_stability_status)] <- "FAIL"
  factors$calibration_status[is.na(factors$calibration_status)] <- "WARN"
  factors$calibration_drift_warn_threshold[is.na(factors$calibration_drift_warn_threshold)] <- thresholds$plate$calibration_drift_warn_od
  factors$notes[is.na(factors$notes)] <- "The candidate calibration control failed stability screening or was missing on this plate."

  calibrated <- merge(cleaned, factors[c("plate_id", "calibration_factor")], by = "plate_id", all.x = TRUE)
  calibrated$calibration_factor[is.na(calibrated$calibration_factor)] <- 0
  calibrated$calibrated_od <- calibrated$blank_corrected_od + calibrated$calibration_factor
  list(factors = factors, calibrated = calibrated)
}

select_stable_control <- function(x, reference_qc, control_types, value_col) {
  plate <- unique(x$plate_id)[1]
  empty <- list(value = NA_real_, control_type = NA_character_, status = "FAIL", score = NA_real_, flags = "missing_control")
  if (is.null(reference_qc) || !nrow(reference_qc)) return(empty)
  rq <- reference_qc[reference_qc$plate_id == plate & reference_qc$control_type %in% control_types, , drop = FALSE]
  rq <- rq[rq$reference_stability_status != "FAIL", , drop = FALSE]
  if (!nrow(rq)) return(empty)
  rq$priority <- match(rq$control_type, control_types)
  rq <- rq[order(rq$priority, -rq$reference_stability_score), , drop = FALSE]
  selected <- rq[1, ]
  value <- safe_median(x[[value_col]][x$control_type == selected$control_type])
  if (is.na(value)) return(empty)
  list(
    value = value,
    control_type = selected$control_type,
    status = selected$reference_stability_status,
    score = selected$reference_stability_score,
    flags = selected$stability_flags
  )
}

safe_ratio <- function(numerator, denominator) {
  ifelse(is.na(denominator) | abs(denominator) < .Machine$double.eps | denominator <= 0 | is.na(numerator) | numerator <= 0, NA_real_, numerator / denominator)
}

normalization_status <- function(...) {
  statuses <- c(...)
  if (any(statuses == "FAIL")) return("FAIL")
  if (any(statuses == "WARN")) return("WARN")
  "PASS"
}

normalization_flags <- function(named_statuses) {
  bad <- names(named_statuses)[named_statuses == "FAIL"]
  warn <- names(named_statuses)[named_statuses == "WARN"]
  flags <- c(
    if (length(bad)) paste0("dropped_failed_or_missing_", bad),
    if (length(warn)) paste0("warning_", warn)
  )
  ifelse(length(flags), paste(flags, collapse = ";"), "0")
}

normalize_responses <- function(cleaned, reference_qc = NULL) {
  if (is.null(cleaned) || !nrow(cleaned)) return(data.frame())
  plates <- split(cleaned, cleaned$plate_id)
  out <- lapply(plates, function(x) {
    value_col <- if ("calibrated_od" %in% names(x)) "calibrated_od" else "blank_corrected_od"
    neg <- select_stable_control(x, reference_qc, "negative_control", value_col)
    activation_ref <- select_stable_control(x, reference_qc, c("positive_control", "agonist_challenge_control"), value_col)
    agonist <- select_stable_control(x, reference_qc, "agonist_challenge_control", value_col)
    antagonist <- select_stable_control(x, reference_qc, "known_antagonist_control", value_col)
    denom_activation <- activation_ref$value - neg$value
    denom_inhibition <- agonist$value - antagonist$value
    x$percent_activation <- if (is.na(denom_activation) || abs(denom_activation) < .Machine$double.eps) NA_real_ else 100 * (x[[value_col]] - neg$value) / denom_activation
    x$percent_inhibition <- if (is.na(denom_inhibition) || abs(denom_inhibition) < .Machine$double.eps) NA_real_ else 100 * (agonist$value - x[[value_col]]) / denom_inhibition
    x$fold_change_vs_negative <- safe_ratio(x[[value_col]], neg$value)
    x$log2_fold_change_vs_negative <- log2(x$fold_change_vs_negative)
    x$fold_change_vs_activation_reference <- safe_ratio(x[[value_col]], activation_ref$value)
    x$log2_fold_change_vs_activation_reference <- log2(x$fold_change_vs_activation_reference)
    x$fold_change_vs_agonist_challenge <- safe_ratio(x[[value_col]], agonist$value)
    x$log2_fold_change_vs_agonist_challenge <- log2(x$fold_change_vs_agonist_challenge)
    x$negative_control_type <- neg$control_type
    x$negative_control_stability_status <- neg$status
    x$negative_control_stability_score <- neg$score
    x$activation_reference_control_type <- activation_ref$control_type
    x$activation_reference_stability_status <- activation_ref$status
    x$activation_reference_stability_score <- activation_ref$score
    x$agonist_challenge_control_stability_status <- agonist$status
    x$known_antagonist_control_stability_status <- antagonist$status
    x$activation_normalization_status <- normalization_status(negative_control = neg$status, activation_reference = activation_ref$status)
    x$inhibition_normalization_status <- normalization_status(agonist_challenge_control = agonist$status, known_antagonist_control = antagonist$status)
    x$fold_change_negative_status <- normalization_status(negative_control = neg$status)
    x$fold_change_activation_reference_status <- normalization_status(activation_reference = activation_ref$status)
    x$fold_change_agonist_challenge_status <- normalization_status(agonist_challenge_control = agonist$status)
    x$activation_normalization_flags <- normalization_flags(c(negative_control = neg$status, activation_reference = activation_ref$status))
    x$inhibition_normalization_flags <- normalization_flags(c(agonist_challenge_control = agonist$status, known_antagonist_control = antagonist$status))
    x$fold_change_flags <- normalization_flags(c(negative_control = neg$status, activation_reference = activation_ref$status, agonist_challenge_control = agonist$status))
    x$normalization_control_status <- ifelse(x$assay_mode == "antagonist", x$inhibition_normalization_status, ifelse(x$assay_mode == "agonist", x$activation_normalization_status, normalization_status(negative_control = neg$status, activation_reference = activation_ref$status)))
    x$normalization_control_flags <- ifelse(x$assay_mode == "antagonist", x$inhibition_normalization_flags, ifelse(x$assay_mode == "agonist", x$activation_normalization_flags, x$fold_change_flags))
    x
  })
  do.call(rbind, out)
}

summarize_primary <- function(normalized, thresholds = default_qc_thresholds()) {
  test <- normalized[normalized$control_type == "test_sample", ]
  if (!nrow(test)) return(data.frame())
  if (!"target_id" %in% names(test)) test$target_id <- NA_character_
  if (!"assay_stage" %in% names(test)) test$assay_stage <- "unknown"
  if (!"expected_activity" %in% names(test)) {
    test$expected_activity <- ifelse("compound_role" %in% names(test), test$compound_role, ifelse("cpd_role" %in% names(test), test$cpd_role, "unknown"))
  }
  test$expected_activity[is.na(test$expected_activity) | test$expected_activity == ""] <- "unknown"
  keys <- unique(test[c("plate_id", "assay_stage", "assay_mode", "target_id", "peptide_id", "expected_activity", "concentration_uM")])
  result <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    x <- test[
        test$plate_id == key$plate_id &
        test$assay_stage == key$assay_stage &
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
  result$primary_hit[is.na(result$primary_hit)] <- "NO_HIT"
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
  fit <- suppressWarnings(try(stats::nls(
    y ~ bottom + (top - bottom) / (1 + (ec50 / x)^hill),
    start = start,
    control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
  ), silent = TRUE))
  if (inherits(fit, "try-error")) return(NULL)
  pred <- stats::predict(fit)
  list(coef = coef(fit), rmse = sqrt(mean((y - pred)^2, na.rm = TRUE)), residual_max = max(abs(y - pred), na.rm = TRUE))
}

fit_dose_response <- function(primary, thresholds = default_qc_thresholds()) {
  if (!nrow(primary)) return(list(results = data.frame(), qc = data.frame()))
  if (!"target_id" %in% names(primary)) primary$target_id <- NA_character_
  if (!"assay_stage" %in% names(primary)) primary$assay_stage <- "unknown"
  primary <- primary[primary$assay_stage %in% c("secondary", "unknown"), , drop = FALSE]
  if (!nrow(primary)) return(list(results = data.frame(), qc = data.frame()))
  keys <- unique(primary[c("plate_id", "assay_stage", "assay_mode", "target_id", "peptide_id")])
  rows <- list()
  qc <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    x <- primary[primary$plate_id == key$plate_id & primary$assay_stage == key$assay_stage & primary$assay_mode == key$assay_mode & primary$target_id == key$target_id & primary$peptide_id == key$peptide_id, ]
    response_col <- if (key$assay_mode == "antagonist") "inhibition_mean" else "activation_mean"
    fit <- fit_one_curve(x, response_col)
    response <- x[[response_col]]
    n_dose_points <- length(unique(x$concentration_uM[is.finite(x$concentration_uM)]))
    response_min <- safe_min(response)
    response_max <- safe_max(response)
    dynamic_range <- response_max - response_min
    ordered_x <- x[order(x$concentration_uM), ]
    ordered_response <- ordered_x[[response_col]]
    monotonic_violations <- if (key$assay_mode == "antagonist") {
      sum(diff(ordered_response) < -10, na.rm = TRUE)
    } else {
      sum(diff(ordered_response) < -10, na.rm = TRUE)
    }
    response_for_cv <- x[[response_col]]
    cv_col <- paste0(ifelse(key$assay_mode == "antagonist", "inhibition", "activation"), "_cv")
    interpretable_points <- is.finite(response_for_cv) & abs(response_for_cv) >= 10
    max_rep_cv <- max(x[[cv_col]][interpretable_points], na.rm = TRUE)
    if (!is.finite(max_rep_cv)) max_rep_cv <- NA_real_
    if (is.null(fit)) {
      rows[[i]] <- data.frame(key, n_dose_points = n_dose_points, bottom = NA, top = NA, ec50_ic50 = NA, hill = NA, response_min = round(response_min, 4), response_max = round(response_max, 4), dynamic_range = round(dynamic_range, 4), rmse = NA, max_residual = NA, curve_status = "FIT_FAILED")
      qc[[i]] <- data.frame(key, curve_qc_status = "FAIL", n_dose_points = n_dose_points, max_replicate_cv = round(max_rep_cv, 4), monotonic_violations = monotonic_violations, dynamic_range = round(dynamic_range, 4), ec50_in_range = FALSE, top_plateau_observed = FALSE, bottom_plateau_observed = FALSE, dose_points_pass_threshold = thresholds$curve$dose_points_pass, replicate_cv_warn_threshold = thresholds$curve$replicate_cv_warn_percent, dynamic_range_warn_threshold = thresholds$curve$dynamic_range_warn_percent, curve_flags = "FIT_FAILED")
      next
    }
    coef <- fit$coef
    dose_min <- safe_min(x$concentration_uM)
    dose_max <- safe_max(x$concentration_uM)
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
  if (!nrow(dose_qc)) {
    dose_qc <- data.frame(target_id = character(), peptide_id = character(), assay_mode = character(), curve_qc_status = character(), stringsAsFactors = FALSE)
  }
  if (!nrow(counter_qc)) {
    counter_qc <- data.frame(target_id = character(), peptide_id = character(), counter_qc_status = character(), artifact_flags = character(), stringsAsFactors = FALSE)
  }
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
    max_cv <- safe_max(active_cv)
    data.frame(
      key,
      sample_status = sample_status,
      n_doses = length(unique(p_primary$concentration_uM)),
      max_replicate_cv = round(max_cv, 4),
      hit_calls = paste(unique(p_primary$primary_hit[!is.na(p_primary$primary_hit) & p_primary$primary_hit != "NO_HIT"]), collapse = ";"),
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
  if (!nrow(dose_qc)) {
    dose_qc <- data.frame(target_id = character(), peptide_id = character(), curve_qc_status = character(), stringsAsFactors = FALSE)
  }
  if (!nrow(counter_qc)) {
    counter_qc <- data.frame(target_id = character(), peptide_id = character(), counter_qc_status = character(), artifact_flags = character(), stringsAsFactors = FALSE)
  }
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
    true_hits <- unique(p_primary$primary_hit[!is.na(p_primary$primary_hit) & p_primary$primary_hit != "NO_HIT"])
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
  prepared <- prepare_hekblue_inputs(raw_data, plate_map, metadata)
  raw_data <- prepared$raw_data
  plate_map <- prepared$plate_map
  metadata <- prepared$metadata
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
  calibration <- interplate_calibration(cleaned, ref_qc, thresholds)
  normalized <- normalize_responses(calibration$calibrated, ref_qc)
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
    hit_calls = primary[!is.na(primary$primary_hit) & primary$primary_hit != "NO_HIT", ],
    exclusions = exclusions,
    final_qc_table = final
  )
  if (!is.null(output_dir)) export_hekblue_results(results, raw_data, plate_map, metadata, output_dir)
  results
}
