make_wells <- function(format = 96) {
  if (format != 96) stop("The demo generator currently supports 96-well plates.")
  rows <- LETTERS[1:8]
  cols <- sprintf("%02d", 1:12)
  expand.grid(row = rows, col = cols, stringsAsFactors = FALSE) |>
    transform(well = paste0(row, col))
}

clip_od <- function(x) pmax(0.02, pmin(2.8, x))

response_4pl <- function(conc, bottom, top, ec50, hill = 1) {
  bottom + (top - bottom) / (1 + (ec50 / conc)^hill)
}

simulate_od <- function(row, plate_shift = 0, noise_sd = 0.035) {
  ct <- row$control_type
  pep <- row$peptide_id
  conc <- row$concentration_uM
  base <- switch(
    ct,
    blank = 0.075,
    negative_control = 0.18,
    positive_control = 1.32,
    agonist_challenge_control = 1.18,
    known_antagonist_control = 0.34,
    inter_plate_calibrator = 0.76,
    empty = NA_real_,
    0.20
  )
  if (ct == "test_sample" && row$assay_stage == "primary") {
    base <- c("PEP-001" = 1.08, "PEP-002" = 0.24, "PEP-003" = 0.82, "PEP-004" = 0.21, "PEP-005" = 0.52, "PEP-006" = 0.33)[pep]
  }
  if (ct == "test_sample" && row$assay_stage == "secondary" && row$assay_mode == "agonist") {
    top <- c("PEP-001" = 1.20, "PEP-002" = 0.42, "PEP-003" = 0.82, "PEP-004" = 0.25, "PEP-005" = 0.68, "PEP-006" = 0.40)[pep]
    ec50 <- c("PEP-001" = 0.22, "PEP-002" = 1.8, "PEP-003" = 0.75, "PEP-004" = 5.5, "PEP-005" = 1.1, "PEP-006" = 2.8)[pep]
    base <- response_4pl(conc, 0.18, top, ec50, 1.05)
  }
  if (ct == "test_sample" && row$assay_stage == "secondary" && row$assay_mode == "antagonist") {
    max_inhib <- c("PEP-001" = 0.12, "PEP-002" = 0.84, "PEP-003" = 0.44, "PEP-004" = 0.10, "PEP-005" = 0.28, "PEP-006" = 0.18)[pep]
    ic50 <- c("PEP-001" = 4.0, "PEP-002" = 0.18, "PEP-003" = 0.70, "PEP-004" = 7.0, "PEP-005" = 1.4, "PEP-006" = 2.6)[pep]
    base <- 1.18 - max_inhib * conc / (ic50 + conc)
  }
  if (ct == "viability_counter") {
    base <- c("PEP-001" = 0.96, "PEP-002" = 0.91, "PEP-003" = 0.50, "PEP-004" = 0.98, "PEP-005" = 0.78, "PEP-006" = 0.88)[pep]
  }
  if (ct == "no_cell_interference") {
    base <- c("PEP-001" = 0.08, "PEP-002" = 0.09, "PEP-003" = 0.40, "PEP-004" = 0.08, "PEP-005" = 0.15, "PEP-006" = 0.10)[pep]
  }
  if (ct == "unrelated_reporter") {
    base <- c("PEP-001" = 0.18, "PEP-002" = 0.16, "PEP-003" = 0.64, "PEP-004" = 0.17, "PEP-005" = 0.31, "PEP-006" = 0.18)[pep]
  }
  if (ct == "null_cell_reporter") {
    base <- c("PEP-001" = 0.14, "PEP-002" = 0.15, "PEP-003" = 0.48, "PEP-004" = 0.13, "PEP-005" = 0.24, "PEP-006" = 0.15)[pep]
  }
  if (is.na(base)) return(NA_real_)
  edge <- ifelse(row$row %in% c("A", "H"), 0.020, 0) + ifelse(row$col %in% c("01", "12"), 0.015, 0)
  clip_od(round(base + plate_shift + edge + rnorm(1, 0, noise_sd), 4))
}

simulate_hekblue_demo <- function(seed = 42) {
  set.seed(seed)
  target <- data.frame(
    target_id = "TARGET_TLR8_DEMO",
    target_name = "TLR8 reporter pathway demo",
    species = "Human",
    cell_line = "HEK-Blue TLR8 demo cells",
    reporter = "NF-kB/AP-1 SEAP",
    readout_nm = 655,
    stringsAsFactors = FALSE
  )
  peptides <- data.frame(
    peptide_id = sprintf("PEP-%03d", 1:6),
    peptide_name = c("strong agonist-like peptide", "antagonist-like peptide", "artifact-prone peptide", "inactive peptide", "weak mixed-response peptide", "unknown-direction peptide"),
    expected_activity = c("agonist", "antagonist", "artifact", "unknown", "unknown", "unknown"),
    lot = c("P001-B", "P002-B", "P003-B", "P004-A", "P005-A", "P006-A"),
    purity_percent = c(96.2, 94.1, 88.4, 97.8, 91.6, 93.2),
    molecular_weight_da = c(1820, 1765, 1944, 1688, 2015, 1875),
    stock_concentration_mM = c(10, 10, 8, 10, 5, 10),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    field = c("scientist", "project", "assay_date", "target_id", "cell_line", "protocol_version", "incubation_hours", "readout_nm", "cell_passage", "cell_lot", "quanti_blue_lot", "instrument", "peptide_lot", "vehicle", "reader_settings", "notes"),
    value = c("Demo scientist", "HEKBlueR demo", "2026-07-31", target$target_id, target$cell_line, "v1.1", "18", "655", "P18", "CELL-DEMO-02", "QB-DEMO-08", "Demo plate reader", "mixed demo lots", "0.1% DMSO", "endpoint OD 655 nm", "Simulated multi-plate demo with primary, secondary, and counter-assay modules."),
    required = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  wells <- make_wells()
  concs <- c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
  plate_defs <- data.frame(
    plate_id = c("PRIMARY_AGO_B1", "PRIMARY_AGO_B2", "SECONDARY_AGO_CURVE_B1", "SECONDARY_ANT_CURVE_B1", "COUNTER_ARTIFACT_B1"),
    assay_stage = c("primary", "primary", "secondary", "secondary", "counter"),
    assay_mode = c("agonist", "agonist", "agonist", "antagonist", "counter"),
    biological_replicate = c(1, 2, 1, 1, 1),
    plate_shift = c(0.00, 0.05, -0.02, 0.03, 0.01),
    stringsAsFactors = FALSE
  )

  all_maps <- list()
  all_raw <- list()

  for (plate_i in seq_len(nrow(plate_defs))) {
    plate <- plate_defs[plate_i, ]
    map <- wells
    map$plate_id <- plate$plate_id
    map$assay_stage <- plate$assay_stage
    map$assay_mode <- plate$assay_mode
    map$target_id <- target$target_id
    map$sample_id <- "EMPTY"
    map$peptide_id <- NA_character_
    map$expected_activity <- NA_character_
    map$control_type <- "empty"
    map$concentration_uM <- NA_real_
    map$technical_replicate <- NA_integer_
    map$biological_replicate <- plate$biological_replicate

    assign_block <- function(row_name, col_ids, sample_id, peptide_id, control_type, conc = NA_real_, tech = NULL) {
      idx <- map$row == row_name & map$col %in% sprintf("%02d", col_ids)
      map$sample_id[idx] <<- sample_id
      map$peptide_id[idx] <<- peptide_id
      map$expected_activity[idx] <<- ifelse(is.na(peptide_id), NA_character_, peptides$expected_activity[match(peptide_id, peptides$peptide_id)])
      map$control_type[idx] <<- control_type
      map$concentration_uM[idx] <<- conc
      map$technical_replicate[idx] <<- if (is.null(tech)) seq_len(sum(idx)) else tech[seq_len(sum(idx))]
    }

    if (plate$assay_stage == "primary") {
      for (r in LETTERS[1:2]) assign_block(r, 1:2, "blank", NA, "blank")
      for (r in LETTERS[1:2]) assign_block(r, 3:4, "vehicle", NA, "negative_control")
      for (r in LETTERS[1:2]) assign_block(r, 5:6, "positive_agonist", NA, "positive_control")
      for (r in LETTERS[1:2]) assign_block(r, 7:8, "ipc_reference", "PEP-001", "inter_plate_calibrator", 1)
      layout <- data.frame(
        row = rep(LETTERS[3:8], each = 2),
        cols = I(rep(list(1:3, 4:6), 6)),
        peptide_id = rep(peptides$peptide_id, each = 2),
        concentration = rep(c(1, 10), times = 6),
        stringsAsFactors = FALSE
      )
      for (i in seq_len(nrow(layout))) {
        assign_block(layout$row[i], layout$cols[[i]], paste0(layout$peptide_id[i], "_primary_", layout$concentration[i], "uM"), layout$peptide_id[i], "test_sample", layout$concentration[i])
      }
      for (r in LETTERS[3:8]) assign_block(r, 10:12, "vehicle_control_block", NA, "negative_control")
    }

    if (plate$assay_stage == "secondary") {
      for (r in LETTERS[1:2]) assign_block(r, 1:2, "blank", NA, "blank")
      for (r in LETTERS[1:2]) assign_block(r, 3:4, "vehicle", NA, "negative_control")
      if (plate$assay_mode == "agonist") {
        for (r in LETTERS[1:2]) assign_block(r, 5:6, "positive_agonist", NA, "positive_control")
      } else {
        for (r in LETTERS[1:2]) assign_block(r, 5:6, "agonist_challenge", NA, "agonist_challenge_control")
        for (r in LETTERS[1:2]) assign_block(r, 7:8, "known_antagonist", NA, "known_antagonist_control", 1)
      }
      assign_block("A", 9:12, "ipc_reference", "PEP-001", "inter_plate_calibrator", 1)
      curve_peptides <- c("PEP-001", "PEP-002", "PEP-003")
      rows_for <- c("C", "D", "E")
      for (p in seq_along(curve_peptides)) {
        for (j in seq_along(concs)) {
          cols <- if (j <= 4) ((j - 1) * 3 + 1):((j - 1) * 3 + 3) else ((j - 5) * 3 + 1):min(((j - 5) * 3 + 3), 12)
          row <- if (j <= 4) rows_for[p] else LETTERS[match(rows_for[p], LETTERS) + 3]
          assign_block(row, cols, paste0(curve_peptides[p], "_", concs[j], "uM"), curve_peptides[p], "test_sample", concs[j])
        }
      }
    }

    if (plate$assay_stage == "counter") {
      for (r in LETTERS[1:2]) assign_block(r, 1:4, "blank", NA, "blank")
      for (r in LETTERS[1:2]) assign_block(r, 5:8, "vehicle", NA, "negative_control")
      for (r in LETTERS[1:2]) assign_block(r, 9:12, "positive_agonist", NA, "positive_control")
      rows_for <- LETTERS[3:8]
      for (p in seq_along(peptides$peptide_id)) {
        assign_block(rows_for[p], 1:3, paste0(peptides$peptide_id[p], "_viability"), peptides$peptide_id[p], "viability_counter", 10)
        assign_block(rows_for[p], 4:6, paste0(peptides$peptide_id[p], "_no_cell"), peptides$peptide_id[p], "no_cell_interference", 10)
        assign_block(rows_for[p], 7:9, paste0(peptides$peptide_id[p], "_unrelated_reporter"), peptides$peptide_id[p], "unrelated_reporter", 10)
        assign_block(rows_for[p], 10:12, paste0(peptides$peptide_id[p], "_null_cell"), peptides$peptide_id[p], "null_cell_reporter", 10)
      }
    }

    raw <- map[, c("plate_id", "well", "row", "col", "assay_stage", "assay_mode", "target_id", "sample_id", "peptide_id", "expected_activity", "control_type", "concentration_uM", "technical_replicate", "biological_replicate")]
    raw$raw_od <- vapply(seq_len(nrow(raw)), function(i) simulate_od(raw[i, ], plate$plate_shift, noise_sd = ifelse(plate$plate_id == "PRIMARY_AGO_B2", 0.02, 0.012)), numeric(1))
    raw$raw_od[raw$plate_id == "PRIMARY_AGO_B2" & raw$well == "H12"] <- 2.95
    raw$raw_od[raw$plate_id == "SECONDARY_ANT_CURVE_B1" & raw$well == "E11"] <- NA_real_

    all_maps[[plate$plate_id]] <- map
    all_raw[[plate$plate_id]] <- raw
  }

  samplesheet <- data.frame(
    run_id = "demo_multi_plate",
    raw = "data/simulated/raw_plate_reader.csv",
    plate_map = "data/simulated/plate_map.csv",
    metadata = "data/simulated/run_metadata.csv",
    stringsAsFactors = FALSE
  )

  list(
    target = target,
    peptides = peptides,
    metadata = metadata,
    plate_map = do.call(rbind, all_maps),
    raw_data = do.call(rbind, all_raw),
    samplesheet = samplesheet
  )
}

write_demo_dataset <- function(out_dir = "data/simulated", seed = 42) {
  demo <- simulate_hekblue_demo(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(demo$target, file.path(out_dir, "target_metadata.csv"), row.names = FALSE)
  write.csv(demo$peptides, file.path(out_dir, "peptide_metadata.csv"), row.names = FALSE)
  write.csv(demo$metadata, file.path(out_dir, "run_metadata.csv"), row.names = FALSE)
  write.csv(demo$plate_map, file.path(out_dir, "plate_map.csv"), row.names = FALSE)
  write.csv(demo$raw_data, file.path(out_dir, "raw_plate_reader.csv"), row.names = FALSE)
  write.csv(demo$samplesheet, file.path(out_dir, "samplesheet.csv"), row.names = FALSE, quote = FALSE)
  invisible(demo)
}
