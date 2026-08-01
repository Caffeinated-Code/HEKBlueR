make_wells <- function(format = 96) {
  if (format != 96) {
    stop("The demo generator currently supports 96-well plates.")
  }
  rows <- LETTERS[1:8]
  cols <- sprintf("%02d", 1:12)
  expand.grid(row = rows, col = cols, stringsAsFactors = FALSE) |>
    transform(well = paste0(row, col))
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
    peptide_id = c("PEP-001", "PEP-002", "PEP-003"),
    peptide_name = c("Demo agonist-like peptide", "Demo antagonist-like peptide", "Demo artifact-prone peptide"),
    expected_role = c("agonist", "antagonist", "artifact"),
    lot = c("P001-A", "P002-A", "P003-A"),
    purity_percent = c(95, 93, 88),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    field = c(
      "scientist", "project", "assay_date", "target_id", "cell_line",
      "cell_passage", "cell_lot", "quanti_blue_lot", "instrument",
      "protocol_version", "incubation_hours", "readout_nm"
    ),
    value = c(
      "Demo scientist", "HEKBlueR demo", "2026-07-31", target$target_id,
      target$cell_line, "P18", "CELL-DEMO-01", "QB-DEMO-07",
      "Demo plate reader", "v1.0", "18", "655"
    ),
    required = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  wells <- make_wells()
  concentrations <- c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)

  plate_defs <- data.frame(
    plate_id = c("PLATE_AGONIST_001", "PLATE_ANTAGONIST_001", "PLATE_COUNTER_001"),
    assay_mode = c("agonist", "antagonist", "counter"),
    batch_id = c("BATCH_001", "BATCH_001", "BATCH_001"),
    stringsAsFactors = FALSE
  )

  all_maps <- list()
  all_raw <- list()

  for (plate_i in seq_len(nrow(plate_defs))) {
    plate <- plate_defs[plate_i, ]
    map <- wells
    map$plate_id <- plate$plate_id
    map$assay_mode <- plate$assay_mode
    map$target_id <- target$target_id
    map$sample_id <- "EMPTY"
    map$peptide_id <- NA_character_
    map$control_type <- "empty"
    map$concentration_uM <- NA_real_
    map$technical_replicate <- NA_integer_
    map$biological_replicate <- NA_integer_

    assign_block <- function(row_name, col_ids, sample_id, peptide_id, control_type, conc = NA_real_, bio_rep = NA_integer_) {
      idx <- map$row == row_name & map$col %in% sprintf("%02d", col_ids)
      map$sample_id[idx] <<- sample_id
      map$peptide_id[idx] <<- peptide_id
      map$control_type[idx] <<- control_type
      map$concentration_uM[idx] <<- conc
      map$technical_replicate[idx] <<- seq_len(sum(idx))
      map$biological_replicate[idx] <<- bio_rep
    }

    if (plate$assay_mode == "agonist") {
      for (r in LETTERS[1:2]) assign_block(r, 1:4, "blank", NA, "blank")
      for (r in LETTERS[1:2]) assign_block(r, 5:8, "vehicle", NA, "negative_control")
      for (r in LETTERS[1:2]) assign_block(r, 9:12, "positive_agonist", NA, "positive_control")
      row_assign <- c("C", "D", "E")
      for (p in seq_along(peptides$peptide_id)) {
        for (j in seq_along(concentrations)) {
          cols <- if (j <= 4) ((j - 1) * 3 + 1):((j - 1) * 3 + 3) else ((j - 5) * 3 + 1):((j - 5) * 3 + 3)
          row <- if (j <= 4) row_assign[p] else LETTERS[match(row_assign[p], LETTERS) + 3]
          assign_block(row, cols, paste0(peptides$peptide_id[p], "_", concentrations[j], "uM"), peptides$peptide_id[p], "test_sample", concentrations[j], p)
        }
      }
      assign_block("H", 10:12, "ipc_reference", "PEP-001", "inter_plate_calibrator", 1, 1)
    }

    if (plate$assay_mode == "antagonist") {
      for (r in LETTERS[1:2]) assign_block(r, 1:4, "blank", NA, "blank")
      for (r in LETTERS[1:2]) assign_block(r, 5:8, "vehicle", NA, "negative_control")
      for (r in LETTERS[1:2]) assign_block(r, 9:12, "agonist_challenge", NA, "agonist_challenge_control")
      row_assign <- c("C", "D", "E")
      for (p in seq_along(peptides$peptide_id)) {
        for (j in seq_along(concentrations)) {
          cols <- if (j <= 4) ((j - 1) * 3 + 1):((j - 1) * 3 + 3) else ((j - 5) * 3 + 1):((j - 5) * 3 + 3)
          row <- if (j <= 4) row_assign[p] else LETTERS[match(row_assign[p], LETTERS) + 3]
          assign_block(row, cols, paste0(peptides$peptide_id[p], "_", concentrations[j], "uM_plus_agonist"), peptides$peptide_id[p], "test_sample", concentrations[j], p)
        }
      }
      assign_block("H", 9:12, "known_antagonist", NA, "known_antagonist_control", 1)
    }

    if (plate$assay_mode == "counter") {
      for (r in LETTERS[1:2]) assign_block(r, 1:4, "blank", NA, "blank")
      for (r in LETTERS[1:2]) assign_block(r, 5:8, "vehicle", NA, "negative_control")
      for (r in LETTERS[1:2]) assign_block(r, 9:12, "positive_agonist", NA, "positive_control")
      rows_for <- c("C", "D", "E")
      for (p in seq_along(peptides$peptide_id)) {
        assign_block(rows_for[p], 1:3, paste0(peptides$peptide_id[p], "_viability"), peptides$peptide_id[p], "viability_counter", 10, p)
        assign_block(rows_for[p], 4:6, paste0(peptides$peptide_id[p], "_no_cell"), peptides$peptide_id[p], "no_cell_interference", 10, p)
        assign_block(rows_for[p], 7:9, paste0(peptides$peptide_id[p], "_unrelated_reporter"), peptides$peptide_id[p], "unrelated_reporter", 10, p)
        assign_block(rows_for[p], 10:12, paste0(peptides$peptide_id[p], "_null_cell"), peptides$peptide_id[p], "null_cell_reporter", 10, p)
      }
      assign_block("H", 10:12, "ipc_reference", "PEP-001", "inter_plate_calibrator", 1, 1)
    }

    raw <- map[, c("plate_id", "well", "row", "col", "assay_mode", "sample_id", "peptide_id", "control_type", "concentration_uM", "technical_replicate", "biological_replicate")]
    raw$raw_od <- vapply(seq_len(nrow(raw)), function(i) {
      ct <- raw$control_type[i]
      pep <- raw$peptide_id[i]
      conc <- raw$concentration_uM[i]
      base <- switch(
        ct,
        blank = 0.08,
        negative_control = 0.18,
        positive_control = 1.35,
        agonist_challenge_control = 1.20,
        known_antagonist_control = 0.35,
        inter_plate_calibrator = 0.75,
        empty = NA_real_,
        0.20
      )
      if (ct == "test_sample" && plate$assay_mode == "agonist") {
        max_resp <- c("PEP-001" = 1.05, "PEP-002" = 0.35, "PEP-003" = 0.75)[pep]
        ec50 <- c("PEP-001" = 0.25, "PEP-002" = 1.4, "PEP-003" = 0.8)[pep]
        base <- 0.18 + max_resp * conc / (ec50 + conc)
      }
      if (ct == "test_sample" && plate$assay_mode == "antagonist") {
        max_inhib <- c("PEP-001" = 0.15, "PEP-002" = 0.82, "PEP-003" = 0.45)[pep]
        ic50 <- c("PEP-001" = 3.0, "PEP-002" = 0.18, "PEP-003" = 0.65)[pep]
        base <- 1.20 - max_inhib * conc / (ic50 + conc)
      }
      if (ct == "viability_counter") {
        base <- c("PEP-001" = 0.95, "PEP-002" = 0.90, "PEP-003" = 0.52)[pep]
      }
      if (ct == "no_cell_interference") {
        base <- c("PEP-001" = 0.08, "PEP-002" = 0.09, "PEP-003" = 0.38)[pep]
      }
      if (ct == "unrelated_reporter") {
        base <- c("PEP-001" = 0.18, "PEP-002" = 0.16, "PEP-003" = 0.62)[pep]
      }
      if (ct == "null_cell_reporter") {
        base <- c("PEP-001" = 0.14, "PEP-002" = 0.15, "PEP-003" = 0.47)[pep]
      }
      if (is.na(base)) return(NA_real_)
      row_bias <- ifelse(raw$row[i] %in% c("A", "H"), 0.025, 0)
      col_bias <- ifelse(raw$col[i] %in% c("01", "12"), 0.015, 0)
      round(base + row_bias + col_bias + rnorm(1, 0, 0.035), 4)
    }, numeric(1))

    all_maps[[plate$plate_id]] <- map
    all_raw[[plate$plate_id]] <- raw
  }

  list(
    target = target,
    peptides = peptides,
    metadata = metadata,
    plate_map = do.call(rbind, all_maps),
    raw_data = do.call(rbind, all_raw)
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
  invisible(demo)
}
