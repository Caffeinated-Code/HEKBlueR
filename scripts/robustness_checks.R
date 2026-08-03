suppressPackageStartupMessages({
  source("R/qc_metrics.R")
  source("R/qc_thresholds.R")
  source("R/export.R")
  source("R/analysis.R")
  source("R/simulate_data.R")
})

ok <- function(name, expr) {
  message("Checking: ", name)
  force(expr)
  message("  ok")
}

expect_error <- function(name, expr, pattern) {
  message("Checking: ", name)
  err <- tryCatch({
    force(expr)
    NULL
  }, error = function(e) conditionMessage(e))
  if (is.null(err)) stop(name, " did not fail as expected.", call. = FALSE)
  if (!grepl(pattern, err, ignore.case = TRUE)) stop(name, " failed with unexpected message: ", err, call. = FALSE)
  message("  ok: ", err)
}

demo <- simulate_hekblue_demo(seed = 2031)
raw <- demo$raw_data
plate_map <- demo$plate_map
metadata <- demo$metadata

ok("full multi-plate demo", {
  res <- run_hekblue_analysis(raw, plate_map, metadata)
  stopifnot(nrow(res$final_qc_table) == 6)
})

ok("small upload demo subset", {
  keep <- c("PRIMARY_AGO_B1", "SECONDARY_ANT_CURVE_B1", "COUNTER_ARTIFACT_B1")
  res <- run_hekblue_analysis(raw[raw$plate_id %in% keep, ], plate_map[plate_map$plate_id %in% keep, ], metadata)
  stopifnot(nrow(res$final_qc_table) >= 3)
})

ok("primary-only upload", {
  keep <- grepl("^PRIMARY", raw$plate_id)
  res <- run_hekblue_analysis(raw[keep, ], plate_map[plate_map$plate_id %in% unique(raw$plate_id[keep]), ], metadata)
  stopifnot(nrow(res$primary_results) > 0, nrow(res$dose_response_qc) == 0)
})

ok("secondary-only upload", {
  keep <- grepl("^SECONDARY", raw$plate_id)
  res <- run_hekblue_analysis(raw[keep, ], plate_map[plate_map$plate_id %in% unique(raw$plate_id[keep]), ], metadata)
  stopifnot(nrow(res$dose_response_qc) > 0)
})

ok("counter-only upload", {
  keep <- grepl("^COUNTER", raw$plate_id)
  res <- run_hekblue_analysis(raw[keep, ], plate_map[plate_map$plate_id %in% unique(raw$plate_id[keep]), ], metadata)
  stopifnot(nrow(res$counter_assay_qc) > 0, nrow(res$final_qc_table) > 0)
})

ok("raw-only annotated upload", {
  res <- run_hekblue_analysis(raw[raw$plate_id == "PRIMARY_AGO_B1", ], NULL, metadata)
  stopifnot(nrow(res$design_qc) == 1)
})

ok("column aliases and A1 well style", {
  alias <- raw[raw$plate_id == "PRIMARY_AGO_B1", c("plate_id", "well", "assay_stage", "assay_mode", "sample_id", "peptide_id", "expected_activity", "control_type", "concentration_uM", "technical_replicate", "biological_replicate", "raw_od")]
  names(alias) <- c("Plate", "Well Position", "Stage", "Direction", "Sample", "Compound ID", "Expected Activity", "Role", "Dose uM", "Tech Replicate", "Bio Replicate", "OD655")
  alias$`Well Position` <- ifelse(grepl("^[A-H]0[1-9]$", alias$`Well Position`), sub("0", "", alias$`Well Position`), alias$`Well Position`)
  res <- run_hekblue_analysis(alias, NULL, metadata)
  stopifnot(all(grepl("^[A-H][0-9]{2}$", res$cleaned_well_data$well)))
})

ok("wide one-row metadata", {
  wide <- data.frame(scientist = "Test scientist", project = "Robustness", assay_date = "2026-08-03", target_id = "TLR8_TEST")
  res <- run_hekblue_analysis(raw[raw$plate_id == "PRIMARY_AGO_B1", ], plate_map[plate_map$plate_id == "PRIMARY_AGO_B1", ], wide)
  stopifnot(res$assay_manifest$target_id == "TLR8_TEST")
})

ok("missing blank controls warns instead of crashing", {
  no_blank_raw <- raw[raw$plate_id == "PRIMARY_AGO_B1" & raw$control_type != "blank", ]
  no_blank_map <- plate_map[plate_map$plate_id == "PRIMARY_AGO_B1" & plate_map$control_type != "blank", ]
  res <- run_hekblue_analysis(no_blank_raw, no_blank_map, metadata)
  stopifnot(any(res$design_qc$design_status == "FAIL"))
})

ok("nonnumeric OD becomes missing QC", {
  bad <- raw[raw$plate_id == "PRIMARY_AGO_B1", ]
  bad$raw_od <- as.character(bad$raw_od)
  bad$raw_od[1] <- "OVER"
  res <- run_hekblue_analysis(bad, plate_map[plate_map$plate_id == "PRIMARY_AGO_B1", ], metadata)
  stopifnot(any(res$cleaned_well_data$missing_flag))
})

expect_error("duplicate well rows", {
  dup <- rbind(raw[raw$plate_id == "PRIMARY_AGO_B1", ], raw[raw$plate_id == "PRIMARY_AGO_B1", ][1, ])
  run_hekblue_analysis(dup, plate_map[plate_map$plate_id == "PRIMARY_AGO_B1", ], metadata)
}, "duplicate")

expect_error("invalid well IDs", {
  bad <- raw[raw$plate_id == "PRIMARY_AGO_B1", ]
  bad$well[1] <- "Z99"
  run_hekblue_analysis(bad, plate_map[plate_map$plate_id == "PRIMARY_AGO_B1", ], metadata)
}, "invalid well")

expect_error("missing raw OD column", {
  bad <- raw[raw$plate_id == "PRIMARY_AGO_B1", setdiff(names(raw), "raw_od")]
  run_hekblue_analysis(bad, plate_map[plate_map$plate_id == "PRIMARY_AGO_B1", ], metadata)
}, "raw_od")

message("All robustness checks passed.")
