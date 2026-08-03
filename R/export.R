qc_report_lines <- function(results) {
  assay_id <- results$assay_manifest$assay_identifier[1]
  collect_flags <- function(df, status_col, label, flag_col = NULL) {
    if (is.null(df) || !nrow(df) || !status_col %in% names(df)) return(character())
    bad <- df[df[[status_col]] %in% c("WARN", "FAIL"), , drop = FALSE]
    if (!nrow(bad)) return(character())
    apply(bad, 1, function(row) {
      id_cols <- intersect(c("plate_id", "peptide_id", "target_id"), names(df))
      id <- paste(na.omit(unname(row[id_cols])), collapse = " | ")
      flags <- if (!is.null(flag_col) && flag_col %in% names(df)) row[[flag_col]] else row[["notes"]]
      paste0("- ", label, ": ", row[[status_col]], ifelse(nchar(id), paste0(" | ", id), ""), ifelse(!is.na(flags) && nzchar(flags), paste0(" | ", flags), ""))
    })
  }
  warnings <- c(
    collect_flags(results$design_qc, "design_status", "Design QC", "flagged_metrics"),
    collect_flags(results$plate_qc, "plate_qc_status", "Plate QC", "notes"),
    collect_flags(results$intraplate_variability_qc, "intraplate_status", "Intra-plate QC", "notes"),
    collect_flags(results$dose_response_qc, "curve_qc_status", "Dose-response QC", "curve_flags"),
    collect_flags(results$counter_assay_qc, "counter_qc_status", "Counter-assay QC", "artifact_flags"),
    collect_flags(results$final_qc_table, "final_status", "Final QC", "final_action")
  )
  if (!length(warnings)) warnings <- "- No WARN or FAIL QC findings were detected."
  c(
    paste0("# HEKBlueR QC Report: ", assay_id),
    "",
    "## Major WARN/FAIL Findings",
    warnings,
    "",
    "## Run Summary",
    paste0("- Assay identifier: ", assay_id),
    paste0("- Target: ", results$assay_manifest$target_id[1]),
    paste0("- Assay type: ", results$assay_manifest$assay_type[1]),
    paste0("- Peptides: ", results$assay_manifest$peptide_ids[1]),
    paste0("- Thresholds changed from default: ", results$assay_manifest$threshold_changed_from_default[1]),
    paste0("- Threshold change note: ", ifelse(nzchar(results$assay_manifest$threshold_change_note[1]), results$assay_manifest$threshold_change_note[1], "None")),
    "",
    "## Folder Guide",
    "- documentation: assay manifest, run documentation, metadata, active QC thresholds, and analysis config",
    "- tables: raw data, plate map, normalized results, QC tables, final calls, and exclusions",
    "- figures: review plots generated from this analysis"
  )
}

save_export_figures <- function(results, figure_dir) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  source_if_exists <- function(path) if (file.exists(path)) source(path)
  source_if_exists("R/plots.R")
  plots <- list(
    raw_od_heatmap = try(plate_heatmap_plot(results$cleaned_well_data, "raw_od", "Raw OD by well"), silent = TRUE),
    normalized_heatmap = try(plate_heatmap_plot(results$normalized_results, "percent_activation", "Percent activation by well"), silent = TRUE),
    zprime_qc = try(qc_bar_plot(results$plate_qc), silent = TRUE),
    raw_od_distribution = try(raw_distribution_plot(results$cleaned_well_data), silent = TRUE),
    primary_waterfall = try(waterfall_plot(results$primary_results), silent = TRUE),
    dose_response = try(dose_response_plot(results$primary_results, results$dose_response_results), silent = TRUE),
    replicate_cv = try(replicate_cv_plot(results$primary_results), silent = TRUE),
    control_behavior = try(control_boxplot(results$cleaned_well_data), silent = TRUE),
    edge_effect = try(edge_plot(results$cleaned_well_data), silent = TRUE),
    interplate_calibration = try(calibration_plot(results$interplate_calibration), silent = TRUE)
  )
  for (name in names(plots)) {
    if (!inherits(plots[[name]], "try-error")) {
      ggplot2::ggsave(file.path(figure_dir, paste0(name, ".png")), plot = plots[[name]], width = 11, height = 7, dpi = 180)
    }
  }
  invisible(TRUE)
}

export_hekblue_results <- function(results, raw_data, plate_map, metadata, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  write.csv(raw_data, file.path(output_dir, "raw_data.csv"), row.names = FALSE)
  write.csv(plate_map, file.path(output_dir, "plate_map.csv"), row.names = FALSE)
  write.csv(metadata, file.path(output_dir, "metadata.csv"), row.names = FALSE)
  for (name in names(results)) {
    if (is.data.frame(results[[name]])) {
      write.csv(results[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
    }
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      list(
        assay_manifest = results$assay_manifest,
        run_documentation = results$run_documentation,
        metadata = metadata,
        final_qc_table = results$final_qc_table,
        plate_qc = results$plate_qc
      ),
      file.path(output_dir, "run_summary.json"),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }
  writeLines(
    c(
      "assay_type: HEK-Blue SEAP reporter",
      "normalization: blank correction plus control scaling",
      "threshold_source: R/qc_thresholds.R",
      "curve_model: four_parameter_logistic_base_R",
      "outputs_are_database_ready: true"
    ),
    file.path(output_dir, "analysis_config.yml")
  )
  invisible(output_dir)
}

export_hekblue_results_structured <- function(results, raw_data, plate_map, metadata, output_root) {
  assay_id <- results$assay_manifest$assay_identifier[1]
  output_dir <- file.path(output_root, paste0(assay_id, "_results"))
  doc_dir <- file.path(output_dir, "documentation")
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(results$assay_manifest, file.path(doc_dir, "assay_manifest.csv"), row.names = FALSE)
  write.csv(results$run_documentation, file.path(doc_dir, "run_documentation.csv"), row.names = FALSE)
  write.csv(metadata, file.path(doc_dir, "metadata.csv"), row.names = FALSE)
  write.csv(results$qc_thresholds, file.path(doc_dir, "qc_thresholds.csv"), row.names = FALSE)
  writeLines(qc_report_lines(results), file.path(output_dir, "qc_report.md"))
  writeLines(
    c(
      "assay_type: HEK-Blue SEAP reporter",
      "normalization: blank correction plus control scaling",
      "threshold_source: R/qc_thresholds.R",
      "curve_model: four_parameter_logistic_base_R",
      "outputs_are_database_ready: true"
    ),
    file.path(doc_dir, "analysis_config.yml")
  )

  write.csv(raw_data, file.path(table_dir, "raw_data.csv"), row.names = FALSE)
  write.csv(plate_map, file.path(table_dir, "plate_map.csv"), row.names = FALSE)
  for (name in names(results)) {
    if (is.data.frame(results[[name]]) && !name %in% c("assay_manifest", "run_documentation", "qc_thresholds")) {
      write.csv(results[[name]], file.path(table_dir, paste0(name, ".csv")), row.names = FALSE)
    }
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      list(assay_manifest = results$assay_manifest, final_qc_table = results$final_qc_table, qc_report = qc_report_lines(results)),
      file.path(doc_dir, "run_summary.json"),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }
  save_export_figures(results, figure_dir)
  invisible(output_dir)
}
