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
      "z_prime_pass: 0.5",
      "control_cv_warn_percent: 20",
      "replicate_cv_warn_percent: 20",
      "curve_model: four_parameter_logistic_base_R",
      "outputs_are_database_ready: true"
    ),
    file.path(output_dir, "analysis_config.yml")
  )
  invisible(output_dir)
}

