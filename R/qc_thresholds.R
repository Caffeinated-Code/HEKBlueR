default_qc_thresholds <- function() {
  list(
    metadata = list(
      completeness_pass_percent = 80,
      completeness_warn_percent = 50
    ),
    design = list(
      technical_replicates_pass = 3,
      biological_replicates_pass = 3,
      dose_points_pass = 8,
      dose_points_warn = 5,
      control_wells_pass = 8,
      control_wells_fail = 4
    ),
    plate = list(
      z_prime_pass = 0.5,
      z_prime_warn = 0.3,
      control_cv_warn_percent = 20,
      edge_effect_warn = 0.15,
      calibration_drift_warn_od = 0.15
    ),
    primary = list(
      agonist_hit_percent = 50,
      antagonist_hit_percent = 50,
      replicate_cv_warn_percent = 25,
      replicate_cv_fail_percent = 50
    ),
    curve = list(
      dose_points_pass = 8,
      dose_points_warn = 5,
      dynamic_range_warn_percent = 30,
      replicate_cv_warn_percent = 30,
      hill_min = 0.2,
      hill_max = 4,
      plateau_tolerance_percent = 20,
      monotonic_violations_warn = 1
    ),
    counter = list(
      viability_min_od = 0.7,
      no_cell_max_od = 0.25,
      unrelated_reporter_max_od = 0.4,
      null_cell_max_od = 0.35
    )
  )
}

qc_threshold_table <- function(thresholds = default_qc_thresholds()) {
  data.frame(
    qc_layer = c(
      "Metadata", "Metadata",
      "Design", "Design", "Design", "Design",
      "Plate", "Plate", "Plate", "Plate",
      "Primary", "Primary", "Primary",
      "Dose-response", "Dose-response", "Dose-response", "Dose-response", "Dose-response",
      "Counter-assay", "Counter-assay", "Counter-assay", "Counter-assay"
    ),
    metric = c(
      "Completeness pass", "Completeness warning",
      "Technical replicates", "Dose points pass", "Dose points warning", "Control wells",
      "Z-prime pass", "Z-prime warning", "Control CV warning", "Calibration drift warning",
      "Agonist hit", "Antagonist hit", "Replicate CV warning",
      "Dose points pass", "Dynamic range warning", "Replicate CV warning", "Hill slope range", "Monotonic reversals",
      "Viability minimum", "No-cell interference maximum", "Unrelated reporter maximum", "Null-cell reporter maximum"
    ),
    threshold = c(
      paste0(">=", thresholds$metadata$completeness_pass_percent, "%"),
      paste0(">=", thresholds$metadata$completeness_warn_percent, "%"),
      paste0(">=", thresholds$design$technical_replicates_pass),
      paste0(">=", thresholds$design$dose_points_pass),
      paste0(thresholds$design$dose_points_warn, " to ", thresholds$design$dose_points_pass - 1),
      paste0("PASS >= ", thresholds$design$control_wells_pass, "; FAIL < ", thresholds$design$control_wells_fail),
      paste0(">=", thresholds$plate$z_prime_pass),
      paste0(thresholds$plate$z_prime_warn, " to ", thresholds$plate$z_prime_pass),
      paste0(">", thresholds$plate$control_cv_warn_percent, "%"),
      paste0("absolute OD >", thresholds$plate$calibration_drift_warn_od),
      paste0("activation >= ", thresholds$primary$agonist_hit_percent, "%"),
      paste0("inhibition >= ", thresholds$primary$antagonist_hit_percent, "%"),
      paste0(">", thresholds$primary$replicate_cv_warn_percent, "%"),
      paste0(">=", thresholds$curve$dose_points_pass),
      paste0("< ", thresholds$curve$dynamic_range_warn_percent, "%"),
      paste0(">", thresholds$curve$replicate_cv_warn_percent, "%"),
      paste0(thresholds$curve$hill_min, " to ", thresholds$curve$hill_max),
      paste0(">", thresholds$curve$monotonic_violations_warn),
      paste0("< ", thresholds$counter$viability_min_od, " OD"),
      paste0("> ", thresholds$counter$no_cell_max_od, " OD"),
      paste0("> ", thresholds$counter$unrelated_reporter_max_od, " OD"),
      paste0("> ", thresholds$counter$null_cell_max_od, " OD")
    ),
    interpretation = c(
      "Run documentation is complete enough for reproducible review.",
      "Run documentation is sparse and should be improved.",
      "Technical replicate structure supports stable summary statistics.",
      "Dose series is strong enough for curve review.",
      "Dose series can be reviewed, but potency estimates are less stable.",
      "Required controls have enough wells for reliable control summaries.",
      "Positive and negative controls are well separated.",
      "Control separation is marginal and needs review.",
      "Control variation may reflect pipetting, reagent, timing, or reader issues.",
      "Shared reference indicates plate-to-plate drift.",
      "Primary activation is large enough to call a screen hit if replicate noise is acceptable.",
      "Primary inhibition is large enough to call a screen hit if replicate noise is acceptable.",
      "Dose-level replicate noise may weaken hit confidence.",
      "Curve has enough tested concentrations for secondary review.",
      "Observed response window is weak.",
      "Technical replicate noise is high at one or more doses.",
      "Extreme slopes may indicate poor fit, noisy data, or non-standard behavior.",
      "Dose trend is not consistent with expected response direction.",
      "Low viability can confound reporter interpretation.",
      "Signal in no-cell wells suggests assay interference.",
      "Activity in unrelated reporter suggests reporter or pathway artifact.",
      "Activity in null cells suggests target-independent biology."
    ),
    stringsAsFactors = FALSE
  )
}
