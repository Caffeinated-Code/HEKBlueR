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
      z_prime_warn = 0,
      control_cv_warn_percent = 20,
      control_cv_fail_percent = 35,
      edge_effect_warn = 0.15,
      calibration_drift_warn_od = 0.15,
      calibration_drift_fail_od = 0.30,
      intraplate_cv_warn_percent = 20,
      spatial_bias_warn = 0.15,
      spatial_bias_fail = 0.30,
      outlier_rate_warn_percent = 5,
      outlier_rate_fail_percent = 15,
      reference_stability_pass_score = 80,
      reference_stability_warn_score = 50
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
      "Plate", "Plate", "Plate", "Plate", "Plate", "Plate", "Plate",
      "Reference control", "Reference control", "Reference control", "Reference control", "Reference control",
      "Primary", "Primary", "Primary",
      "Dose-response", "Dose-response", "Dose-response", "Dose-response", "Dose-response",
      "Counter-assay", "Counter-assay", "Counter-assay", "Counter-assay"
    ),
    metric = c(
      "Completeness pass", "Completeness warning",
      "Technical replicates", "Dose points pass", "Dose points warning", "Control wells",
      "Z-prime pass", "Z-prime warning", "Control CV warning", "Calibration drift warning", "Intra-plate CV warning", "Spatial bias warning", "Outlier rate warning",
      "Reference stability pass", "Reference stability warning", "Control CV fail", "Calibration drift fail", "Outlier rate fail",
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
      paste0(">", thresholds$plate$intraplate_cv_warn_percent, "%"),
      paste0(">", thresholds$plate$spatial_bias_warn, " OD ratio"),
      paste0(">", thresholds$plate$outlier_rate_warn_percent, "% wells"),
      paste0(">=", thresholds$plate$reference_stability_pass_score),
      paste0(thresholds$plate$reference_stability_warn_score, " to ", thresholds$plate$reference_stability_pass_score - 1),
      paste0(">", thresholds$plate$control_cv_fail_percent, "%"),
      paste0("absolute OD >", thresholds$plate$calibration_drift_fail_od),
      paste0(">", thresholds$plate$outlier_rate_fail_percent, "% wells"),
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
      "Within-plate controls are too variable for confident normalization.",
      "Edge, row, or column behavior suggests spatial artifacts.",
      "A high flagged-well fraction suggests plate handling, readout, or dispense problems.",
      "Combined control stability score is high enough for unflagged normalization.",
      "Combined control stability score supports use with a warning.",
      "Control variation is too high for stable reference use.",
      "Shared reference drift is too high for stable reference use.",
      "A high flagged-well fraction weakens reference confidence.",
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

qc_threshold_recommendations <- function() {
  data.frame(
    input_id = c(
      "z_prime_pass", "z_prime_warn", "control_cv_warn_percent", "edge_effect_warn",
      "calibration_drift_warn_od", "intraplate_cv_warn_percent", "spatial_bias_warn",
      "outlier_rate_warn_percent", "reference_stability_pass_score", "reference_stability_warn_score",
      "replicate_cv_warn_percent", "replicate_cv_fail_percent",
      "dose_points_pass", "dose_points_warn", "dynamic_range_warn_percent",
      "curve_replicate_cv_warn_percent", "hill_min", "hill_max",
      "viability_min_od", "no_cell_max_od", "unrelated_reporter_max_od", "null_cell_max_od"
    ),
    label = c(
      "Z-prime pass", "Z-prime warn floor", "Control CV warning (%)", "Edge effect warning",
      "Inter-plate drift warning (OD)", "Intra-plate CV warning (%)", "Spatial bias warning",
      "Outlier rate warning (%)", "Reference stability PASS score", "Reference stability WARN floor",
      "Primary replicate CV warning (%)", "Primary replicate CV fail (%)",
      "Dose points pass", "Dose points warning", "Dynamic range warning (%)",
      "Curve replicate CV warning (%)", "Hill slope minimum", "Hill slope maximum",
      "Viability minimum OD", "No-cell interference max OD", "Unrelated reporter max OD", "Null-cell reporter max OD"
    ),
    default = c(0.5, 0, 20, 0.15, 0.15, 20, 0.15, 5, 80, 50, 25, 50, 8, 5, 30, 30, 0.2, 4, 0.7, 0.25, 0.4, 0.35),
    safe_min = c(0.3, -0.5, 5, 0.05, 0.03, 5, 0.05, 1, 70, 35, 10, 25, 6, 4, 10, 10, 0.1, 2, 0.3, 0.05, 0.05, 0.05),
    safe_max = c(0.8, 0.4, 35, 0.35, 0.5, 35, 0.35, 15, 95, 70, 40, 80, 12, 7, 60, 50, 0.5, 6, 1.5, 1, 1, 1),
    rationale = c(
      "Classic HTS practice uses 0.5 as a strong assay-quality target, while complex cell assays may need contextual review.",
      "Values above 0 but below 0.5 are marginal but can be usable for cell-based biology.",
      "Assay Guidance Manual recommends control CV below 20 percent for max and min plates.",
      "Flags edge-center shift large enough to suggest evaporation or handling artifacts.",
      "Flags shared reference drift large enough to change normalized OD interpretation.",
      "Uses the same 20 percent control CV principle inside each plate.",
      "Flags row, column, or edge bias large enough to affect local wells.",
      "Flags a plate when many wells require manual review.",
      "Requires strong combined control performance before unflagged normalization.",
      "Allows flagged control use when the combined score remains reviewable.",
      "Primary screen hit calls need stricter replicate agreement than exploratory plots.",
      "Large replicate CV can make a dose-level hit uninterpretable.",
      "Eight or more concentrations is a practical default for secondary curve review.",
      "Five to seven concentrations can be reviewed with caution.",
      "Weak response windows make potency estimates unstable.",
      "Curve-level replicate noise above 30 percent needs review.",
      "Very shallow fitted slopes often indicate weak or unstable fits.",
      "Very steep fitted slopes often indicate fitting instability or non-standard biology.",
      "Low viability can confound reporter interpretation.",
      "Signal without cells suggests optical or reagent interference.",
      "Unrelated reporter activity suggests reporter or pathway artifact.",
      "Null-cell activity suggests target-independent biology."
    ),
    stringsAsFactors = FALSE
  )
}
