plate_heatmap_plot <- function(df, value_col = "raw_od", title = "Plate heatmap") {
  ggplot2::ggplot(df, ggplot2::aes(x = as.integer(col), y = row, fill = .data[[value_col]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_y_discrete(limits = rev(LETTERS[1:8])) +
    ggplot2::scale_x_continuous(breaks = 1:12) +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90") +
    ggplot2::facet_wrap(~ plate_id) +
    ggplot2::labs(x = "Column", y = "Row", fill = value_col, title = title) +
    ggplot2::theme_minimal(base_size = 12)
}

control_boxplot <- function(df) {
  controls <- df[df$control_type != "test_sample" & df$control_type != "empty", ]
  ggplot2::ggplot(controls, ggplot2::aes(x = control_type, y = blank_corrected_od, fill = control_type)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.5) +
    ggplot2::facet_wrap(~ plate_id, scales = "free_y") +
    ggplot2::labs(x = "Control", y = "Blank-corrected OD", title = "Control behavior") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
}

qc_bar_plot <- function(plate_qc) {
  ggplot2::ggplot(plate_qc, ggplot2::aes(x = plate_id, y = z_prime, fill = plate_qc_status)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "#1f7a4d") +
    ggplot2::geom_hline(yintercept = 0.3, linetype = "dotted", color = "#b45309") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Plate", y = "Z-prime", title = "Plate QC summary") +
    ggplot2::theme_minimal(base_size = 12)
}

waterfall_plot <- function(primary) {
  if (!nrow(primary)) return(ggplot2::ggplot() + ggplot2::theme_void())
  metric <- ifelse(primary$assay_mode == "antagonist", primary$inhibition_mean, primary$activation_mean)
  plot_df <- primary
  plot_df$response <- metric
  ggplot2::ggplot(plot_df, ggplot2::aes(x = stats::reorder(paste(peptide_id, concentration_uM, sep = " "), response), y = response, fill = primary_hit)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ assay_mode, scales = "free_x") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Peptide dose", y = "Normalized response", title = "Primary screen waterfall") +
    ggplot2::theme_minimal(base_size = 12)
}

dose_response_plot <- function(primary, dose_results) {
  if (!nrow(primary)) return(ggplot2::ggplot() + ggplot2::theme_void())
  plot_df <- primary
  plot_df$response <- ifelse(plot_df$assay_mode == "antagonist", plot_df$inhibition_mean, plot_df$activation_mean)
  ggplot2::ggplot(plot_df, ggplot2::aes(x = concentration_uM, y = response, color = peptide_id)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line() +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~ assay_mode) +
    ggplot2::labs(x = "Concentration (uM)", y = "Response", title = "Dose-response summary") +
    ggplot2::theme_minimal(base_size = 12)
}

edge_plot <- function(df) {
  z <- df
  z$plate_region <- ifelse(z$row %in% c("A", "H") | z$col %in% c("01", "12"), "edge", "center")
  ggplot2::ggplot(z, ggplot2::aes(x = plate_region, y = blank_corrected_od, fill = plate_region)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~ plate_id, scales = "free_y") +
    ggplot2::labs(x = "Plate region", y = "Blank-corrected OD", title = "Edge effect review") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")
}

