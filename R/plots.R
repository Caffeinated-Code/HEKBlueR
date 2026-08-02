hek_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", color = "#17324d", margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(color = "#263849"),
      axis.text = ggplot2::element_text(color = "#334155"),
      strip.text = ggplot2::element_text(face = "bold", color = "#17324d"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.key.height = ggplot2::unit(0.55, "cm"),
      plot.margin = ggplot2::margin(14, 18, 14, 14)
    )
}

plate_heatmap_plot <- function(df, value_col = "raw_od", title = "Plate heatmap") {
  plot_df <- df
  plot_df$hover_text <- paste0(
    "Plate: ", plot_df$plate_id,
    "<br>Well: ", plot_df$well,
    "<br>Sample: ", plot_df$sample_id,
    "<br>Peptide/compound: ", plot_df$peptide_id,
    "<br>Control type: ", plot_df$control_type,
    "<br>Assay mode: ", plot_df$assay_mode,
    "<br>Concentration uM: ", plot_df$concentration_uM,
    "<br>", value_col, ": ", round(plot_df[[value_col]], 2)
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = as.integer(col), y = row, fill = .data[[value_col]], text = hover_text)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45) +
    ggplot2::scale_y_discrete(limits = rev(LETTERS[1:8])) +
    ggplot2::scale_x_continuous(breaks = 1:12) +
    ggplot2::scale_fill_viridis_c(option = "D", na.value = "grey92", guide = ggplot2::guide_colorbar(barheight = ggplot2::unit(4, "cm"), barwidth = ggplot2::unit(0.45, "cm"))) +
    ggplot2::facet_wrap(~ plate_id, ncol = 1) +
    ggplot2::coord_fixed(ratio = 1.05) +
    ggplot2::labs(x = "Column", y = "Row", fill = value_col, title = title) +
    hek_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12, color = "#16211f", face = "bold"),
      axis.title = ggplot2::element_text(size = 13, color = "#16211f", face = "bold"),
      strip.text = ggplot2::element_text(size = 12, face = "bold", color = "#16211f")
    )
}

plate_heatmap_plotly <- function(df, value_col = "raw_od", title = "Plate heatmap") {
  if (!requireNamespace("plotly", quietly = TRUE)) stop("plotly is required for interactive plate heatmaps")
  rows <- LETTERS[1:8]
  cols <- 1:12
  plates <- unique(df$plate_id)
  plots <- lapply(plates, function(pid) {
    plate <- df[df$plate_id == pid, , drop = FALSE]
    z <- matrix(NA_real_, nrow = length(rows), ncol = length(cols), dimnames = list(rows, cols))
    text <- matrix("", nrow = length(rows), ncol = length(cols), dimnames = list(rows, cols))
    for (i in seq_len(nrow(plate))) {
      r <- as.character(plate$row[i])
      c <- as.integer(plate$col[i])
      if (!r %in% rows || is.na(c) || !c %in% cols) next
      val <- plate[[value_col]][i]
      z[r, as.character(c)] <- val
      text[r, as.character(c)] <- paste0(
        "Plate: ", plate$plate_id[i],
        "<br>Well: ", plate$well[i],
        "<br>Sample: ", plate$sample_id[i],
        "<br>Peptide/compound: ", plate$peptide_id[i],
        "<br>Control type: ", plate$control_type[i],
        "<br>Assay mode: ", plate$assay_mode[i],
        "<br>Concentration uM: ", plate$concentration_uM[i],
        "<br>", value_col, ": ", round(val, 2)
      )
    }
    y_order <- rows
    plotly::plot_ly(
      x = cols,
      y = y_order,
      z = z[y_order, , drop = FALSE],
      text = text[y_order, , drop = FALSE],
      type = "heatmap",
      colors = c("#440154", "#31688e", "#35b779", "#fde725"),
      hoverinfo = "text",
      colorbar = list(title = value_col, len = 0.7)
    ) |>
      plotly::layout(
        title = list(text = pid, font = list(size = 15)),
        xaxis = list(title = "Column", tickmode = "array", tickvals = cols, ticktext = sprintf("%02d", cols), side = "bottom", tickfont = list(size = 13), titlefont = list(size = 14)),
        yaxis = list(title = "Row", tickmode = "array", tickvals = y_order, ticktext = y_order, tickfont = list(size = 13), titlefont = list(size = 14), autorange = "reversed"),
        margin = list(l = 70, r = 40, t = 55, b = 55)
      )
  })
  plotly::subplot(plots, nrows = length(plots), shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE, margin = 0.04) |>
    plotly::layout(title = list(text = title, font = list(size = 18)), showlegend = FALSE)
}

control_boxplot <- function(df) {
  controls <- df[df$control_type != "test_sample" & df$control_type != "empty", ]
  ggplot2::ggplot(controls, ggplot2::aes(x = control_type, y = blank_corrected_od, fill = control_type)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.5) +
    ggplot2::facet_wrap(~ plate_id, scales = "free_y") +
    ggplot2::labs(x = "Control", y = "Blank-corrected OD", title = "Control behavior") +
    hek_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none")
}

qc_bar_plot <- function(plate_qc) {
  ggplot2::ggplot(plate_qc, ggplot2::aes(x = plate_id, y = z_prime, fill = plate_qc_status)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "#1f7a4d") +
    ggplot2::geom_hline(yintercept = 0.3, linetype = "dotted", color = "#b45309") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Plate", y = "Z-prime", title = "Plate QC summary") +
    ggplot2::scale_fill_manual(values = c(PASS = "#2f8f6b", WARN = "#d89b2b", FAIL = "#c94c4c")) +
    hek_theme()
}

waterfall_plot <- function(primary) {
  if (!nrow(primary)) return(ggplot2::ggplot() + ggplot2::theme_void())
  metric <- ifelse(primary$assay_mode == "antagonist", primary$inhibition_mean, primary$activation_mean)
  plot_df <- primary
  plot_df$response <- metric
  plot_df$label <- paste(plot_df$peptide_id, paste0(plot_df$concentration_uM, " uM"), sep = " | ")
  ggplot2::ggplot(plot_df, ggplot2::aes(x = stats::reorder(label, response), y = response, fill = primary_hit)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(AGONIST_HIT = "#2f8f6b", ANTAGONIST_HIT = "#3867b7", NO_HIT = "#9aa7b2")) +
    ggplot2::facet_wrap(~ assay_mode, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Peptide dose", y = "Normalized response", title = "Primary screen waterfall") +
    hek_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}

dose_response_plot <- function(primary, dose_results) {
  if (!nrow(primary)) return(ggplot2::ggplot() + ggplot2::theme_void())
  plot_df <- primary
  plot_df$response <- ifelse(plot_df$assay_mode == "antagonist", plot_df$inhibition_mean, plot_df$activation_mean)
  fit_df <- data.frame()
  if (!is.null(dose_results) && nrow(dose_results)) {
    usable <- dose_results[is.finite(dose_results$bottom) & is.finite(dose_results$top) & is.finite(dose_results$ec50_ic50) & is.finite(dose_results$hill), ]
    fit_df <- do.call(rbind, lapply(seq_len(nrow(usable)), function(i) {
      row <- usable[i, ]
      source <- plot_df[plot_df$plate_id == row$plate_id & plot_df$assay_mode == row$assay_mode & plot_df$peptide_id == row$peptide_id, ]
      if (!nrow(source)) return(data.frame())
      xgrid <- exp(seq(log(min(source$concentration_uM, na.rm = TRUE)), log(max(source$concentration_uM, na.rm = TRUE)), length.out = 120))
      ygrid <- row$bottom + (row$top - row$bottom) / (1 + (row$ec50_ic50 / xgrid)^row$hill)
      data.frame(plate_id = row$plate_id, assay_mode = row$assay_mode, peptide_id = row$peptide_id, concentration_uM = xgrid, response = ygrid)
    }))
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = concentration_uM, y = response, color = peptide_id)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~ assay_mode) +
    ggplot2::labs(x = "Concentration (uM)", y = "Response", title = "Dose-response summary") +
    hek_theme()
  if (nrow(fit_df)) {
    p <- p + ggplot2::geom_line(data = fit_df, linewidth = 0.9)
  } else {
    p <- p + ggplot2::geom_line()
  }
  p
}

edge_plot <- function(df) {
  z <- df
  z$plate_region <- ifelse(z$row %in% c("A", "H") | z$col %in% c("01", "12"), "edge", "center")
  ggplot2::ggplot(z, ggplot2::aes(x = plate_region, y = blank_corrected_od, fill = plate_region)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~ plate_id, scales = "free_y") +
    ggplot2::labs(x = "Plate region", y = "Blank-corrected OD", title = "Edge effect review") +
    hek_theme() +
    ggplot2::theme(legend.position = "none")
}

raw_distribution_plot <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = raw_od, fill = assay_mode)) +
    ggplot2::geom_histogram(bins = 32, alpha = 0.8, color = "white") +
    ggplot2::facet_wrap(~ plate_id, scales = "free_y") +
    ggplot2::labs(x = "Raw OD", y = "Well count", title = "Raw OD distribution") +
    hek_theme()
}

replicate_cv_plot <- function(primary) {
  if (!nrow(primary)) return(ggplot2::ggplot() + ggplot2::theme_void())
  z <- primary
  z$cv <- ifelse(z$assay_mode == "antagonist", z$inhibition_cv, z$activation_cv)
  z$cv_display <- pmin(z$cv, 100)
  ggplot2::ggplot(z, ggplot2::aes(x = concentration_uM, y = cv_display, color = peptide_id)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "#d89b2b") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~ assay_mode) +
    ggplot2::labs(x = "Concentration (uM)", y = "Technical replicate CV (%), capped at 100 for display", title = "Replicate noise by dose") +
    hek_theme()
}

calibration_plot <- function(calibration) {
  if (!nrow(calibration)) return(ggplot2::ggplot() + ggplot2::theme_void())
  ggplot2::ggplot(calibration, ggplot2::aes(x = plate_id, y = calibration_factor, fill = calibration_status)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = c(-0.15, 0.15), linetype = "dashed", color = "#d89b2b") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(PASS = "#2f8f6b", WARN = "#d89b2b", FAIL = "#c94c4c")) +
    ggplot2::labs(x = "Plate", y = "Calibration factor", title = "Inter-plate calibration drift") +
    hek_theme()
}
