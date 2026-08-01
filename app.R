suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(plotly)
})

source("R/qc_metrics.R")
source("R/export.R")
source("R/analysis.R")
source("R/plots.R")
source("R/simulate_data.R")

APP_VERSION <- "0.2.0"

metric_help <- list(
  "Z-prime" = "Z-prime measures separation between positive and negative controls. Values above 0.5 are preferred. Values from 0.3 to 0.5 need review. Lower values usually mean the assay window is weak.",
  "Robust Z-prime" = "Robust Z-prime uses medians and median absolute deviation. It is less sensitive to individual outlier wells.",
  "SSMD" = "SSMD is a standardized difference between control groups. Larger absolute values indicate stronger separation.",
  "Control CV" = "Control coefficient of variation measures control stability. Values above 20 percent often indicate pipetting, reagent, incubation, or reader variability.",
  "Replicate CV" = "Replicate CV measures technical replicate noise at a peptide dose. High CV can make hit calls and curve fits unreliable.",
  "Edge effect" = "Edge effect compares edge wells to center wells. A large difference can reflect evaporation, temperature gradients, or plate handling effects.",
  "Inter-plate calibration" = "Inter-plate calibration aligns plates using a shared calibrator or shared positive control. Drift above 0.15 OD is flagged for review.",
  "EC50 or IC50" = "EC50 is the concentration giving half-maximal activation. IC50 is the concentration giving half-maximal inhibition. Estimates are strongest when they sit inside the tested dose range.",
  "Hill slope" = "Hill slope describes curve steepness. Very shallow or very steep slopes can indicate weak biology, noisy data, or fitting instability.",
  "Metadata completeness" = "Metadata completeness tracks how much run documentation was supplied. Complete metadata makes results easier to reproduce, search, and audit."
)

metric_link <- function(label) {
  actionLink(paste0("help_", gsub("[^A-Za-z0-9]", "_", label)), label, class = "metric-help")
}

read_uploaded_csv <- function(file_input, demo_path) {
  if (is.null(file_input)) read.csv(demo_path, stringsAsFactors = FALSE) else read.csv(file_input$datapath, stringsAsFactors = FALSE)
}

status_badge <- function(status) {
  color <- switch(status, PASS = "ok", WARN = "warn", FAIL = "fail", "neutral")
  span(class = paste("status-badge", color), status)
}

required_header <- function(label, required = TRUE) {
  tagList(label, span(class = ifelse(required, "field-badge required", "field-badge optional"), ifelse(required, "Required", "Optional")))
}

clean_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x), x, "hekblue")
}

status_datatable <- function(df, status_col = NULL, page_length = 15) {
  if (is.null(df) || !nrow(df)) return(datatable(data.frame(Message = "No rows available.")))
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  dt <- datatable(
    df,
    extensions = "FixedHeader",
    filter = "top",
    options = list(pageLength = page_length, scrollX = TRUE, scrollY = "430px", fixedHeader = TRUE)
  )
  if (length(numeric_cols)) dt <- formatRound(dt, numeric_cols, digits = 4)
  if (!is.null(status_col) && status_col %in% names(df)) {
    dt <- formatStyle(
      dt,
      status_col,
      target = "row",
      backgroundColor = styleEqual(c("PASS", "WARN", "FAIL"), c("#dff3ea", "#fff1cc", "#f9d7d7")),
      color = styleEqual(c("PASS", "WARN", "FAIL"), c("#143f2d", "#593a00", "#5b1111"))
    )
  }
  dt
}

download_plot <- function(plot_expr, prefix, active_inputs) {
  downloadHandler(
    filename = function() {
      dat <- active_inputs()
      meta <- if (is.null(dat)) data.frame(field = character(), value = character()) else dat$metadata
      target <- meta$value[match("target_id", meta$field)]
      project <- meta$value[match("project", meta$field)]
      paste(clean_name(paste(project, target, prefix, Sys.Date(), sep = "_")), ".png", sep = "")
    },
    content = function(file) {
      ggplot2::ggsave(file, plot = plot_expr(), width = 11, height = 7, dpi = 180)
    }
  )
}

schema_table <- data.frame(
  Column = c("plate_id", "well", "assay_mode", "sample_id", "peptide_id", "control_type", "concentration_uM", "raw_od", "technical_replicate", "biological_replicate"),
  Requirement = c("Required", "Required", "Required", "Required", "Required for samples", "Required", "Required for test samples", "Required", "Recommended", "Recommended"),
  Meaning = c("Unique plate identifier", "Well coordinate such as A01", "agonist, antagonist, or counter", "Sample or control label", "Peptide identifier", "Well role", "Peptide dose in uM", "Raw optical density", "Technical replicate number", "Biological replicate number"),
  stringsAsFactors = FALSE
)

ui <- page_navbar(
  title = "HEKBlueR",
  theme = bs_theme(
    version = 5,
    primary = "#215f7a",
    secondary = "#6a7583",
    success = "#2f8f6b",
    warning = "#d89b2b",
    danger = "#c94c4c",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter")
  ),
  header = tags$head(tags$style(HTML("
    body { background: #f4f7f9; color: #1f2d3a; }
    .navbar { box-shadow: 0 2px 10px rgba(23,50,77,0.12); }
    .bslib-card { border: 1px solid #d7e2ea; border-radius: 8px; box-shadow: 0 1px 5px rgba(23,50,77,0.05); }
    .card-header { background: #ffffff; color: #17324d; font-weight: 750; }
    .metric-card { border: 1px solid #d7e2ea; border-radius: 8px; padding: 14px; background: #ffffff; min-height: 108px; }
    .metric-label { font-size: 0.84rem; color: #536577; margin-bottom: 4px; }
    .metric-value { font-size: 1.45rem; font-weight: 760; color: #17324d; }
    .small-note { font-size: 0.88rem; color: #5d6b7a; }
    .status-badge { display: inline-block; border-radius: 999px; padding: 0.22rem 0.62rem; font-weight: 750; font-size: 0.78rem; }
    .status-badge.ok { background: #dff3ea; color: #143f2d; }
    .status-badge.warn { background: #fff1cc; color: #593a00; }
    .status-badge.fail { background: #f9d7d7; color: #5b1111; }
    .status-badge.neutral { background: #e7edf2; color: #263849; }
    .field-badge { margin-left: 8px; border-radius: 999px; padding: 2px 8px; font-size: 0.72rem; font-weight: 750; }
    .field-badge.required { background: #dff3ea; color: #143f2d; }
    .field-badge.optional { background: #e7edf2; color: #314155; }
    .metric-help { font-weight: 750; color: #215f7a; text-decoration: none; }
    .download-row { margin-top: 8px; display: flex; gap: 8px; flex-wrap: wrap; }
    .plot-card .html-widget { width: 100% !important; }
    @media (max-width: 900px) {
      .container-fluid { padding-left: 10px; padding-right: 10px; }
      .bslib-grid { grid-template-columns: 1fr !important; }
      .metric-card { min-height: auto; }
    }
  "))),
  nav_panel(
    "Start",
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Automated HEK-Blue screening review"),
        p("Upload raw plate-reader files, plate maps, and metadata. HEKBlueR runs design QC, plate QC, normalization, inter-plate calibration, dose-response review, counter-assay review, and final decision tables."),
        tags$ol(
          tags$li("Load demo data or upload your files."),
          tags$li("Preview uploaded data and metadata."),
          tags$li("Review EDA and cleaned data."),
          tags$li("Inspect QC, plots, and final decisions."),
          tags$li("Export database-ready results.")
        ),
        actionButton("load_demo", "Load demo data", class = "btn-primary"),
        span(" "),
        actionButton("run_analysis", "Run analysis", class = "btn-success")
      ),
      card(
        card_header("Demo data downloads"),
        p("Download these files to see the expected format."),
        div(class = "download-row",
          downloadButton("download_demo_raw", "Raw OD CSV"),
          downloadButton("download_demo_plate_map", "Plate map CSV"),
          downloadButton("download_demo_metadata", "Metadata CSV"),
          downloadButton("download_demo_peptides", "Peptide metadata CSV"),
          downloadButton("download_demo_target", "Target metadata CSV")
        ),
        hr(),
        uiOutput("run_status")
      )
    )
  ),
  nav_panel(
    "Upload",
    card(
      card_header("Expected raw data schema"),
      DTOutput("schema_table")
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      card(card_header(required_header("Raw plate-reader CSV", TRUE)), p(class = "small-note", "Needs one row per well. Required columns are listed above."), fileInput("raw_file", NULL)),
      card(card_header(required_header("Plate map CSV", TRUE)), p(class = "small-note", "Defines well role, sample, peptide, concentration, and replicate structure."), fileInput("plate_map_file", NULL)),
      card(card_header(required_header("Run metadata CSV", FALSE)), p(class = "small-note", "Optional if using the form below. Strongly recommended for reproducibility."), fileInput("metadata_file", NULL))
    ),
    card(
      card_header("Metadata form"),
      layout_columns(
        textInput("scientist", "Scientist name (Required)", ""),
        textInput("project", "Project (Required)", ""),
        dateInput("assay_date", "Assay date (Required)", value = Sys.Date()),
        textInput("target_id", "Target ID (Required)", "TARGET_TLR8_DEMO"),
        textInput("cell_line", "Cell line (Required)", ""),
        textInput("protocol_version", "Protocol version (Required)", "v1.0"),
        numericInput("incubation_hours", "Incubation hours (Required)", value = 18, min = 0),
        numericInput("readout_nm", "Readout wavelength (Required)", value = 655, min = 400, max = 800),
        textInput("cell_passage", "Cell passage (Optional)", ""),
        textInput("cell_lot", "Cell lot (Optional)", ""),
        textInput("reagent_lot", "QUANTI-Blue lot (Optional)", ""),
        textInput("instrument", "Instrument (Optional)", "")
      ),
      textAreaInput("run_notes", "Notes and protocol deviations (Optional)", "", height = "90px")
    )
  ),
  nav_panel(
    "Uploaded Preview",
    card(card_header("Raw data preview"), DTOutput("raw_preview")),
    card(card_header("Plate map preview"), DTOutput("plate_map_preview")),
    card(card_header("Metadata preview"), DTOutput("metadata_preview"))
  ),
  nav_panel(
    "EDA",
    layout_columns(
      col_widths = c(4, 8),
      card(card_header("Input EDA checks"), DTOutput("eda_table")),
      card(class = "plot-card", card_header("Raw OD distribution"), plotlyOutput("raw_distribution", height = "500px"), div(class = "download-row", downloadButton("download_raw_distribution", "Download raw OD plot")))
    )
  ),
  nav_panel(
    "Cleaned Data",
    card(card_header("Cleaning summary"), DTOutput("cleaning_summary")),
    card(card_header("Cleaned well data"), DTOutput("cleaned_table"))
  ),
  nav_panel(
    "Metadata",
    layout_columns(
      col_widths = c(4, 8),
      card(card_header(tagList("Completeness ", metric_link("Metadata completeness"))), uiOutput("metadata_card")),
      card(card_header("Submitted metadata"), DTOutput("metadata_table"))
    )
  ),
  nav_panel("Design QC", card(card_header("Experimental design review"), DTOutput("design_qc_table"))),
  nav_panel(
    "Plate QC",
    card(card_header(tagList("QC metrics: ", metric_link("Z-prime"), " | ", metric_link("Robust Z-prime"), " | ", metric_link("SSMD"), " | ", metric_link("Control CV"), " | ", metric_link("Edge effect"))), DTOutput("plate_qc_table")),
    card(class = "plot-card", card_header("Z-prime summary"), plotlyOutput("qc_plot", height = "430px"), div(class = "download-row", downloadButton("download_qc_plot", "Download Z-prime plot")))
  ),
  nav_panel(
    "Reference & Calibration",
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Reference control stability"), DTOutput("reference_qc_table")),
      card(card_header(tagList("Inter-plate calibration ", metric_link("Inter-plate calibration"))), DTOutput("calibration_table"))
    ),
    card(class = "plot-card", card_header("Calibration drift"), plotlyOutput("calibration_plot", height = "420px"), div(class = "download-row", downloadButton("download_calibration_plot", "Download calibration plot")))
  ),
  nav_panel(
    "Plate Layout",
    card(class = "plot-card", card_header("Raw OD heatmap"), plotlyOutput("raw_heatmap", height = "760px"), div(class = "download-row", downloadButton("download_raw_heatmap", "Download raw heatmap"))),
    card(class = "plot-card", card_header("Normalized heatmap"), plotlyOutput("normalized_heatmap", height = "760px"), div(class = "download-row", downloadButton("download_norm_heatmap", "Download normalized heatmap")))
  ),
  nav_panel(
    "Primary Results",
    card(card_header(tagList("Primary screen table ", metric_link("Replicate CV"))), DTOutput("primary_table")),
    card(class = "plot-card", card_header("Waterfall plot"), plotlyOutput("waterfall", height = "700px"), div(class = "download-row", downloadButton("download_waterfall", "Download waterfall plot")))
  ),
  nav_panel(
    "Secondary Curves",
    card(card_header(tagList("Dose-response results ", metric_link("EC50 or IC50"), " | ", metric_link("Hill slope"))), DTOutput("dose_table")),
    card(card_header("Detailed curve QC"), DTOutput("dose_qc_table")),
    card(class = "plot-card", card_header("Dose-response plot"), plotlyOutput("dose_plot", height = "560px"), div(class = "download-row", downloadButton("download_dose_plot", "Download dose-response plot"))),
    card(class = "plot-card", card_header("Replicate noise by dose"), plotlyOutput("replicate_cv_plot", height = "520px"), div(class = "download-row", downloadButton("download_replicate_cv_plot", "Download replicate CV plot")))
  ),
  nav_panel(
    "Counter-Assays",
    card(card_header("Counter-assay QC"), DTOutput("counter_table")),
    card(card_header("Artifact-aware hit table"), DTOutput("hit_table"))
  ),
  nav_panel(
    "Plots",
    layout_columns(
      col_widths = c(3, 9),
      card(card_header("Custom plot"), uiOutput("custom_plot_controls"), downloadButton("download_custom_plot", "Download custom plot")),
      card(class = "plot-card", card_header("Selected plot"), plotlyOutput("custom_plot", height = "620px"))
    ),
    card(class = "plot-card", card_header("Control behavior"), plotlyOutput("control_plot", height = "520px"), div(class = "download-row", downloadButton("download_control_plot", "Download control plot"))),
    card(class = "plot-card", card_header("Edge effect review"), plotlyOutput("edge_effect_plot", height = "520px"), div(class = "download-row", downloadButton("download_edge_plot", "Download edge plot")))
  ),
  nav_panel(
    "Final QC",
    card(card_header("Final QC table"), DTOutput("final_qc_table")),
    card(card_header("Export"), div(class = "download-row", downloadButton("download_final_qc", "Final QC CSV"), downloadButton("download_normalized", "Normalized results CSV"), downloadButton("download_package", "Run package ZIP")))
  )
)

server <- function(input, output, session) {
  demo_loaded <- reactiveVal(FALSE)
  analysis_results <- reactiveVal(NULL)
  active_inputs <- reactiveVal(NULL)

  lapply(names(metric_help), function(label) {
    observeEvent(input[[paste0("help_", gsub("[^A-Za-z0-9]", "_", label))]], {
      showModal(modalDialog(title = label, metric_help[[label]], easyClose = TRUE, footer = modalButton("Close")))
    })
  })

  observeEvent(input$load_demo, {
    demo_loaded(TRUE)
    showNotification("Demo data loaded. Click Run analysis.", type = "message")
  })

  metadata_from_form <- reactive({
    data.frame(
      field = c("scientist", "project", "assay_date", "target_id", "cell_line", "protocol_version", "incubation_hours", "readout_nm", "cell_passage", "cell_lot", "quanti_blue_lot", "instrument", "notes"),
      value = c(input$scientist, input$project, as.character(input$assay_date), input$target_id, input$cell_line, input$protocol_version, as.character(input$incubation_hours), as.character(input$readout_nm), input$cell_passage, input$cell_lot, input$reagent_lot, input$instrument, input$run_notes),
      required = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
      stringsAsFactors = FALSE
    )
  })

  get_inputs <- reactive({
    if (!demo_loaded() && is.null(input$raw_file)) return(NULL)
    raw_data <- read_uploaded_csv(input$raw_file, "data/simulated/raw_plate_reader.csv")
    plate_map <- read_uploaded_csv(input$plate_map_file, "data/simulated/plate_map.csv")
    metadata <- if (!is.null(input$metadata_file)) {
      read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE)
    } else if (demo_loaded()) {
      read.csv("data/simulated/run_metadata.csv", stringsAsFactors = FALSE)
    } else {
      metadata_from_form()
    }
    list(raw_data = raw_data, plate_map = plate_map, metadata = metadata)
  })

  observeEvent(input$run_analysis, {
    dat <- get_inputs()
    validate(need(!is.null(dat), "Load demo data or upload raw files first."))
    results <- run_hekblue_analysis(dat$raw_data, dat$plate_map, dat$metadata)
    analysis_results(results)
    active_inputs(dat)
    showNotification("Analysis complete.", type = "message")
  })

  output$schema_table <- renderDT(status_datatable(schema_table, NULL, 10))

  output$run_status <- renderUI({
    res <- analysis_results()
    if (is.null(res)) return(tagList(status_badge("NOT RUN"), p(class = "small-note", paste("Version", APP_VERSION))))
    final <- res$final_qc_table
    tagList(
      layout_columns(
        div(class = "metric-card", div(class = "metric-label", "Peptides reviewed"), div(class = "metric-value", nrow(final))),
        div(class = "metric-card", div(class = "metric-label", "Final actions"), div(class = "small-note", paste(unique(final$final_action), collapse = ", ")))
      )
    )
  })

  output$metadata_card <- renderUI({
    res <- analysis_results()
    if (is.null(res)) return(status_badge("NOT RUN"))
    meta <- res$metadata_completeness
    div(class = "metric-card", div(class = "metric-label", "Metadata completeness"), div(class = "metric-value", paste0(meta$value, "%")), status_badge(meta$status), div(class = "small-note", paste(meta$required_complete, "of", meta$required_total, "required fields complete")))
  })

  output$raw_preview <- renderDT({ dat <- get_inputs(); if (is.null(dat)) return(datatable(data.frame(Message = "Load demo data or upload raw data."))); status_datatable(head(dat$raw_data, 200), NULL, 12) })
  output$plate_map_preview <- renderDT({ dat <- get_inputs(); if (is.null(dat)) return(datatable(data.frame(Message = "Load demo data or upload a plate map."))); status_datatable(head(dat$plate_map, 200), NULL, 12) })
  output$metadata_preview <- renderDT({ dat <- get_inputs(); if (is.null(dat)) return(datatable(data.frame(Message = "Load demo data or enter metadata."))); status_datatable(dat$metadata, NULL, 15) })
  output$metadata_table <- renderDT({ req(active_inputs()); status_datatable(active_inputs()$metadata, NULL, 15) })
  output$eda_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$input_eda, "status", 10) })
  output$cleaned_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$cleaned_well_data, NULL, 20) })
  output$cleaning_summary <- renderDT({
    req(analysis_results())
    z <- as.data.frame(table(analysis_results()$cleaned_well_data$cleaning_action), stringsAsFactors = FALSE)
    names(z) <- c("cleaning_action", "well_count")
    z$status <- ifelse(z$cleaning_action == "KEEP", "PASS", "WARN")
    status_datatable(z, "status", 10)
  })
  output$design_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$design_qc, "design_status", 10) })
  output$plate_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$plate_qc, "plate_qc_status", 10) })
  output$reference_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$reference_control_qc, "status", 15) })
  output$calibration_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$interplate_calibration, "calibration_status", 10) })
  output$primary_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$primary_results, "primary_status", 20) })
  output$dose_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$dose_response_results, NULL, 20) })
  output$dose_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$dose_response_qc, "curve_qc_status", 20) })
  output$counter_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$counter_assay_qc, "counter_qc_status", 10) })
  output$hit_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$hit_calls, NULL, 20) })
  output$final_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$final_qc_table, "final_status", 20) })

  raw_heatmap_obj <- reactive({ req(analysis_results()); plate_heatmap_plot(analysis_results()$cleaned_well_data, "raw_od", "Raw OD by well") })
  norm_heatmap_obj <- reactive({ req(analysis_results()); plate_heatmap_plot(analysis_results()$normalized_results, "percent_activation", "Percent activation by well") })
  qc_plot_obj <- reactive({ req(analysis_results()); qc_bar_plot(analysis_results()$plate_qc) })
  raw_distribution_obj <- reactive({ req(analysis_results()); raw_distribution_plot(analysis_results()$cleaned_well_data) })
  waterfall_obj <- reactive({ req(analysis_results()); waterfall_plot(analysis_results()$primary_results) })
  dose_plot_obj <- reactive({ req(analysis_results()); dose_response_plot(analysis_results()$primary_results, analysis_results()$dose_response_results) })
  control_plot_obj <- reactive({ req(analysis_results()); control_boxplot(analysis_results()$cleaned_well_data) })
  edge_plot_obj <- reactive({ req(analysis_results()); edge_plot(analysis_results()$cleaned_well_data) })
  replicate_cv_obj <- reactive({ req(analysis_results()); replicate_cv_plot(analysis_results()$primary_results) })
  calibration_plot_obj <- reactive({ req(analysis_results()); calibration_plot(analysis_results()$interplate_calibration) })

  output$raw_heatmap <- renderPlotly(plotly::layout(ggplotly(raw_heatmap_obj(), tooltip = c("x", "y", "fill")), margin = list(l = 60, r = 90, t = 60, b = 50)))
  output$normalized_heatmap <- renderPlotly(plotly::layout(ggplotly(norm_heatmap_obj(), tooltip = c("x", "y", "fill")), margin = list(l = 60, r = 90, t = 60, b = 50)))
  output$qc_plot <- renderPlotly(plotly::layout(ggplotly(qc_plot_obj()), margin = list(l = 120, r = 30, t = 60, b = 50)))
  output$raw_distribution <- renderPlotly(ggplotly(raw_distribution_obj()))
  output$waterfall <- renderPlotly(plotly::layout(ggplotly(waterfall_obj()), margin = list(l = 180, r = 30, t = 60, b = 60)))
  output$dose_plot <- renderPlotly(ggplotly(dose_plot_obj()))
  output$control_plot <- renderPlotly(ggplotly(control_plot_obj()))
  output$edge_effect_plot <- renderPlotly(ggplotly(edge_plot_obj()))
  output$replicate_cv_plot <- renderPlotly(ggplotly(replicate_cv_obj()))
  output$calibration_plot <- renderPlotly(ggplotly(calibration_plot_obj()))

  output$custom_plot_controls <- renderUI({
    req(analysis_results())
    cols <- names(analysis_results()$normalized_results)
    tagList(
      selectInput("custom_x", "X", choices = cols, selected = "concentration_uM"),
      selectInput("custom_y", "Y", choices = cols, selected = "percent_activation"),
      selectInput("custom_color", "Color", choices = cols, selected = "peptide_id"),
      selectInput("custom_facet", "Facet", choices = c("None", cols), selected = "assay_mode")
    )
  })

  custom_plot_obj <- reactive({
    req(analysis_results(), input$custom_x, input$custom_y, input$custom_color)
    df <- analysis_results()$normalized_results
    p <- ggplot(df, aes(x = .data[[input$custom_x]], y = .data[[input$custom_y]], color = .data[[input$custom_color]])) + geom_point(alpha = 0.82, size = 2) + hek_theme()
    if (!is.null(input$custom_facet) && input$custom_facet != "None") p <- p + facet_wrap(stats::as.formula(paste("~", input$custom_facet)))
    p
  })
  output$custom_plot <- renderPlotly(ggplotly(custom_plot_obj()))

  output$download_demo_raw <- downloadHandler(filename = function() "hekblue_demo_raw_plate_reader.csv", content = function(file) file.copy("data/simulated/raw_plate_reader.csv", file))
  output$download_demo_plate_map <- downloadHandler(filename = function() "hekblue_demo_plate_map.csv", content = function(file) file.copy("data/simulated/plate_map.csv", file))
  output$download_demo_metadata <- downloadHandler(filename = function() "hekblue_demo_run_metadata.csv", content = function(file) file.copy("data/simulated/run_metadata.csv", file))
  output$download_demo_peptides <- downloadHandler(filename = function() "hekblue_demo_peptide_metadata.csv", content = function(file) file.copy("data/simulated/peptide_metadata.csv", file))
  output$download_demo_target <- downloadHandler(filename = function() "hekblue_demo_target_metadata.csv", content = function(file) file.copy("data/simulated/target_metadata.csv", file))

  output$download_raw_heatmap <- download_plot(raw_heatmap_obj, "raw_heatmap", active_inputs)
  output$download_norm_heatmap <- download_plot(norm_heatmap_obj, "normalized_heatmap", active_inputs)
  output$download_qc_plot <- download_plot(qc_plot_obj, "zprime_qc", active_inputs)
  output$download_raw_distribution <- download_plot(raw_distribution_obj, "raw_od_distribution", active_inputs)
  output$download_waterfall <- download_plot(waterfall_obj, "primary_waterfall", active_inputs)
  output$download_dose_plot <- download_plot(dose_plot_obj, "dose_response", active_inputs)
  output$download_control_plot <- download_plot(control_plot_obj, "control_behavior", active_inputs)
  output$download_edge_plot <- download_plot(edge_plot_obj, "edge_effect", active_inputs)
  output$download_replicate_cv_plot <- download_plot(replicate_cv_obj, "replicate_cv", active_inputs)
  output$download_calibration_plot <- download_plot(calibration_plot_obj, "interplate_calibration", active_inputs)
  output$download_custom_plot <- download_plot(custom_plot_obj, "custom_plot", active_inputs)

  output$download_final_qc <- downloadHandler(filename = function() "final_qc_table.csv", content = function(file) write.csv(analysis_results()$final_qc_table, file, row.names = FALSE))
  output$download_normalized <- downloadHandler(filename = function() "normalized_results.csv", content = function(file) write.csv(analysis_results()$normalized_results, file, row.names = FALSE))
  output$download_package <- downloadHandler(
    filename = function() paste0("hekblue_run_package_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = function(file) {
      dat <- active_inputs()
      tmp <- tempfile("hekblue_export_")
      export_hekblue_results(analysis_results(), dat$raw_data, dat$plate_map, dat$metadata, tmp)
      old <- setwd(dirname(tmp))
      on.exit(setwd(old), add = TRUE)
      export_files <- list.files(basename(tmp), recursive = TRUE, full.names = TRUE)
      utils::zip(file, files = export_files)
    }
  )
}

shinyApp(ui, server)
