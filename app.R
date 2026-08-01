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

read_uploaded_csv <- function(file_input, demo_path) {
  if (is.null(file_input)) {
    read.csv(demo_path, stringsAsFactors = FALSE)
  } else {
    read.csv(file_input$datapath, stringsAsFactors = FALSE)
  }
}

status_badge <- function(status) {
  color <- switch(status, PASS = "success", WARN = "warning", FAIL = "danger", "secondary")
  span(class = paste("badge rounded-pill text-bg", color, sep = "-"), status)
}

ui <- page_navbar(
  title = "HEKBlueR",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  header = tags$head(tags$style(HTML("
    .metric-card {border: 1px solid #d7dee8; border-radius: 8px; padding: 14px; background: #ffffff; min-height: 104px;}
    .metric-label {font-size: 0.86rem; color: #516173; margin-bottom: 4px;}
    .metric-value {font-size: 1.4rem; font-weight: 700; color: #172033;}
    .small-note {font-size: 0.86rem; color: #5d6b7a;}
    .table-title {font-weight: 700; margin: 14px 0 8px;}
  "))),
  nav_panel(
    "Start",
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Workflow"),
        tags$ol(
          tags$li("Upload raw HEK-Blue plate-reader data."),
          tags$li("Upload a plate map and run metadata."),
          tags$li("Choose agonist, antagonist, or counter-assay interpretation."),
          tags$li("Run automated design QC, plate QC, normalization, and curve review."),
          tags$li("Export database-ready tables and a reproducibility package.")
        ),
        actionButton("load_demo", "Load demo data", class = "btn-primary"),
        span(" "),
        actionButton("run_analysis", "Run analysis", class = "btn-success")
      ),
      card(
        card_header("Run status"),
        uiOutput("run_status"),
        div(class = "small-note", "Demo data includes one target, three peptides, primary agonist, antagonist, and counter-assay plates.")
      )
    )
  ),
  nav_panel(
    "Upload",
    layout_columns(
      col_widths = c(4, 4, 4),
      card(card_header("Raw data"), fileInput("raw_file", "Raw plate-reader CSV"), downloadButton("download_demo_raw", "Demo raw CSV")),
      card(card_header("Plate map"), fileInput("plate_map_file", "Plate map CSV"), downloadButton("download_demo_plate_map", "Demo plate map")),
      card(card_header("Metadata"), fileInput("metadata_file", "Run metadata CSV"), downloadButton("download_demo_metadata", "Demo metadata"))
    ),
    card(
      card_header("Metadata helper"),
      layout_columns(
        textInput("scientist", "Scientist name", ""),
        textInput("project", "Project", ""),
        dateInput("assay_date", "Assay date", value = Sys.Date()),
        textInput("target_id", "Target ID", "TARGET_TLR8_DEMO"),
        textInput("cell_line", "Cell line", ""),
        textInput("cell_passage", "Cell passage", ""),
        textInput("cell_lot", "Cell lot", ""),
        textInput("reagent_lot", "QUANTI-Blue lot", ""),
        textInput("instrument", "Instrument", ""),
        textInput("protocol_version", "Protocol version", "v1.0"),
        numericInput("incubation_hours", "Incubation hours", value = 18, min = 0),
        numericInput("readout_nm", "Readout wavelength", value = 655, min = 400, max = 800)
      ),
      textAreaInput("run_notes", "Notes and protocol deviations", "", height = "90px")
    )
  ),
  nav_panel(
    "Metadata",
    layout_columns(
      col_widths = c(4, 8),
      card(card_header("Completeness"), uiOutput("metadata_card")),
      card(card_header("Submitted metadata"), DTOutput("metadata_table"))
    )
  ),
  nav_panel(
    "Design QC",
    card(card_header("Experimental design review"), DTOutput("design_qc_table"))
  ),
  nav_panel(
    "Plate QC",
    layout_columns(
      col_widths = c(12),
      card(card_header("QC metrics"), DTOutput("plate_qc_table")),
      card(card_header("Z-prime summary"), plotlyOutput("qc_plot", height = "360px"))
    )
  ),
  nav_panel(
    "Plate Layout",
    card(card_header("Raw OD heatmap"), plotlyOutput("raw_heatmap", height = "520px")),
    card(card_header("Normalized heatmap"), plotlyOutput("normalized_heatmap", height = "520px"))
  ),
  nav_panel(
    "Primary Results",
    card(card_header("Primary screen table"), DTOutput("primary_table")),
    card(card_header("Waterfall plot"), plotlyOutput("waterfall", height = "520px"))
  ),
  nav_panel(
    "Secondary Curves",
    card(card_header("Dose-response results"), DTOutput("dose_table")),
    card(card_header("Dose-response plot"), plotlyOutput("dose_plot", height = "520px"))
  ),
  nav_panel(
    "Counter-Assays",
    card(card_header("Counter-assay QC"), DTOutput("counter_table")),
    card(card_header("Artifact interpretation"), DTOutput("hit_table"))
  ),
  nav_panel(
    "Plots",
    layout_columns(
      col_widths = c(3, 9),
      card(
        card_header("Custom plot"),
        uiOutput("custom_plot_controls")
      ),
      card(
        card_header("Selected plot"),
        plotlyOutput("custom_plot", height = "560px")
      )
    ),
    card(card_header("Control behavior"), plotlyOutput("control_plot", height = "440px")),
    card(card_header("Edge effect review"), plotlyOutput("edge_effect_plot", height = "440px"))
  ),
  nav_panel(
    "Final QC",
    card(card_header("Final QC table"), DTOutput("final_qc_table")),
    card(
      card_header("Export"),
      downloadButton("download_final_qc", "Final QC CSV"),
      span(" "),
      downloadButton("download_normalized", "Normalized results CSV"),
      span(" "),
      downloadButton("download_package", "Run package ZIP")
    )
  )
)

server <- function(input, output, session) {
  demo_loaded <- reactiveVal(FALSE)
  analysis_results <- reactiveVal(NULL)
  active_inputs <- reactiveVal(NULL)

  observeEvent(input$load_demo, {
    demo_loaded(TRUE)
    showNotification("Demo data loaded. Click Run analysis.", type = "message")
  })

  metadata_from_form <- reactive({
    data.frame(
      field = c("scientist", "project", "assay_date", "target_id", "cell_line", "cell_passage", "cell_lot", "quanti_blue_lot", "instrument", "protocol_version", "incubation_hours", "readout_nm", "notes"),
      value = c(input$scientist, input$project, as.character(input$assay_date), input$target_id, input$cell_line, input$cell_passage, input$cell_lot, input$reagent_lot, input$instrument, input$protocol_version, as.character(input$incubation_hours), as.character(input$readout_nm), input$run_notes),
      required = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
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

  output$run_status <- renderUI({
    res <- analysis_results()
    if (is.null(res)) {
      return(tagList(status_badge("NOT RUN"), tags$p("Load demo data or upload files to begin.")))
    }
    final <- res$final_qc_table
    tagList(
      div(class = "metric-card",
          div(class = "metric-label", "Peptides reviewed"),
          div(class = "metric-value", nrow(final)),
          div(class = "small-note", paste("Actions:", paste(unique(final$final_action), collapse = ", ")))
      )
    )
  })

  output$metadata_card <- renderUI({
    res <- analysis_results()
    if (is.null(res)) return(status_badge("NOT RUN"))
    meta <- res$metadata_completeness
    tagList(
      div(class = "metric-card",
          div(class = "metric-label", "Metadata completeness"),
          div(class = "metric-value", paste0(meta$value, "%")),
          status_badge(meta$status)
      )
    )
  })

  output$metadata_table <- renderDT({
    dat <- active_inputs()
    if (is.null(dat)) return(datatable(data.frame()))
    datatable(dat$metadata, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$design_qc_table <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$design_qc, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$plate_qc_table <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$plate_qc, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$primary_table <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$primary_results, options = list(pageLength = 20, scrollX = TRUE))
  })

  output$dose_table <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$dose_response_results, options = list(pageLength = 20, scrollX = TRUE))
  })

  output$counter_table <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$counter_assay_qc, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$hit_table <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$hit_calls, options = list(pageLength = 20, scrollX = TRUE))
  })

  output$final_qc_table <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$final_qc_table, options = list(pageLength = 20, scrollX = TRUE))
  })

  output$qc_plot <- renderPlotly({
    req(analysis_results())
    ggplotly(qc_bar_plot(analysis_results()$plate_qc))
  })

  output$raw_heatmap <- renderPlotly({
    req(analysis_results())
    ggplotly(plate_heatmap_plot(analysis_results()$cleaned_well_data, "raw_od", "Raw OD by well"))
  })

  output$normalized_heatmap <- renderPlotly({
    req(analysis_results())
    ggplotly(plate_heatmap_plot(analysis_results()$normalized_results, "percent_activation", "Percent activation by well"))
  })

  output$waterfall <- renderPlotly({
    req(analysis_results())
    ggplotly(waterfall_plot(analysis_results()$primary_results))
  })

  output$dose_plot <- renderPlotly({
    req(analysis_results())
    ggplotly(dose_response_plot(analysis_results()$primary_results, analysis_results()$dose_response_results))
  })

  output$control_plot <- renderPlotly({
    req(analysis_results())
    ggplotly(control_boxplot(analysis_results()$cleaned_well_data))
  })

  output$edge_effect_plot <- renderPlotly({
    req(analysis_results())
    ggplotly(edge_plot(analysis_results()$cleaned_well_data))
  })

  output$custom_plot_controls <- renderUI({
    req(analysis_results())
    df <- analysis_results()$normalized_results
    cols <- names(df)
    tagList(
      selectInput("custom_x", "X", choices = cols, selected = "concentration_uM"),
      selectInput("custom_y", "Y", choices = cols, selected = "percent_activation"),
      selectInput("custom_color", "Color", choices = cols, selected = "peptide_id"),
      selectInput("custom_facet", "Facet", choices = c("None", cols), selected = "assay_mode")
    )
  })

  output$custom_plot <- renderPlotly({
    req(analysis_results(), input$custom_x, input$custom_y, input$custom_color)
    df <- analysis_results()$normalized_results
    p <- ggplot(df, aes(x = .data[[input$custom_x]], y = .data[[input$custom_y]], color = .data[[input$custom_color]])) +
      geom_point(alpha = 0.8) +
      theme_minimal(base_size = 12)
    if (!is.null(input$custom_facet) && input$custom_facet != "None") {
      p <- p + facet_wrap(stats::as.formula(paste("~", input$custom_facet)))
    }
    ggplotly(p)
  })

  output$download_demo_raw <- downloadHandler(
    filename = function() "hekblue_demo_raw_plate_reader.csv",
    content = function(file) file.copy("data/simulated/raw_plate_reader.csv", file)
  )
  output$download_demo_plate_map <- downloadHandler(
    filename = function() "hekblue_demo_plate_map.csv",
    content = function(file) file.copy("data/simulated/plate_map.csv", file)
  )
  output$download_demo_metadata <- downloadHandler(
    filename = function() "hekblue_demo_run_metadata.csv",
    content = function(file) file.copy("data/simulated/run_metadata.csv", file)
  )
  output$download_final_qc <- downloadHandler(
    filename = function() "final_qc_table.csv",
    content = function(file) write.csv(analysis_results()$final_qc_table, file, row.names = FALSE)
  )
  output$download_normalized <- downloadHandler(
    filename = function() "normalized_results.csv",
    content = function(file) write.csv(analysis_results()$normalized_results, file, row.names = FALSE)
  )
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
