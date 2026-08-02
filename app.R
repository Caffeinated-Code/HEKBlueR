suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(plotly)
})

source("R/qc_metrics.R")
source("R/qc_thresholds.R")
source("R/export.R")
source("R/analysis.R")
source("R/plots.R")
source("R/simulate_data.R")

APP_VERSION <- "0.4.0"

metric_help <- list(
  "Z-prime" = "Z-prime measures separation between positive and negative controls. Values above 0.5 are preferred. Values from 0.3 to 0.5 need review. Lower values usually mean the assay window is weak.",
  "Robust Z-prime" = "Robust Z-prime uses medians and median absolute deviation. It is less sensitive to individual outlier wells.",
  "SSMD" = "SSMD is a standardized difference between control groups. Larger absolute values indicate stronger separation.",
  "Signal-to-background" = "Signal-to-background compares positive control signal against baseline. Low values mean the assay may have too little dynamic range for reliable screening.",
  "Signal window" = "Signal window is positive control minus negative control. It is the usable assay range after blank correction.",
  "Control CV" = "Control coefficient of variation measures control stability. Values above 20 percent often indicate pipetting, reagent, incubation, or reader variability.",
  "Replicate CV" = "Replicate CV measures technical replicate noise at a peptide dose. High CV can make hit calls and curve fits unreliable.",
  "Edge effect" = "Edge effect compares edge wells to center wells. A large difference can reflect evaporation, temperature gradients, or plate handling effects.",
  "Row and column bias" = "Row and column bias checks spatial drift across the plate. This can indicate dispense order, incubation gradient, or reader position effects.",
  "Inter-plate calibration" = "Inter-plate calibration aligns plates using a shared calibrator or shared positive control. Drift above 0.15 OD is flagged for review.",
  "Reference control stability" = "Reference control stability checks whether shared controls are consistent enough to support normalization and plate comparison.",
  "EC50 or IC50" = "EC50 is the concentration giving half-maximal activation. IC50 is the concentration giving half-maximal inhibition. Estimates are strongest when they sit inside the tested dose range.",
  "Hill slope" = "Hill slope describes curve steepness. Very shallow or very steep slopes can indicate weak biology, noisy data, or fitting instability.",
  "Dynamic range" = "Dynamic range is the distance between the lowest and highest observed response. Weak dynamic range limits confidence in potency estimates.",
  "RMSE" = "RMSE summarizes residual error between observed responses and the fitted curve. Lower values indicate a curve that follows the observed data more closely.",
  "Plateau checks" = "Plateau checks ask whether the tested dose range captures the low and high ends of the curve. Missing plateaus make potency estimates more extrapolated.",
  "Monotonicity" = "Monotonicity checks whether dose increases generally move response in the expected direction. Strong reversals can indicate noise, toxicity, solubility limits, or mixed mechanisms.",
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
    filter = "top",
    rownames = FALSE,
    options = list(pageLength = page_length, scrollX = TRUE, scrollY = "48vh", scrollCollapse = TRUE)
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

analysis_signature <- function(dat) {
  tmp <- tempfile(fileext = ".rds")
  saveRDS(dat, tmp, version = 2)
  on.exit(unlink(tmp), add = TRUE)
  unname(tools::md5sum(tmp))
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

dose_glossary <- data.frame(
  QC_metric = c("n_dose_points", "dynamic_range", "rmse", "max_residual", "max_replicate_cv", "monotonic_violations", "ec50_in_range", "top_plateau_observed", "bottom_plateau_observed", "curve_flags"),
  What_it_means = c(
    "Number of unique concentrations used for a peptide curve.",
    "Observed response range across the tested concentration series.",
    "Average mismatch between observed data and fitted curve.",
    "Largest observed-vs-fitted curve mismatch.",
    "Highest technical replicate CV across dose points.",
    "Number of dose steps that move opposite the expected trend.",
    "Whether EC50 or IC50 falls inside the tested dose range.",
    "Whether the high-response plateau appears covered by the data.",
    "Whether the low-response plateau appears covered by the data.",
    "Compact list of curve warnings."
  ),
  How_to_interpret = c(
    "Eight or more points is preferred for robust secondary screening.",
    "Low dynamic range means potency estimates may be unstable.",
    "Lower RMSE is better. High RMSE means the curve does not explain the data well.",
    "A high max residual means one dose may be driving the fit.",
    "Values above 20 to 30 percent indicate noisy technical replicates.",
    "Multiple violations can indicate noisy data, toxicity, precipitation, or mixed biology.",
    "TRUE is preferred. FALSE means the potency estimate is extrapolated.",
    "TRUE is preferred for confident efficacy and potency estimates.",
    "TRUE is preferred for confident baseline and potency estimates.",
    "GOOD_CURVE is strongest. Other flags require review before advancing."
  ),
  stringsAsFactors = FALSE
)

threshold_glossary <- qc_threshold_table()

ui <- page_navbar(
  title = "HEKBlueR",
  theme = bs_theme(
    version = 5,
    primary = "#0f766e",
    secondary = "#5b6d68",
    success = "#0f766e",
    warning = "#9a5f22",
    danger = "#b42318",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter")
  ),
  header = tags$head(tags$style(HTML("
    :root { --bg: #f8fbfa; --bg-soft: #edf6f2; --surface: #ffffff; --surface-strong: #e4f0ec; --text: #16211f; --muted: #5b6d68; --line: #d5e4de; --brand: #0f766e; --brand-dark: #115e59; --accent: #9a5f22; --danger: #b42318; --warn-bg: #fff4dc; --pass-bg: #e1f4ee; --fail-bg: #fde2df; --shadow: 0 18px 52px rgba(22,33,31,0.10); }
    body { background: var(--bg); color: var(--text); line-height: 1.55; }
    .navbar { background: color-mix(in srgb, var(--bg) 92%, white); border-bottom: 1px solid var(--line); box-shadow: 0 8px 30px rgba(22,33,31,0.08); }
    .navbar-brand { font-weight: 850; letter-spacing: 0; color: var(--brand-dark) !important; }
    .nav-link { font-weight: 720; color: var(--muted) !important; }
    .nav-link.active { color: var(--brand-dark) !important; }
    .bslib-card { border: 1px solid var(--line); border-radius: 8px; box-shadow: var(--shadow); background: var(--surface); margin-bottom: 1rem; }
    .card-header { background: var(--surface); color: var(--text); font-weight: 800; border-bottom: 1px solid var(--line); }
    .hero-band { padding: 1.3rem; border: 1px solid var(--line); border-radius: 8px; background: radial-gradient(circle at top left, color-mix(in srgb, var(--brand) 13%, transparent), transparent 30rem), linear-gradient(180deg, var(--bg-soft), var(--surface)); }
    .eyebrow { color: var(--brand-dark); font-size: 0.76rem; font-weight: 850; letter-spacing: 0.08em; text-transform: uppercase; margin-bottom: 0.35rem; }
    .metric-card { border: 1px solid var(--line); border-radius: 8px; padding: 14px; background: var(--surface); min-height: 112px; }
    .metric-label { font-size: 0.82rem; color: var(--muted); margin-bottom: 4px; font-weight: 680; }
    .metric-value { font-size: 1.45rem; font-weight: 820; color: var(--text); }
    .small-note { font-size: 0.88rem; color: var(--muted); }
    .section-note { padding: 0.8rem 1rem; border-left: 4px solid var(--brand); background: var(--bg-soft); border-radius: 6px; color: var(--muted); margin-bottom: 1rem; }
    .status-badge { display: inline-block; border-radius: 999px; padding: 0.22rem 0.62rem; font-weight: 800; font-size: 0.76rem; }
    .status-badge.ok { background: var(--pass-bg); color: #143f2d; }
    .status-badge.warn { background: var(--warn-bg); color: #5e3a05; }
    .status-badge.fail { background: var(--fail-bg); color: #651b13; }
    .status-badge.neutral { background: #e7efeb; color: #263934; }
    .field-badge { margin-left: 8px; border-radius: 999px; padding: 2px 8px; font-size: 0.70rem; font-weight: 820; }
    .field-badge.required { background: var(--pass-bg); color: #143f2d; }
    .field-badge.optional { background: #e7efeb; color: #314842; }
    .metric-help { font-weight: 800; color: var(--brand-dark); text-decoration: none; border-bottom: 1px dotted var(--brand); }
    .global-search-wrap { min-width: min(28rem, 82vw); padding: 0.35rem 0.6rem; }
    .global-search-wrap input { border: 1px solid var(--line); border-radius: 999px; padding: 0.45rem 0.9rem; font-size: 0.9rem; }
    .page-search-hit { outline: 2px solid color-mix(in srgb, var(--brand) 45%, transparent); outline-offset: 2px; }
    .page-search-dim { opacity: 0.38; }
    .download-row { margin-top: 8px; display: flex; gap: 8px; flex-wrap: wrap; }
    .btn { border-radius: 999px; font-weight: 760; }
    .btn-default, .btn-secondary { border-color: var(--line); background: var(--surface); color: var(--text); }
    .plot-card { overflow: visible; }
    .plot-card .card-body { min-height: 360px; }
    .plot-card .html-widget, .plot-card .shiny-plot-output, .plot-card .js-plotly-plot { width: 100% !important; min-height: 460px !important; height: 460px !important; }
    .plot-card.plot-tall .card-body { min-height: 720px; }
    .plot-card.plot-tall .html-widget, .plot-card.plot-tall .shiny-plot-output, .plot-card.plot-tall .js-plotly-plot { min-height: 700px !important; height: 700px !important; }
    .plot-card.plot-wide .card-body { min-height: 560px; }
    .plot-card.plot-wide .html-widget, .plot-card.plot-wide .shiny-plot-output, .plot-card.plot-wide .js-plotly-plot { min-height: 540px !important; height: 540px !important; }
    .upload-compact .bslib-card { box-shadow: none; }
    .upload-compact .card-body { padding: 0.9rem; }
    .upload-compact .form-group { margin-bottom: 0.65rem; }
    .dataTables_scrollBody { border-bottom: 1px solid var(--line); }
    .dataTables_scrollHead th { position: sticky; top: 0; z-index: 3; }
    .form-section-title { font-size: 0.78rem; font-weight: 850; letter-spacing: 0.08em; text-transform: uppercase; color: var(--accent); margin: 0.5rem 0 0.75rem; }
    .tab-content { padding-top: 0.75rem; }
    table.dataTable thead th { background: var(--surface-strong); color: var(--text); }
    @media (max-width: 900px) {
      .container-fluid { padding-left: 10px; padding-right: 10px; }
      .bslib-grid { grid-template-columns: 1fr !important; }
      .metric-card { min-height: auto; }
      .navbar-nav { gap: 0.15rem; }
      .global-search-wrap { width: 100%; padding-left: 0; }
    }
  ")),
  tags$script(HTML("
    $(document).on('input', '#page_search', function() {
      var q = String(this.value || '').toLowerCase().trim();
      var activePane = $('.tab-pane.active');
      var searchable = activePane.find('.bslib-card, .metric-card, .section-note, .tab-content');
      searchable.removeClass('page-search-hit page-search-dim');
      if (!q) return;
      searchable.each(function() {
        var hit = $(this).text().toLowerCase().indexOf(q) >= 0;
        $(this).toggleClass('page-search-hit', hit);
        $(this).toggleClass('page-search-dim', !hit);
      });
    });
    $(document).on('shown.bs.tab shown.bs.collapse', function() {
      setTimeout(function() {
        window.dispatchEvent(new Event('resize'));
        if (window.Plotly) {
          $('.tab-pane.active .js-plotly-plot:visible').each(function() {
            Plotly.Plots.resize(this);
          });
        }
      }, 200);
    });
  "))),
  nav_item(tags$div(class = "global-search-wrap", tags$input(id = "page_search", type = "search", class = "form-control", placeholder = "Search current page..."))),
  nav_panel(
    "Demo",
    layout_columns(
      col_widths = c(7, 5),
      card(
        class = "hero-band",
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
    div(class = "section-note", "Upload raw OD data and plate maps first. Use the metadata form when a metadata CSV is not available. Required fields support reproducibility. Optional fields improve search, audit, and troubleshooting later."),
    tabsetPanel(
      tabPanel(
        "Files",
        layout_columns(
          class = "upload-compact",
          col_widths = c(4, 4, 4),
          card(card_header(required_header("Raw plate-reader CSV", TRUE)), p(class = "small-note", "One row per well. Required columns are listed above."), fileInput("raw_file", NULL, buttonLabel = "Browse", placeholder = "No file selected")),
          card(card_header(required_header("Plate map CSV", TRUE)), p(class = "small-note", "Well role, sample, peptide, concentration, and replicate structure."), fileInput("plate_map_file", NULL, buttonLabel = "Browse", placeholder = "No file selected")),
          card(card_header(required_header("Run metadata CSV", FALSE)), p(class = "small-note", "Optional if using the form. Strongly recommended for reproducibility."), fileInput("metadata_file", NULL, buttonLabel = "Browse", placeholder = "No file selected"))
        )
      ),
      tabPanel(
        "Required Metadata",
        layout_columns(
          class = "upload-compact",
          col_widths = c(4, 4, 4),
          textInput("scientist", "Scientist name", ""),
          textInput("project", "Project", ""),
          dateInput("assay_date", "Assay date", value = Sys.Date()),
          textInput("target_id", "Target ID", "TARGET_TLR8_DEMO"),
          textInput("cell_line", "Cell line", ""),
          textInput("protocol_version", "Protocol version", "v1.0"),
          numericInput("incubation_hours", "Incubation hours", value = 18, min = 0),
          numericInput("readout_nm", "Readout wavelength", value = 655, min = 400, max = 800)
        )
      ),
      tabPanel(
        "Optional Metadata",
        layout_columns(
          class = "upload-compact",
          col_widths = c(4, 4, 4),
          textInput("cell_passage", "Cell passage", ""),
          textInput("cell_lot", "Cell lot", ""),
          textInput("reagent_lot", "QUANTI-Blue lot", ""),
          textInput("instrument", "Instrument", ""),
          textInput("peptide_lot", "Peptide lot", ""),
          textInput("peptide_purity", "Peptide purity", ""),
          textInput("vehicle", "Vehicle", ""),
          textInput("reader_settings", "Plate reader settings", "")
        )
      ),
      tabPanel(
        "Notes",
        card(textAreaInput("run_notes", "Notes and protocol deviations (Optional)", "", height = "150px"))
      ),
      tabPanel(
        "Expected Schema",
        card(
          card_header("Expected raw data schema"),
          DTOutput("schema_table")
        )
      )
    )
  ),
  nav_panel(
    "Uploaded Preview",
    div(class = "section-note", "Use these tabs to inspect exactly what HEKBlueR parsed before analysis. This keeps large tables navigable."),
    tabsetPanel(
      tabPanel("Raw data", card(card_header("Raw data preview"), DTOutput("raw_preview"))),
      tabPanel("Plate map", card(card_header("Plate map preview"), DTOutput("plate_map_preview"))),
      tabPanel("Metadata", card(card_header("Metadata preview"), DTOutput("metadata_preview"))),
      tabPanel("Expected schema", card(card_header("Required and recommended fields"), DTOutput("schema_table_preview")))
    )
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
    div(class = "section-note", "Review cleaning actions first, then inspect well-level rows only when a flag needs follow-up."),
    tabsetPanel(
      tabPanel("Cleaning summary", card(card_header("Cleaning action counts"), DTOutput("cleaning_summary"))),
      tabPanel("Well-level review", card(card_header("Cleaned well data"), DTOutput("cleaned_table")))
    )
  ),
  nav_panel(
    "Metadata",
    layout_columns(
      col_widths = c(4, 4, 4),
      card(card_header(tagList("Overall ", metric_link("Metadata completeness"))), uiOutput("metadata_card")),
      card(card_header("Required score"), uiOutput("required_metadata_card")),
      card(card_header("Optional score"), uiOutput("optional_metadata_card"))
    ),
    tabsetPanel(
      tabPanel("All metadata", card(card_header("Submitted metadata"), DTOutput("metadata_table"))),
      tabPanel("Required fields", card(card_header("Required metadata"), DTOutput("required_metadata_table"))),
      tabPanel("Optional fields", card(card_header("Optional metadata"), DTOutput("optional_metadata_table")))
    )
  ),
  nav_panel("Design QC", card(card_header("Experimental design review"), DTOutput("design_qc_table"))),
  nav_panel(
    "QC Thresholds",
    div(class = "section-note", "These thresholds drive the automated QC status columns. Exported result packages include the same threshold table for traceability."),
    card(card_header("Active QC thresholds"), DTOutput("threshold_table"))
  ),
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
    tabsetPanel(
      tabPanel("Raw OD", card(class = "plot-card plot-tall", card_header("Raw OD heatmap"), plotlyOutput("raw_heatmap", height = "700px"), div(class = "download-row", downloadButton("download_raw_heatmap", "Download raw heatmap")))),
      tabPanel("Normalized", card(class = "plot-card plot-tall", card_header("Normalized heatmap"), plotlyOutput("normalized_heatmap", height = "700px"), div(class = "download-row", downloadButton("download_norm_heatmap", "Download normalized heatmap"))))
    )
  ),
  nav_panel(
    "Primary Results",
    div(class = "section-note", "Primary results are split by analysis task so reviewers can move from tables to plots without scrolling through one long page."),
    tabsetPanel(
      tabPanel("All results", card(card_header(tagList("Primary screen table ", metric_link("Replicate CV"))), DTOutput("primary_table"))),
      tabPanel("Agonist", card(card_header("Agonist results"), DTOutput("primary_agonist_table"))),
      tabPanel("Antagonist", card(card_header("Antagonist results"), DTOutput("primary_antagonist_table"))),
      tabPanel("Sample QC", card(card_header("Sample-level QC summary"), DTOutput("sample_qc_table"))),
      tabPanel("Waterfall", card(class = "plot-card plot-tall", card_header("Interactive waterfall plot"), plotlyOutput("waterfall", height = "700px"), div(class = "download-row", downloadButton("download_waterfall", "Download waterfall plot")))),
      tabPanel("Replicate noise", card(class = "plot-card plot-wide", card_header("Replicate CV by dose"), plotlyOutput("primary_replicate_cv_plot", height = "540px"), div(class = "download-row", downloadButton("download_replicate_cv_plot_primary", "Download replicate CV plot"))))
    )
  ),
  nav_panel(
    "Secondary Curves",
    div(class = "section-note", "Dose-response is the main review workspace. The fitted curve plot is interactive. QC tables cover potency range, residuals, plateaus, monotonicity, and replicate noise."),
    tabsetPanel(
      tabPanel("Interactive curves", card(class = "plot-card plot-tall", card_header(tagList("Dose-response plot ", metric_link("EC50 or IC50"), " | ", metric_link("Hill slope"))), plotlyOutput("dose_plot", height = "700px"), div(class = "download-row", downloadButton("download_dose_plot", "Download dose-response plot")))),
      tabPanel("Curve fit table", card(card_header(tagList("Dose-response results ", metric_link("Dynamic range"), " | ", metric_link("RMSE"))), DTOutput("dose_table"))),
      tabPanel("Curve QC table", card(card_header(tagList("Detailed curve QC ", metric_link("Plateau checks"), " | ", metric_link("Monotonicity"))), DTOutput("dose_qc_table"))),
      tabPanel("Replicate noise", card(class = "plot-card plot-wide", card_header(tagList("Replicate noise by dose ", metric_link("Replicate CV"))), plotlyOutput("replicate_cv_plot", height = "540px"), div(class = "download-row", downloadButton("download_replicate_cv_plot", "Download replicate CV plot")))),
      tabPanel("QC glossary", card(card_header("Dose-response QC explanations"), DTOutput("dose_glossary_table")))
    )
  ),
  nav_panel(
    "Counter-Assays",
    card(card_header("Counter-assay QC"), DTOutput("counter_table")),
    card(card_header("Artifact-aware hit table"), DTOutput("hit_table"))
  ),
  nav_panel(
    "Plots",
    tabsetPanel(
      tabPanel(
        "Custom",
        layout_columns(
          col_widths = c(3, 9),
          card(card_header("Custom plot"), uiOutput("custom_plot_controls"), downloadButton("download_custom_plot", "Download custom plot")),
          card(class = "plot-card plot-wide", card_header("Selected plot"), plotlyOutput("custom_plot", height = "540px"))
        )
      ),
      tabPanel("Controls", card(class = "plot-card plot-wide", card_header("Control behavior"), plotlyOutput("control_plot", height = "540px"), div(class = "download-row", downloadButton("download_control_plot", "Download control plot")))),
      tabPanel("Edge Effects", card(class = "plot-card plot-wide", card_header("Edge effect review"), plotlyOutput("edge_effect_plot", height = "540px"), div(class = "download-row", downloadButton("download_edge_plot", "Download edge plot"))))
    )
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
  analysis_cache <- reactiveValues(signature = NULL, results = NULL, inputs = NULL, message = "No analysis has been run in this session.")

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
      field = c("scientist", "project", "assay_date", "target_id", "cell_line", "protocol_version", "incubation_hours", "readout_nm", "cell_passage", "cell_lot", "quanti_blue_lot", "instrument", "peptide_lot", "peptide_purity", "vehicle", "reader_settings", "notes"),
      value = c(input$scientist, input$project, as.character(input$assay_date), input$target_id, input$cell_line, input$protocol_version, as.character(input$incubation_hours), as.character(input$readout_nm), input$cell_passage, input$cell_lot, input$reagent_lot, input$instrument, input$peptide_lot, input$peptide_purity, input$vehicle, input$reader_settings, input$run_notes),
      required = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
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
    sig <- analysis_signature(dat)
    if (!is.null(analysis_cache$signature) && identical(sig, analysis_cache$signature)) {
      results <- analysis_cache$results
      analysis_cache$message <- "Inputs unchanged. Cached analysis was reused."
      showNotification("Inputs unchanged. Cached analysis reused.", type = "message")
    } else {
      results <- run_hekblue_analysis(dat$raw_data, dat$plate_map, dat$metadata, thresholds = default_qc_thresholds())
      analysis_cache$signature <- sig
      analysis_cache$results <- results
      analysis_cache$inputs <- dat
      analysis_cache$message <- "Inputs changed. Analysis was run from the current data."
      showNotification("Analysis complete.", type = "message")
    }
    analysis_results(results)
    active_inputs(dat)
  })

  output$schema_table <- renderDT(status_datatable(schema_table, NULL, 10))
  output$schema_table_preview <- renderDT(status_datatable(schema_table, NULL, 10))

  output$run_status <- renderUI({
    res <- analysis_results()
    if (is.null(res)) return(tagList(status_badge("NOT RUN"), p(class = "small-note", paste("Version", APP_VERSION))))
    final <- res$final_qc_table
    tagList(
      layout_columns(
        div(class = "metric-card", div(class = "metric-label", "Peptides reviewed"), div(class = "metric-value", nrow(final))),
        div(class = "metric-card", div(class = "metric-label", "Final actions"), div(class = "small-note", paste(unique(final$final_action), collapse = ", "))),
        div(class = "metric-card", div(class = "metric-label", "Analysis cache"), div(class = "small-note", analysis_cache$message))
      )
    )
  })

  output$metadata_card <- renderUI({
    res <- analysis_results()
    if (is.null(res)) return(status_badge("NOT RUN"))
    meta <- res$metadata_completeness
    div(class = "metric-card", div(class = "metric-label", "Metadata completeness"), div(class = "metric-value", paste0(meta$value, "%")), status_badge(meta$status), div(class = "small-note", paste(meta$required_complete, "of", meta$required_total, "required fields complete")))
  })

  output$required_metadata_card <- renderUI({
    res <- analysis_results()
    if (is.null(res)) return(status_badge("NOT RUN"))
    meta <- res$metadata_completeness
    div(class = "metric-card", div(class = "metric-label", "Required fields"), div(class = "metric-value", paste0(meta$required_percent, "%")), status_badge(ifelse(meta$required_complete == meta$required_total, "PASS", "FAIL")), div(class = "small-note", paste(meta$required_complete, "of", meta$required_total, "complete")))
  })

  output$optional_metadata_card <- renderUI({
    res <- analysis_results()
    if (is.null(res)) return(status_badge("NOT RUN"))
    meta <- res$metadata_completeness
    div(class = "metric-card", div(class = "metric-label", "Optional fields"), div(class = "metric-value", paste0(meta$optional_percent, "%")), status_badge(ifelse(meta$optional_percent >= 80, "PASS", ifelse(meta$optional_percent >= 40, "WARN", "FAIL"))), div(class = "small-note", paste(meta$optional_complete, "of", meta$optional_total, "complete")))
  })

  output$raw_preview <- renderDT({ dat <- get_inputs(); if (is.null(dat)) return(datatable(data.frame(Message = "Load demo data or upload raw data."))); status_datatable(head(dat$raw_data, 200), NULL, 12) })
  output$plate_map_preview <- renderDT({ dat <- get_inputs(); if (is.null(dat)) return(datatable(data.frame(Message = "Load demo data or upload a plate map."))); status_datatable(head(dat$plate_map, 200), NULL, 12) })
  output$metadata_preview <- renderDT({ dat <- get_inputs(); if (is.null(dat)) return(datatable(data.frame(Message = "Load demo data or enter metadata."))); status_datatable(dat$metadata, NULL, 15) })
  output$metadata_table <- renderDT({ req(active_inputs()); status_datatable(active_inputs()$metadata, NULL, 15) })
  output$required_metadata_table <- renderDT({
    req(active_inputs())
    md <- active_inputs()$metadata
    status_datatable(md[md$required %in% c(TRUE, "TRUE", "true", "1"), , drop = FALSE], NULL, 15)
  })
  output$optional_metadata_table <- renderDT({
    req(active_inputs())
    md <- active_inputs()$metadata
    status_datatable(md[!(md$required %in% c(TRUE, "TRUE", "true", "1")), , drop = FALSE], NULL, 15)
  })
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
  output$threshold_table <- renderDT({ status_datatable(threshold_glossary, NULL, 20) })
  output$plate_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$plate_qc, "plate_qc_status", 10) })
  output$reference_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$reference_control_qc, "status", 15) })
  output$calibration_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$interplate_calibration, "calibration_status", 10) })
  output$primary_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$primary_results, "primary_status", 20) })
  output$sample_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$sample_qc_table, "sample_status", 20) })
  output$primary_agonist_table <- renderDT({
    req(analysis_results())
    z <- analysis_results()$primary_results
    status_datatable(z[z$assay_mode == "agonist", , drop = FALSE], "primary_status", 20)
  })
  output$primary_antagonist_table <- renderDT({
    req(analysis_results())
    z <- analysis_results()$primary_results
    status_datatable(z[z$assay_mode == "antagonist", , drop = FALSE], "primary_status", 20)
  })
  output$dose_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$dose_response_results, NULL, 20) })
  output$dose_qc_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$dose_response_qc, "curve_qc_status", 20) })
  output$dose_glossary_table <- renderDT(status_datatable(dose_glossary, NULL, 10))
  output$counter_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$counter_assay_qc, "counter_qc_status", 10) })
  output$hit_table <- renderDT({ req(analysis_results()); status_datatable(analysis_results()$hit_calls, "primary_status", 20) })
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
  output$primary_replicate_cv_plot <- renderPlotly(ggplotly(replicate_cv_obj()))
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
  output$download_replicate_cv_plot_primary <- download_plot(replicate_cv_obj, "primary_replicate_cv", active_inputs)
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
