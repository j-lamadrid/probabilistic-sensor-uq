#!/usr/bin/env Rscript
library(shiny)
library(plotly)
library(pracma)
dt_installed <- requireNamespace("DT", quietly = TRUE)
if (!dt_installed) message("Package 'DT' not installed; features tab will use a basic table.")
library(dplyr)

# ── 1 GB upload limit ──────────────────────────────────────────────────────────
options(shiny.maxRequestSize = 1024 * 1024 * 1024)

# Source helper functions
if (file.exists("R/utils.R"))   source("R/utils.R")
if (file.exists("R/plot_3d.R")) source("R/plot_3d.R")

# list available processed RDS files (if any)
processed_files <- character(0)
proc_dir <- file.path("data", "processed")
if (dir.exists(proc_dir)) {
  processed_files <- list.files(proc_dir, pattern = "\\.rds$", full.names = TRUE)
}
processed_choices <- if (length(processed_files) > 0)
  setNames(processed_files, basename(processed_files)) else NULL

# ── Aggregate per-spectrum peak rows into _agg style summary ──────────────────
# Mirrors the groupPeaks + summarize pipeline in model.R.
# Input:  per-spectrum data.frame with columns:
#           peak_mz, peak_height, peak_to_noise, quality_score,
#           retention_time  (at minimum)
# Output: aggregated data.frame grouped by mz cluster
aggregateFeatures <- function(df, mz_tolerance = 0.05, n_spectra_total = NULL) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  # Need groupPeaks from utils.R
  if (!exists("groupPeaks")) {
    warning("groupPeaks() not found; returning raw features.")
    return(df)
  }

  df <- df %>%
    mutate(peak_group_id = groupPeaks(peak_mz, mz_tolerance))

  agg <- df %>%
    group_by(peak_group_id) %>%
    summarize(
      n_spectra       = n(),
      mz_center       = mean(peak_mz),
      mz_min          = min(peak_mz),
      mz_max          = max(peak_mz),
      avg_peak_height = mean(peak_height),
      sd_peak_height  = sd(peak_height),
      avg_peak_to_noise = mean(peak_to_noise),
      avg_quality_score = mean(quality_score),
      rt_first        = min(retention_time),
      rt_last         = max(retention_time),
      .groups = "drop"
    ) %>%
    arrange(desc(avg_quality_score))

  # Optional artifact filtering mirroring model.R thresholds
  if (!is.null(n_spectra_total) && n_spectra_total > 0) {
    agg <- agg %>%
      filter(n_spectra > n_spectra_total * 0.01,
             n_spectra < n_spectra_total * 0.80)
  }

  agg
}

# ── 8-bit CSS ──────────────────────────────────────────────────────────────────
retro_css <- "
  @import url('https://fonts.googleapis.com/css2?family=Press+Start+2P&family=VT323:wght@400&display=swap');

  :root {
    --bg:        #1b1b1b;
    --bg2:       #242424;
    --panel:     #242424;
    --border:    #00ff41;
    --accent1:   #00ff41;
    --accent2:   #ff00ff;
    --accent3:   #00e5ff;
    --warn:      #ffdd00;
    --text:      #c8ffc8;
    --text-dim:  #4a8c4a;
    --pixel:     2px;
  }

  html, body {
    background: var(--bg) !important;
    color: var(--text) !important;
    font-family: 'VT323', monospace !important;
    font-size: 18px;
    overflow-x: hidden;
  }

  /* scanline overlay */
  body::after {
    content: '';
    position: fixed;
    inset: 0;
    background: repeating-linear-gradient(
      0deg,
      transparent,
      transparent 2px,
      rgba(0,0,0,0.18) 2px,
      rgba(0,0,0,0.18) 4px
    );
    pointer-events: none;
    z-index: 9999;
  }

  /* CRT vignette + flicker */
  body::before {
    content: '';
    position: fixed;
    inset: 0;
    background: radial-gradient(ellipse at center,
      rgba(0,255,65,0.04) 0%,
      transparent 70%,
      rgba(0,0,0,0.6) 100%);
    pointer-events: none;
    z-index: 9998;
    animation: flicker 8s infinite;
  }

  @keyframes flicker {
    0%,100% { opacity: 1; }
    92%      { opacity: 1; }
    93%      { opacity: 0.85; }
    94%      { opacity: 1; }
    96%      { opacity: 0.9; }
    97%      { opacity: 1; }
  }

  h2, h2.shiny-app-title {
    font-family: 'Press Start 2P', monospace !important;
    font-size: 14px !important;
    color: var(--accent1) !important;
    text-shadow: 0 0 8px var(--accent1), 0 0 20px var(--accent1) !important;
    letter-spacing: 2px;
    padding: 16px 0 12px !important;
    border-bottom: 2px solid var(--border);
    margin-bottom: 16px !important;
  }

  .well {
    background: var(--panel) !important;
    border: 2px solid var(--border) !important;
    border-radius: 0 !important;
    box-shadow: 4px 4px 0 var(--accent1) !important;
    padding: 16px !important;
  }

  label, .control-label {
    font-family: 'Press Start 2P', monospace !important;
    font-size: 7px !important;
    color: var(--accent3) !important;
    letter-spacing: 1px;
    text-transform: uppercase;
    margin-bottom: 6px !important;
  }

  input[type='number'], input[type='text'], select, .form-control {
    background: #000 !important;
    color: var(--accent1) !important;
    border: 2px solid var(--border) !important;
    border-radius: 0 !important;
    font-family: 'VT323', monospace !important;
    font-size: 18px !important;
    padding: 4px 8px !important;
    box-shadow: 2px 2px 0 var(--text-dim) !important;
    outline: none !important;
  }
  input[type='number']:focus, input[type='text']:focus,
  select:focus, .form-control:focus {
    box-shadow: 2px 2px 0 var(--accent1), 0 0 10px var(--accent1) !important;
  }

  input[type='checkbox'] { accent-color: var(--accent1); width:16px; height:16px; }

  .btn, .btn-default {
    font-family: 'Press Start 2P', monospace !important;
    font-size: 8px !important;
    background: #000 !important;
    color: var(--accent1) !important;
    border: 2px solid var(--accent1) !important;
    border-radius: 0 !important;
    padding: 10px 14px !important;
    box-shadow: 4px 4px 0 var(--accent1) !important;
    text-transform: uppercase;
    letter-spacing: 1px;
    transition: all 0.05s !important;
  }
  .btn:hover, .btn-default:hover {
    background: var(--accent1) !important;
    color: #000 !important;
    box-shadow: 2px 2px 0 var(--text-dim) !important;
    transform: translate(2px, 2px);
  }
  .btn:active, .btn-default:active {
    transform: translate(4px, 4px) !important;
    box-shadow: none !important;
  }

  .btn-file {
    font-family: 'Press Start 2P', monospace !important;
    font-size: 7px !important;
    background: #000 !important;
    color: var(--accent2) !important;
    border: 2px solid var(--accent2) !important;
    border-radius: 0 !important;
    box-shadow: 3px 3px 0 var(--accent2) !important;
  }
  .btn-file:hover {
    background: var(--accent2) !important;
    color: #000 !important;
    transform: translate(2px,2px);
    box-shadow: 1px 1px 0 #880088 !important;
  }
  .form-group .input-group .form-control[readonly] {
    font-family: 'VT323', monospace !important;
    font-size: 16px !important;
    color: var(--accent2) !important;
    border-color: var(--accent2) !important;
    box-shadow: none !important;
  }

  hr {
    border: none !important;
    border-top: 2px dashed var(--text-dim) !important;
    margin: 14px 0 !important;
  }

  .help-block {
    font-family: 'VT323', monospace !important;
    font-size: 14px !important;
    color: var(--text-dim) !important;
  }

  .nav-tabs { border-bottom: 2px solid var(--border) !important; }
  .nav-tabs > li > a {
    font-family: 'Press Start 2P', monospace !important;
    font-size: 7px !important;
    color: var(--text-dim) !important;
    background: #000 !important;
    border: 2px solid var(--text-dim) !important;
    border-bottom: none !important;
    border-radius: 0 !important;
    margin-right: 4px !important;
    padding: 8px 12px !important;
    letter-spacing: 1px;
  }
  .nav-tabs > li > a:hover {
    color: var(--accent1) !important;
    border-color: var(--accent1) !important;
    background: #0a1a0a !important;
    text-decoration: none;
  }
  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus,
  .nav-tabs > li.active > a:hover {
    color: #000 !important;
    background: var(--accent1) !important;
    border-color: var(--accent1) !important;
  }
  .tab-content {
    border: 2px solid var(--border) !important;
    border-top: none !important;
    background: var(--panel) !important;
    padding: 12px !important;
  }

  .progress {
    background: #000 !important;
    border: 2px solid var(--border) !important;
    border-radius: 0 !important;
    height: 20px !important;
  }
  .progress-bar {
    background-image: repeating-linear-gradient(
      90deg, var(--accent1) 0px, var(--accent1) 8px,
      #000 8px, #000 10px) !important;
    transition: width 0.1s steps(10) !important;
    box-shadow: 0 0 10px var(--accent1) !important;
  }

  pre {
    background: #000 !important;
    color: var(--accent1) !important;
    border: 2px solid var(--border) !important;
    border-radius: 0 !important;
    font-family: 'VT323', monospace !important;
    font-size: 18px !important;
    padding: 12px !important;
    box-shadow: 3px 3px 0 var(--text-dim) !important;
  }
  pre::before { content: '> '; color: var(--accent2); }

  /* DT table */
  .dataTables_wrapper {
    font-family: 'VT323', monospace !important;
    font-size: 16px !important;
    color: var(--text) !important;
  }
  table.dataTable {
    background: #000 !important;
    border: 2px solid var(--border) !important;
    border-radius: 0 !important;
  }
  table.dataTable thead th {
    font-family: 'Press Start 2P', monospace !important;
    font-size: 7px !important;
    color: var(--accent1) !important;
    background: #0a1a0a !important;
    border-bottom: 2px solid var(--border) !important;
    padding: 10px 8px !important;
    letter-spacing: 1px;
  }
  table.dataTable tbody tr   { background: #000 !important; color: var(--text) !important; }
  table.dataTable tbody tr:hover { background: #0a1a0a !important; color: var(--accent1) !important; }
  table.dataTable tbody tr.odd  { background: #000 !important; }
  table.dataTable tbody tr.even { background: #05100a !important; }
  table.dataTable tbody td { border-top: 1px solid var(--text-dim) !important; padding: 6px 8px !important; }
  .dataTables_filter input {
    background: #242424 !important; color: var(--accent1) !important;
    border: 2px solid var(--border) !important; border-radius: 0 !important;
    font-family: 'VT323', monospace !important; font-size: 16px !important;
  }
  .dataTables_paginate .paginate_button {
    font-family: 'Press Start 2P', monospace !important; font-size: 7px !important;
    color: var(--text-dim) !important; background: #000 !important;
    border: 2px solid var(--text-dim) !important; border-radius: 0 !important; margin: 2px !important;
  }
  .dataTables_paginate .paginate_button.current,
  .dataTables_paginate .paginate_button:hover {
    color: #000 !important; background: var(--accent1) !important; border-color: var(--accent1) !important;
  }
  .dataTables_info, .dataTables_length label {
    font-family: 'VT323', monospace !important; font-size: 16px !important; color: var(--text-dim) !important;
  }
  .dataTables_length select {
    background: #000 !important; color: var(--accent1) !important;
    border: 2px solid var(--border) !important; border-radius: 0 !important;
    font-family: 'VT323', monospace !important; font-size: 16px !important;
  }

  .js-plotly-plot, .plot-container {
    border: 2px solid var(--border) !important;
    background: #242424 !important;
    box-shadow: 4px 4px 0 var(--text-dim);
  }

  .selectize-control .selectize-input {
    background: #242424 !important; border: 2px solid var(--border) !important;
    border-radius: 0 !important; color: var(--accent1) !important;
    font-family: 'VT323', monospace !important; font-size: 18px !important;
    box-shadow: 2px 2px 0 var(--text-dim) !important;
  }
  .selectize-control .selectize-dropdown {
    background: #242424 !important; border: 2px solid var(--border) !important;
    border-radius: 0 !important; color: var(--text) !important;
    font-family: 'VT323', monospace !important; font-size: 18px !important;
  }
  .selectize-dropdown .option:hover,
  .selectize-dropdown .option.active { background: var(--accent1) !important; color: #000 !important; }
"

# ── UI ─────────────────────────────────────────────────────────────────────────
tabs_list <- list(
  tabPanel(">> 2D (FINGERPRINT)", plotlyOutput("plot2d", height = "680px")),
  tabPanel(">> 3D (TEMPORAL)", plotlyOutput("plot3d", height = "680px"))
)
if (dt_installed) {
  tabs_list[[length(tabs_list) + 1]] <- tabPanel(">> FEATURES", DT::DTOutput("features_table"))
} else {
  tabs_list[[length(tabs_list) + 1]] <- tabPanel(">> FEATURES", tableOutput("features_table_base"))
}
tabs_list[[length(tabs_list) + 1]] <- tabPanel(">> STATUS", verbatimTextOutput("status"))

ui <- fluidPage(
  tags$head(tags$style(HTML(retro_css))),

  titlePanel(
    tags$span(
      style = "font-family:'Press Start 2P',monospace; font-size:12px;
               color:#00ff41; text-shadow:0 0 10px #00ff41, 0 0 25px #00ff41;
               letter-spacing:3px;",
      "▌LC-MS ANALYZER▐"
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$div(
        style = "font-family:'Press Start 2P',monospace; font-size:6px;
                 color:#000; background:#00ff41; display:inline-block;
                 padding:3px 6px; margin-bottom:10px; letter-spacing:1px;",
        "MAX: 1 GB"
      ),
      fileInput("mzml", "UPLOAD .mzML FILE", accept = c('.mzML', '.mzML.gz')),
      tags$hr(),
      numericInput("bin_width",   "M/Z BIN WIDTH",          value = 0.01, step = 0.001),
      numericInput("max_spectra", "MAX SPECTRA TO PROCESS", value = 200,  min = 1),
      numericInput("window_size", "PEAK WINDOW SIZE (M/Z)", value = 0.8,  step = 0.1),
      numericInput("mz_tolerance","MZ GROUP TOLERANCE",     value = 0.05, step = 0.005),
      checkboxInput("log_scale",  "LOG INTENSITY (LOG1P)",  value = FALSE),
      tags$div(style = "margin-top:14px;"),
      actionButton("analyze", "[ ANALYZE ]"),
      tags$hr(),
      selectInput("processed_rds", "SELECT PROCESSED RDS",
                  choices = processed_choices, selected = NULL),
      helpText("LOAD .mzML + ANALYZE, OR SELECT A PRE-PROCESSED RDS.")
    ),

    mainPanel(do.call(tabsetPanel, tabs_list))
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  features_tbl <- reactiveVal(NULL)
  current_plot_2d <- reactiveVal(NULL)
  current_plot_3d <- reactiveVal(NULL)
  status_msg   <- reactiveVal("IDLE - AWAITING INPUT...")
  n_spectra_processed <- reactiveVal(0L)

  # ── Load pre-processed RDS ──────────────────────────────────────────────────
  observeEvent(input$processed_rds, {
    req(input$processed_rds)
    pth <- input$processed_rds
    if (!file.exists(pth)) return(NULL)
    status_msg(paste("LOADING RDS:", basename(pth)))
    df <- tryCatch(readRDS(pth),
                   error = function(e) { status_msg(paste("RDS ERROR:", e$message)); NULL })
    if (!is.null(df)) {
      features_tbl(df)
      status_msg(paste("RDS LOADED:", nrow(df), "ROWS"))
    }
  }, ignoreNULL = TRUE)

  # ── Analyze button ──────────────────────────────────────────────────────────
  observeEvent(input$analyze, {
    req(input$mzml)
    fileinfo  <- input$mzml
    mzml_path <- paste0(fileinfo$datapath, ".mzML")
    file.rename(fileinfo$datapath, mzml_path)
    mzml_name <- fileinfo$name

    if (!exists("loadMZML") || !exists("binData") || !exists("findRegions") || !exists("generatePlotsFromFeatures")) {
      status_msg("ERROR: HELPER FUNCTIONS NOT FOUND. CHECK R/utils.R AND R/plot_3d.R")
      return(NULL)
    }

    status_msg(paste("LOADING:", mzml_name))

    withProgress(message = "PROCESSING mzML...", value = 0, {
      incProgress(0.05)

      features_obj <- tryCatch(
        loadMZML(mzml_path, tools::file_path_sans_ext(mzml_name)),
        error = function(e) { status_msg(paste("LOAD ERROR:", e$message)); NULL }
      )
      if (is.null(features_obj)) return(NULL)

      n_spec <- nrow(features_obj)
      n_use  <- min(n_spec, as.integer(input$max_spectra))
      n_spectra_processed(n_use)
      status_msg(paste("SPECTRA FOUND:", n_spec, "| PROCESSING:", n_use))
      incProgress(0.1)

      # ── Per-spectrum peak detection ────────────────────────────────────────
      res_list <- list()
      for (i in seq_len(n_use)) {
        binned <- tryCatch(
          binData(features_obj, i, bin_width = input$bin_width),
          error = function(e) NULL
        )
        if (is.null(binned) || nrow(binned) == 0) next

        regs <- tryCatch(
          findRegions(binned, window_size = input$window_size),
          error = function(e) NULL
        )
        if (is.null(regs)) next

        for (r in seq_along(regs)) {
          reg <- regs[[r]]
          res_list[[length(res_list) + 1]] <- data.frame(
            spectrum_id    = i,
            retention_time = features_obj$retention_time[i],
            peak_mz        = reg$peak_mz,
            peak_height    = reg$peak_height,
            mz_start       = reg$mz_range[1],
            mz_end         = reg$mz_range[2],
            n_points       = reg$n_points,
            peak_to_noise  = reg$peak_to_noise,
            quality_score  = reg$quality_score,
            stringsAsFactors = FALSE
          )
        }
        if (i %% 10 == 0) incProgress(0.5 * i / n_use)
      }

      if (length(res_list) == 0) {
        status_msg("NO FEATURES DETECTED")
        features_tbl(data.frame())
      } else {
        raw_df <- dplyr::bind_rows(res_list)
        status_msg(sprintf("DETECTED %d RAW PEAKS - AGGREGATING BY M/Z GROUP...", nrow(raw_df)))

        # ── Aggregate to _agg style ──────────────────────────────────────────
        agg_df <- tryCatch(
          aggregateFeatures(raw_df,
                            mz_tolerance    = input$mz_tolerance,
                            n_spectra_total = n_use),
          error = function(e) {
            status_msg(paste("AGG ERROR:", e$message))
            raw_df  # fall back to raw
          }
        )
        features_tbl(agg_df)
        status_msg(sprintf(
          "COMPLETE: %d RAW PEAKS -> %d M/Z GROUPS",
          nrow(raw_df), nrow(agg_df)
        ))
      }

      incProgress(0.1)

      # ── Build 2D & 3D plots from already-loaded features_obj ────────
      status_msg(paste(status_msg(), "| BUILDING PLOTS..."))
      
      plots <- tryCatch(
        generatePlotsFromFeatures(
          features_obj,
          bin_width   = input$bin_width,
          max_spectra = input$max_spectra,
          log_scale   = input$log_scale
        ),
        error   = function(e) { status_msg(paste("PLOT ERROR:",   e$message)); NULL },
        warning = function(w) { status_msg(paste("PLOT WARNING:", w$message)); NULL }
      )
      
      if (!is.null(plots)) {
        current_plot_2d(plots$plot2d)
        current_plot_3d(plots$plot3d)
      }

      final_status <- if (is.null(plots)) {
        "WARNING: PLOTS NULL - TRY INCREASING BIN WIDTH"
      } else {
        sprintf("DONE: %d M/Z GROUPS | PLOTS OK", nrow(features_tbl()))
      }
      status_msg(final_status)
      incProgress(0.25)
    })
  })

  # ── Render 2D Plot ──────────────────────────────────────────────────────────
  output$plot2d <- renderPlotly({
    p <- current_plot_2d()
    empty_layout <- list(
      paper_bgcolor = "#242424", plot_bgcolor  = "#242424",
      font          = list(color = "#00ff41", family = "monospace"),
      xaxis         = list(title = "M/Z",       color = "#00ff41", gridcolor = "#1a3a1a"),
      yaxis         = list(title = "INTENSITY", color = "#00ff41", gridcolor = "#1a3a1a")
    )
    
    apply_retro_layout <- function(plt, extra_args = list()) {
      do.call(plotly::layout, c(list(plt), empty_layout, extra_args))
    }
    
    if (is.null(p)) {
      apply_retro_layout(
        plotly::plot_ly(type = "scatter", mode = "lines"), 
        list(title = list(text = "UPLOAD AN mzML FILE AND CLICK [ ANALYZE ]", font = list(family = "monospace", color = "#00ff41", size = 14)))
      )
    } else {
      apply_retro_layout(p)
    }
  })

  # ── Render 3D Plot ────────────────────────────────────────────────────────────
  output$plot3d <- renderPlotly({
    p <- current_plot_3d()
    empty_layout <- list(
      paper_bgcolor = "#242424", plot_bgcolor  = "#242424",
      font          = list(color = "#00ff41", family = "monospace"),
      scene = list(
        bgcolor = "#242424",
        xaxis   = list(title = "M/Z",       color = "#00ff41", gridcolor = "#1a3a1a"),
        yaxis   = list(title = "TIME (RT)", color = "#00ff41", gridcolor = "#1a3a1a"),
        zaxis   = list(title = "INTENSITY", color = "#00ff41", gridcolor = "#1a3a1a")
      )
    )
    
    apply_retro_layout <- function(plt, extra_args = list()) {
      do.call(plotly::layout, c(list(plt), empty_layout, extra_args))
    }
    
    if (is.null(p)) {
      apply_retro_layout(
        plotly::plot_ly(type = "scatter3d", mode = "lines"),
        list(title = list(text = "UPLOAD AN mzML FILE AND CLICK [ ANALYZE ]", font = list(family = "monospace", color = "#00ff41", size = 14)))
      )
    } else {
      apply_retro_layout(p)
    }
  })

  # ── Render features table (_agg style) ─────────────────────────────────────
  fmt_cols <- function(df) {
    # Round numeric columns for readability
    num_cols <- sapply(df, is.numeric)
    df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
    df
  }

  empty_agg_msg <- data.frame(
    MESSAGE = "NO FEATURES. UPLOAD mzML + ANALYZE OR SELECT AN RDS."
  )

  if (dt_installed) {
    output$features_table <- DT::renderDT({
      df <- features_tbl()
      tbl <- if (is.null(df) || nrow(df) == 0) empty_agg_msg else fmt_cols(df)
      DT::datatable(
        tbl,
        rownames = FALSE,
        options  = list(
          pageLength = 15,
          scrollX    = TRUE,
          scrollY    = "450px",
          initComplete = DT::JS(
            "function(settings,json){",
            "  $(this.api().table().header()).css({'background-color':'#0a1a0a','color':'#00ff41'});",
            "}"
          )
        )
      )
    })
  } else {
    output$features_table_base <- renderTable({
      df <- features_tbl()
      if (is.null(df) || nrow(df) == 0) empty_agg_msg else fmt_cols(df)
    }, striped = TRUE, hover = TRUE)
  }

  output$status <- renderText({ status_msg() })
}

shinyApp(ui, server)
