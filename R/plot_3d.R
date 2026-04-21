# Helper: 3D LC-MS plotting (x = m/z, y = intensity, z = time)

generatePlotsFromFeatures <- function(features,
                                      mz_min      = NULL,
                                      mz_max      = NULL,
                                      bin_width   = 0.01,
                                      max_spectra = NULL,
                                      log_scale   = FALSE,
                                      point_size  = 2, 
                                      opacity     = 0.8) {

  if (!exists("binData"))
    stop("Function 'binData' not found. Source utils.R first.")

  n_spec <- nrow(features)
  if (n_spec == 0) stop("No spectra in features data frame.")

  n_use <- if (!is.null(max_spectra)) min(n_spec, as.integer(max_spectra)) else n_spec
  features <- features[seq_len(n_use), , drop = FALSE]

  # Process and bin data exactly ONCE for both plots
  rows <- vector("list", n_use)
  for (i in seq_len(n_use)) {
    binned <- tryCatch(binData(features, i, bin_width = bin_width), error = function(e) NULL)
    if (is.null(binned) || nrow(binned) == 0L) next
    binned$time <- features$retention_time[i]
    rows[[i]]   <- binned
  }

  non_null <- rows[!vapply(rows, is.null, logical(1))]
  if (length(non_null) == 0L) {
    warning("No non-zero intensities after binning.")
    return(NULL)
  }

  df_long <- do.call(rbind, non_null)

  if (!is.null(mz_min)) df_long <- df_long[df_long$mz >= mz_min, ]
  if (!is.null(mz_max)) df_long <- df_long[df_long$mz <= mz_max, ]

  if (nrow(df_long) == 0L) {
    warning("No data in specified mz range.")
    return(NULL)
  }

  # Apply log scale globally if requested
  if (log_scale) df_long$intensity <- log1p(df_long$intensity)

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Please install the 'plotly' package.")

  col_pal <- if (requireNamespace("viridis", quietly = TRUE)) viridis::viridis(256) else NULL

  # --- 1. GENERATE 2D PROFILE PLOT ---
  df_profile <- aggregate(intensity ~ mz, data = df_long, FUN = max)
  df_profile <- df_profile[order(df_profile$mz), ]

  p2d <- plotly::plot_ly(
    df_profile,
    x      = ~mz,
    y      = ~intensity,
    type   = "scatter",
    mode   = "lines",
    line   = list(color = "#00ff41", width = 1.5)
  ) %>%
    plotly::layout(
      xaxis = list(title = "m/z"),
      yaxis = list(title = "Intensity")
    )

  # --- 2. GENERATE 3D WATERFALL PLOT ---
  df_long <- df_long[order(df_long$time, df_long$mz), ]

  p3d <- plotly::plot_ly(
    df_long,
    x      = ~mz,
    y      = ~time,           
    z      = ~intensity,      
    split  = ~time,           
    color  = ~intensity,      
    colors = col_pal,
    type   = "scatter3d",
    mode   = "lines",
    line   = list(width = 2),
    hoverinfo = "x+y+z"
  ) %>%
    plotly::layout(
      showlegend = FALSE,     
      scene = list(
        xaxis = list(title = "m/z"),
        yaxis = list(title = "Time (RT)"),
        zaxis = list(title = "Intensity")
      )
    )

  # Return both plots safely
  return(list(plot2d = p2d, plot3d = p3d))
}

plot3DLCMS <- function(raw_mzml_path,
                       mz_min      = NULL,
                       mz_max      = NULL,
                       bin_width   = 0.01,
                       max_spectra = NULL,
                       log_scale   = FALSE,
                       point_size  = 2,
                       opacity     = 0.8) {

  if (!file.exists(raw_mzml_path))
    stop(paste("File not found:", raw_mzml_path))
  if (!exists("loadMZML"))
    stop("Function 'loadMZML' not found. Source utils.R first.")

  cultivar_name <- tools::file_path_sans_ext(basename(raw_mzml_path))
  features      <- loadMZML(raw_mzml_path, cultivar_name)

  plot3DFromFeatures(features,
                     mz_min      = mz_min,
                     mz_max      = mz_max,
                     bin_width   = bin_width,
                     max_spectra = max_spectra,
                     log_scale   = log_scale,
                     point_size  = point_size,
                     opacity     = opacity)
}
