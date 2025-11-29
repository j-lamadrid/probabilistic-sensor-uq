library(rstan)
library(dplyr)
library(ggplot2)
library(pracma)
library(here)
library(ggrepel)
library(GGally)
library(DiagrammeR)

source("R/utils.R")


plotMSAgg <- function(name) {
  if (!file.exists(name)) {
    warning(paste("File not found:", name))
    next
  }
  
  final_fingerprint_table = readRDS(name)
  plot_data = final_fingerprint_table %>%
    filter(mz_r2 > 0.75)
  
  target_cultivar = plot_data$cultivar[1]
  
  SIGMA = 0.2
  BASELINE = 100
  
  x_fine = seq(min(plot_data$mz_position) - 10,
               max(plot_data$mz_position) + 10,
               length.out = 30000)
  
  peak_components_list = list()
  
  for (i in 1:nrow(plot_data)) {
    mu = plot_data$mz_position[i]
    amp = plot_data$mz_amplitude[i]
    r2 = plot_data$mz_r2[i]
    rho = plot_data$mz_precision[i]
    region_label = plot_data$mz_center[i]
    
    peak_signal = amp * exp(-0.5 * rho * (x_fine - mu)^2)
    
    peak_components_list[[i]] = data.frame(
      mz = x_fine,
      intensity = peak_signal + BASELINE,
      mz_r2 = r2,
      peak_region = region_label
    )
  }
  
  synthetic_peaks_df = bind_rows(peak_components_list)
  
  label_data = plot_data %>%
    arrange(desc(mz_amplitude)) %>%
    head(10)
  
  p = ggplot() +
    
    geom_line(
      data = synthetic_peaks_df,
      aes(
        x = mz,
        y = intensity,
        group = peak_region,
        color = mz_r2
      ),
      linewidth = 0.8
    ) +
    
    geom_text_repel(
      data = label_data,
      aes(
        x = mz_position,
        y = mz_amplitude,
        label = round(mz_center, 2)
      ),
      size = 3.0,
      fontface = "bold",
      direction = "y",
      segment.size = 0.5,
      segment.color = "grey50",
      nudge_y = 200,
    ) +
    
    scale_color_viridis_c(option = "mako",
                          name = "Confidence (R2)",
                          direction = -1) +
    
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    
    labs(
      title = paste("Reconstructed Mass Spectrum:", target_cultivar),
      subtitle = "Probabilistic Collapsed Signal Reconstruction",
      x = "m/z",
      y = "Intensity",
      alpha = "R2"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold")
    )
  
  print(p)
}



plotPrec <- function(raw_name, n = 10) {
  if (!file.exists(raw_name)) {
    warning(paste("File not found:", raw_name))
    next
  }
  
  raw_data = readRDS(raw_name)
  MZ_TOLERANCE = 0.05
  data = raw_data %>%
    mutate(peak_group_id = groupPeaks(peak_mz_prior, MZ_TOLERANCE))
  
  significant_peaks_summary <- data %>%
    group_by(peak_group_id) %>%
    summarize(
      Total_Amplitude = sum(model_amplitude_a),
      Avg_R2 = mean(model_r2),
      Mz_Center = mean(model_position_mu),
      .groups = 'drop'
    ) %>%
    filter(Avg_R2 > 0.8, Mz_Center > 80) %>%
    arrange(desc(Total_Amplitude)) %>%
    slice_head(n = n)
  
  plot_data <- data %>%
    inner_join(significant_peaks_summary, by = "peak_group_id") %>%
    mutate(Peak_Label = paste0("m/z ", sprintf("%.3f", Mz_Center))) %>%
    mutate(Peak_Label = reorder(Peak_Label, Mz_Center))
  
  drift_data <- plot_data %>%
    group_by(peak_group_id) %>%
    mutate(
      avg_group_mz = mean(model_position_mu),
      ppm_shift = (model_position_mu - avg_group_mz) / avg_group_mz * 1e6
    ) %>%
    ungroup()
  
  p_drift <- ggplot(drift_data, aes(x = ppm_shift, y = Peak_Label, fill = Avg_R2)) +
    
    geom_density_ridges(
      alpha = 0.8,
      scale = 0.9,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "black",
      linewidth = 0.5
    ) +
    geom_vline(
      xintercept = c(-5, 5),
      linetype = "dashed",
      color = "darkgreen",
      alpha = 0.6
    ) +
    geom_vline(
      xintercept = c(-20, 20),
      linetype = "dotted",
      color = "red",
      alpha = 0.6
    ) +
    
    xlim(-30, 30) +
    scale_fill_viridis_c(option = "mako",
                         name = "Confidence (R2)",
                         direction = -1) +
    
    labs(
      title = "Sensor Precision",
      subtitle = "Most Significant Peaks Sorted by m/z Center",
      x = "Mass Drift (ppm)",
      y = "Peak Region"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  print(p_drift)
}



plotEvidence <- function(raw_name, n = 10) {
  if (!file.exists(raw_name)) {
    warning(paste("File not found:", raw_name))
    next
  }
  
  raw_data = readRDS(raw_name)
  MZ_TOLERANCE = 0.05
  data = raw_data %>%
    mutate(peak_group_id = groupPeaks(peak_mz_prior, MZ_TOLERANCE))
  
  significant_peaks_summary <- data %>%
    group_by(peak_group_id) %>%
    summarize(
      Total_Amplitude = sum(model_amplitude_a),
      Avg_R2 = mean(model_r2),
      Mz_Center = mean(model_position_mu),
      .groups = 'drop'
    ) %>%
    filter(Avg_R2 > 0.8, Mz_Center > 80) %>%
    arrange(desc(Total_Amplitude)) %>%
    slice_head(n = n)
  
  plot_data <- data %>%
    inner_join(significant_peaks_summary, by = "peak_group_id") %>%
    mutate(Peak_Label = paste0("m/z ", sprintf("%.3f", Mz_Center))) %>%
    mutate(Peak_Label = reorder(Peak_Label, Mz_Center))
  
  p_evidence <- ggplot(plot_data,
                       aes(x = bayes_factor_evidence, y = Peak_Label, fill = Avg_R2)) +
    
    # Draw the densities
    geom_density_ridges(
      alpha = 0.8,
      scale = 0.9,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    
    xlim(-10, 150) +
    
    geom_vline(
      xintercept = 2,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    ) +
    
    scale_fill_viridis_c(option = "mako",
                         name = "Confidence (R2)",
                         direction = -1) +
    
    labs(
      title = "Bayesian Model Selection: Evidence Distribution",
      subtitle = "BIC Difference > 2 indicates strong evidence for Peak vs. Noise",
      x = "BIC(Noise) - BIC(Peak)",
      y = "Peak Region"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  print(p_evidence)
}



plotAmp <- function(raw_name, n = 10) {
  if (!file.exists(raw_name)) {
    warning(paste("File not found:", raw_name))
    next
  }
  
  raw_data = readRDS(raw_name)
  MZ_TOLERANCE = 0.05
  data = raw_data %>%
    mutate(peak_group_id = groupPeaks(peak_mz_prior, MZ_TOLERANCE))
  
  significant_peaks_summary <- data %>%
    group_by(peak_group_id) %>%
    summarize(
      Total_Amplitude = sum(model_amplitude_a),
      Avg_R2 = mean(model_r2),
      Mz_Center = mean(model_position_mu),
      .groups = 'drop'
    ) %>%
    filter(Avg_R2 > 0.8, Mz_Center > 80) %>%
    arrange(desc(Total_Amplitude)) %>%
    slice_head(n = n)
  
  plot_data <- data %>%
    inner_join(significant_peaks_summary, by = "peak_group_id") %>%
    mutate(Peak_Label = paste0("m/z ", sprintf("%.3f", Mz_Center))) %>%
    mutate(Peak_Label = reorder(Peak_Label, Mz_Center))
  
  p_amp <- ggplot(plot_data,
                  aes(x = model_amplitude_a, y = Peak_Label, fill = Avg_R2)) +
    
    geom_density_ridges(
      scale = 0.9,
      quantile_lines = TRUE,
      quantiles = 2,
      alpha = 0.8
    ) +
    
    scale_x_log10() +
    
    scale_fill_viridis_c(option = "mako",
                         name = "Confidence (R2)",
                         direction = -1) +
    
    labs(
      title = "Posterior Amplitude Distributions",
      subtitle = "Abundance quantification (Log Scale). Color indicates Model Confidence.",
      x = "Amplitude (Log Scale)",
      y = "Peak Region (Ranked by Total Abundance)"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  print(p_amp)
}



plotPairs <- function(processed_rds_path,
                      raw_mzml_path,
                      model_object) {
  if (!file.exists(processed_rds_path))
    stop(paste("Processed file not found:", processed_rds_path))
  if (!file.exists(raw_mzml_path))
    stop(paste("Raw .mzML file not found:", raw_mzml_path))
  
  raw_df <- readRDS(processed_rds_path)
  cultivar_name <- tools::file_path_sans_ext(basename(raw_mzml_path))
  
  if (!exists("loadMZML"))
    stop("Function 'loadMZML' not found. Source utils.R first.")
  features_obj <- loadMZML(raw_mzml_path, cultivar_name)
  
  target_row <- raw_df %>%
    filter(model_r2 > 0.9) %>%
    arrange(desc(model_amplitude_a)) %>%
    slice_head(n = 1)
  
  if (nrow(target_row) == 0) {
    target_row <- raw_df %>% arrange(desc(model_r2)) %>% slice_head(n = 1)
  }
  
  cat(
    sprintf(
      "Generating Pairs Plot for Peak m/z %.4f (R2=%.2f)...\n",
      target_row$peak_mz_prior,
      target_row$model_r2
    )
  )
  
  spec <- binData(features_obj, target_row$spectrum_id, bin_width = 0.015) %>%
    filter(mz > target_row$peak_mz_prior - 0.4,
           mz < target_row$peak_mz_prior + 0.4)
  
  height_prior <- if ("peak_height_prior" %in% colnames(target_row)) {
    target_row$peak_height_prior
  } else {
    max(spec$intensity)
  }
  
  sdata <- list(
    N = nrow(spec),
    x = as.array(spec$mz),
    y = as.array(spec$intensity),
    x_min = min(spec$mz),
    x_max = max(spec$mz),
    peak_mz_prior = target_row$peak_mz_prior,
    peak_height_prior = height_prior,
    baseline_prior = median(spec$intensity),
    noise_prior = sd(spec$intensity),
    K = 1
  )
  
  fit <- suppressWarnings(vb(
    model_object,
    data = sdata,
    seed = 123,
    refresh = 0,
    output_samples = 5000
  ))
  
  post_df <- as.data.frame(fit)
  cols_to_plot <- c("a", "mu", "rho")
  
  if (!all(cols_to_plot %in% colnames(post_df))) {
    warning("Parameters a, mu, rho not found in output.")
    return(NULL)
  }
  
  plot_data <- post_df[, cols_to_plot]
  colnames(plot_data) <- c("Amplitude (a)", "Position (μ)", "Precision (ρ)")
  
  p <- ggpairs(
    plot_data,
    
    lower = list(continuous = wrap(
      "density", color = "black", alpha = 0.5
    )),
    
    diag = list(continuous = wrap(
      "barDiag",
      fill = "black",
      alpha = 0.3,
      bins = 30
    )),
    
    upper = list(continuous = wrap("cor", size = 4)),
    
    title = paste(
      "Joint Posterior Distributions: m/z",
      round(target_row$peak_mz_prior, 3)
    )
  ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 11),
      # Axis labels
      panel.grid.minor = element_blank()
    )
  
  print(p)
}



plotPPC <- function(processed_rds_path, raw_mzml_path) {
  if (!file.exists(processed_rds_path))
    stop(paste("Processed file not found:", processed_rds_path))
  if (!file.exists(raw_mzml_path))
    stop(paste("Raw .mzML file not found:", raw_mzml_path))
  
  raw_df <- readRDS(processed_rds_path)
  cultivar_name <- tools::file_path_sans_ext(basename(raw_mzml_path))
  
  if (!exists("loadMZML"))
    stop("Function 'loadMZML' not found. Source utils.R first.")
  
  cat(sprintf("Loading raw spectral data from %s...\n", raw_mzml_path))
  features_obj <- loadMZML(raw_mzml_path, cultivar_name)
  
  
  best <- raw_df %>%
    filter(model_r2 > 0.9) %>%
    arrange(desc(model_r2)) %>%
    slice_head(n = 1)
  
  worst <- raw_df %>%
    filter(model_r2 < 0.7, model_amplitude_a > 1000) %>%
    arrange(model_r2) %>%
    slice_head(n = 1)
  
  if (nrow(worst) == 0) {
    worst <- raw_df %>%
      filter(model_amplitude_a > 500) %>%
      arrange(model_r2) %>%
      slice_head(n = 1)
  }
  
  targets <- bind_rows(best, worst)
  
  if (nrow(targets) == 0) {
    warning("No valid peaks found to plot.")
    return(NULL)
  }
  
  
  plot_list <- list()
  
  for (i in seq_len(nrow(targets))) {
    row <- targets[i, ]
    
    raw_bin <- binData(features_obj, row$spectrum_id, bin_width = 0.015) %>%
      filter(mz > row$peak_mz_prior - 0.2, mz < row$peak_mz_prior + 0.2)
    
    if (nrow(raw_bin) > 0) {
      sim <- raw_bin %>%
        mutate(
          y_hat = row$model_baseline + row$model_amplitude_a * exp(
            -0.5 * row$model_precision_rho * (mz - row$model_position_mu)^2
          ),
          
          sigma = sqrt(row$noise_variance),
          
          ymin = pmax(y_hat - 2 * sigma, 0),
          
          ymax = y_hat + 2 * sigma,
          
          Label = ifelse(
            row$model_r2 > 0.8,
            paste0("High Confidence (R2=", round(row$model_r2, 2), ")"),
            paste0("High Uncertainty (R2=", round(row$model_r2, 2), ")")
          )
        )
      plot_list[[i]] <- sim
    }
  }
  
  plot_final <- bind_rows(plot_list)
  if (nrow(plot_final) == 0)
    return(NULL)
  
  p <- ggplot(plot_final, aes(x = mz)) +
    
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
                fill = "black",
                alpha = 0.2) +
    geom_point(aes(y = intensity), alpha = 0.5, size = 1.5) +
    geom_line(aes(y = y_hat),
              color = "darkgrey",
              linewidth = 1) +
    
    facet_wrap( ~ Label, scales = "free") +
    
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    
    labs(
      title = "Posterior Predictive Check",
      subtitle = "Ribbon = 95% Uncertainty.",
      x = "m/z",
      y = "Intensity"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 10)
    )
  
  print(p)
}



plotGraph <- function() {
  
    graph = grViz("
      digraph {
        graph [rankdir = TB, layout = dot, splines=spline, nodesep=0.6, ranksep=0.6]
    
        node [fontname = 'Helvetica', fontsize=12, penwidth=1.0]
        edge [arrowsize=0.8]
    
        subgraph cluster_peak {
          label = 'Peak Model\n(5 Parameters)';
          fontname = 'Helvetica-Bold';
          style = rounded;
          color = black;
          bgcolor = 'white';
          
          node [shape = ellipse, style = filled, fillcolor = white]
          a [label = 'a\n~\nLogNormal']
          mu [label = 'μ\n~\nNormal']
          rho [label = 'ρ\n~\nLogNormal']
          baseline_p [label = 'baseline\n~\nLogNormal']
          eps_p [label = 'ε\n~\nLogNormal']
    
          subgraph cluster_N1 {
            label = 'N'; labeljust = r; labelloc = b; style = rounded;
            
            node [shape = box, style = filled, fillcolor = white]
            y_hat_p [label = 'ŷ\n~\nGaussian(a, μ, ρ)']
            
            node [shape = ellipse, style = filled, fillcolor = 'grey85']
            y_p [label = 'y\n~\nNormal']
            
            y_hat_p -> y_p
          }
          
          {a, mu, rho, baseline_p} -> y_hat_p
          eps_p -> y_p
        }
    
        subgraph cluster_null {
          label = 'Noise Model\n(2 Parameters)';
          fontname = 'Helvetica-Bold';
          style = rounded;
          color = black;
          
          node [shape = ellipse, style = filled, fillcolor = white]
          baseline_n [label = 'baseline\n~\nLogNormal']
          eps_n [label = 'ε\n~\nLogNormal']
    
          subgraph cluster_N2 {
            label = 'N'; labeljust = r; labelloc = b; style = rounded;
            
            node [shape = box, style = filled, fillcolor = white]
            y_hat_n [label = 'ŷ\n~\nConstant']
            
            node [shape = ellipse, style = filled, fillcolor = 'grey85']
            y_n [label = 'y\n~\nNormal']
            
            y_hat_n -> y_n
          }
    
          baseline_n -> y_hat_n
          eps_n -> y_n
        }
    
        node [shape = box, style = 'filled,rounded']
        BIC_calc [label = 'Bayesian Model Selection\nEvidence = BIC(Noise) - BIC(Peak)']
    
        # Invisible edges to enforce layout
        y_p -> BIC_calc [style=dashed, label='LogLik']
        y_n -> BIC_calc [style=dashed, label='LogLik']
      }
    ")
    
  return(graph)
}