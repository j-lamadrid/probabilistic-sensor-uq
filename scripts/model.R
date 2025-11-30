library(rstan)
library(dplyr)
library(ggplot2)
library(pracma)
library(here)
library(HDInterval)

source(here("R", "utils.R"))



CULTIVAR = "Bantianyao"
FILEPATH = "C:/Users/jlama/Documents/probabilistic-sensor-uq/data/raw/ST003061/BTY-1.mzML"  # adjust as needed

cat(sprintf("Loading data for %s from %s...\n", CULTIVAR, FILEPATH))
if (!file.exists(FILEPATH)) stop("File not found!")
features = loadMZML(FILEPATH, CULTIVAR)

cat("Compiling models...\n")
model_peak = stan_model(file = here("models", "peak_model.stan"))
model_noise = stan_model(file = here("models", "noise_model.stan"))

all_peak_results = list()
N_SPECTRA = nrow(features)
N_SPECTRA_TO_PROCESS = N_SPECTRA # modify for testing

cat(sprintf("Starting analysis on %d spectra...\n", N_SPECTRA_TO_PROCESS))


for (i in 1:N_SPECTRA_TO_PROCESS) {
  
  binned_spectra = binData(features, i)
  
  if (nrow(binned_spectra) == 0) next
  
  regions = findRegions(binned_spectra, window_size = 0.8)
  
  if (is.null(regions)) next
  
  
  for (j in 1:length(regions)) {
    
    current_region = regions[[j]]
    
    spectrum = binned_spectra %>%
      filter(mz >= current_region$mz_range[1],
             mz <= current_region$mz_range[2])
    
    if (nrow(spectrum) < 3) next
    
    stan_data = list(
      N = nrow(spectrum),
      K = 1,
      x = as.array(spectrum$mz),
      y = as.array(spectrum$intensity),
      x_min = min(spectrum$mz),
      x_max = max(spectrum$mz),
      peak_mz_prior = current_region$peak_mz,
      peak_height_prior = current_region$peak_height,
      baseline_prior = median(spectrum$intensity),
      noise_prior = sd(spectrum$intensity)
    )
    
    bayes_result = tryCatch({
      
      init_peak = function() list(a = current_region$peak_height, 
                                  mu = current_region$peak_mz, 
                                  baseline = median(spectrum$intensity), 
                                  epsilon = sd(spectrum$intensity), 
                                  rho = 1000.0)
      
      fit_peak = suppressWarnings(vb(model_peak, 
                                     data = stan_data, 
                                     init = init_peak, 
                                     seed = 123, 
                                     refresh = 0, 
                                     tol_rel_obj = 0.005))
      
      init_noise = function() list(baseline = median(spectrum$intensity), 
                                   epsilon = sd(spectrum$intensity))
      
      fit_noise = suppressWarnings(vb(model_noise, 
                                      data = stan_data, 
                                      init = init_noise, 
                                      seed = 123, 
                                      refresh = 0, 
                                      tol_rel_obj = 0.005))
      
      res_peak = as.data.frame(fit_peak)
      res_noise = as.data.frame(fit_noise)
      
      bic_peak = mean(res_peak$bic)
      bic_noise = mean(res_noise$bic)
      
      evidence = bic_noise - bic_peak
      
      list(status = "success", fit = fit_peak, evidence = evidence)
      
    }, error = function(e) return(list(status = "error")))
    
    if (bayes_result$status == "success") {
      
      if (bayes_result$evidence > 2) {
        
        posterior = as.data.frame(bayes_result$fit)
        
        all_peak_results[[length(all_peak_results) + 1]] = data.frame(
          spectrum_id = i,
          cultivar = CULTIVAR,
          region_idx = j,
          peak_mz_prior = current_region$peak_mz,
          
          bayes_factor_evidence = bayes_result$evidence,
          noise_variance = mean(posterior$epsilon)^2,
          
          model_r2 = mean(posterior$R_squared),
          model_amplitude_a = mean(posterior$a),
          model_position_mu = mean(posterior$mu),
          model_precision_rho = mean(posterior$rho),
          model_baseline = mean(posterior$baseline)
        )
      }
    }
  }
  
  if (i %% 100 == 0) cat(sprintf("Processed Spectrum %d / %d\n", 
                                 i, 
                                 N_SPECTRA_TO_PROCESS))
}


if (length(all_peak_results) > 0) {
  
  final_table = bind_rows(all_peak_results)
  
  base_filename = tools::file_path_sans_ext(basename(FILEPATH))
  saveRDS(final_table, here("data", "processed", paste0(base_filename, ".rds")))
  
  MZ_TOLERANCE = 0.05
  final_table = final_table %>%
    mutate(peak_group_id = groupPeaks(peak_mz_prior, MZ_TOLERANCE))
  
  # average + HDI for evidence, noise, amplitude, position, precision 
  final_agg_table = final_table %>%
    group_by(cultivar, peak_group_id) %>%
    summarize(
      n_spectra = n(),
      
      mz_center = mean(peak_mz_prior), 
      mz_min = min(peak_mz_prior),
      mz_max = max(peak_mz_prior),
      
      avg_evidence = mean(bayes_factor_evidence),
      evidence_lower = hdi(bayes_factor_evidence, credMass = 0.95)["lower"],
      evidence_upper = hdi(bayes_factor_evidence, credMass = 0.95)["upper"],
      
      avg_noise_var = mean(noise_variance),
      noise_var_lower = hdi(noise_variance, credMass = 0.95)["lower"],
      noise_var_upper = hdi(noise_variance, credMass = 0.95)["upper"],
      
      mz_amplitude = mean(model_amplitude_a),
      amp_lower = hdi(model_amplitude_a, credMass = 0.95)["lower"],
      amp_upper = hdi(model_amplitude_a, credMass = 0.95)["upper"],
      
      mz_position = mean(model_position_mu),
      pos_lower = hdi(model_position_mu, credMass = 0.95)["lower"],
      pos_upper = hdi(model_position_mu, credMass = 0.95)["upper"],
      
      mz_precision = mean(model_precision_rho),
      prec_lower = hdi(model_precision_rho, credMass = 0.95)["lower"],
      prec_upper = hdi(model_precision_rho, credMass = 0.95)["upper"],
      
      mz_r2 = mean(model_r2),
      .groups = 'drop'
    )

  N_SPECTRA_THRESHOLD = N_SPECTRA_TO_PROCESS * 0.01 # remove sensor artifacts
  N_SPECTRA_MAX = N_SPECTRA_TO_PROCESS * 0.8 # remove sensor noise
  
  final_agg_table = final_agg_table %>%
    filter(n_spectra > N_SPECTRA_THRESHOLD) %>%
    filter(n_spectra < N_SPECTRA_MAX) %>%
    arrange(cultivar, peak_group_id)

  saveRDS(final_agg_table, here("data", "processed", 
                                paste0(base_filename, "_agg.rds")))

} else {
  cat("No peaks found.\n")
}