library(MSnbase) 
library(mzR) 
library(BiocParallel) 
library(e1071) 
library(progressr) 
library(dplyr) 
library(tidyr) 

register(SerialParam()) 



loadMZML <- function(fp, cultivar) {
  ms_data = readMSData(fp, mode = 'inMemory')
  print(ms_data)
  
  extract_features <- function(ms_data, group_label) {
    intensities = intensity(ms_data)
    mz_vals = mz(ms_data)
    rts = rtime(ms_data)
    
    df = data.frame(
      spectrum_id = seq_along(intensities),
      n_peaks = lengths(intensities),
      total_ion_current = sapply(intensities, sum, na.rm = TRUE),
      retention_time = rts,
      group = group_label,
      intensity_skewness = sapply(intensities, function(x) {
        if (length(x) > 2)
          e1071::skewness(x, na.rm = TRUE)
        else
          NA
      }),
      intensity_kurtosis = sapply(intensities, function(x) {
        if (length(x) > 3)
          e1071::kurtosis(x, na.rm = TRUE)
        else
          NA
      })
    )
    
    df$intensities = intensities
    df$mz = mz_vals
    return(df)
  }
  
  spectra_df = extract_features(ms_data, cultivar)
  
  return(spectra_df)
}



findRegions <- function(binned_spectrum,
                        window_size = 3,
                        n_peaks = 5) {
  valid_intensities = binned_spectrum$intensity[binned_spectrum$intensity > 0]
  noise_baseline = median(valid_intensities)
  noise_mad = mad(valid_intensities) # median absolute deviation
  min_peak_height = noise_baseline + (3 * noise_mad)
  
  # pracma
  peaks = findpeaks(
    binned_spectrum$intensity,
    minpeakheight = min_peak_height,
    minpeakdistance = 5,
    npeaks = n_peaks,
    sortstr = TRUE
  )
  
  if (is.null(peaks)) {
    cat("No peaks found. Try lower min_peak_height.\n")
    return(NULL)
  }
  
  regions = list()
  for (i in 1:nrow(peaks)) {
    peak_mz = binned_spectrum$mz[peaks[i, 2]]
    peak_height = peaks[i, 1]
    
    mz_start = peak_mz - window_size / 2
    mz_end = peak_mz + window_size / 2
    
    region_idx = which(binned_spectrum$mz >= mz_start &
                         binned_spectrum$mz <= mz_end)
    region_int = binned_spectrum$intensity[region_idx]
    
    peak_to_noise = peak_height / median(region_int[region_int > 0])
    n_points = sum(region_int > 0)
    
    regions[[i]] = list(
      peak_mz = peak_mz,
      # mass-to-charge
      peak_height = peak_height,
      # amplitude
      mz_range = c(mz_start, mz_end),
      n_points = n_points,
      peak_to_noise = peak_to_noise,
      # snr
      quality_score = peak_to_noise * sqrt(n_points)
    )
  }
  
  quality_scores = sapply(regions, function(r)
    r$quality_score)
  regions = regions[order(quality_scores, decreasing = TRUE)]
  
  return(regions)
}



binData <- function(features, idx, bin_width = 0.015) {
  mz_vals = features$mz[[idx]]
  intensity_vals = features$intensities[[idx]]
  
  keep_idx = intensity_vals > 0
  if (sum(keep_idx) == 0) {
    return(data.frame(mz = numeric(0), intensity = numeric(0)))
  }
  
  mz_vals = mz_vals[keep_idx]
  intensity_vals = intensity_vals[keep_idx]
  
  mz_binned = round(mz_vals / bin_width) * bin_width
  
  temp_df = data.frame(mz = mz_binned, intensity = intensity_vals)
  
  binned_spectra = aggregate(intensity ~ mz, data = temp_df, FUN = max)
  
  return(binned_spectra)
}



groupPeaks <- function(mz_values, tolerance = 0.05) {
  ord = order(mz_values)
  sorted_mz = mz_values[ord]
  
  gaps = c(TRUE, diff(sorted_mz) > tolerance)
  cluster_ids = cumsum(gaps)
  
  final_ids = integer(length(mz_values))
  final_ids[ord] = cluster_ids
  
  return(final_ids)
}