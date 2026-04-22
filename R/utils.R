# Libraries are loaded in app.R. This file defines helper functions only.
# (mzR, e1071, dplyr, tidyr must be loaded before these functions are called)

# loadMZML: uses mzR directly (no MSnbase inMemory) so RAM stays flat.
# mzR opens the file as a handle and reads spectra one at a time on demand.
# The returned data frame stores mz/intensity as list-columns exactly as
# the rest of the app expects, but the raw MSnbase object is never created.
loadMZML <- function(fp, cultivar) {

  if (!file.exists(fp))
    stop(paste("File not found:", fp))

  handle  <- mzR::openMSfile(fp)
  on.exit(mzR::close(handle), add = TRUE)

  hdr <- mzR::header(handle)

  # keep only MS1
  ms1_idx <- which(hdr$msLevel == 1)
  if (length(ms1_idx) == 0) {
    # fall back to all scans if msLevel column absent or all NA
    ms1_idx <- seq_len(nrow(hdr))
  }

  hdr_ms1 <- hdr[ms1_idx, ]
  n        <- nrow(hdr_ms1)

  # pre-allocate output vectors / lists
  intensities      <- vector("list", n)
  mz_vals          <- vector("list", n)
  retention_times  <- numeric(n)
  n_peaks_vec      <- integer(n)
  tic_vec          <- numeric(n)
  skew_vec         <- numeric(n)
  kurt_vec         <- numeric(n)

  for (i in seq_len(n)) {
    raw_idx <- ms1_idx[i]
    peaks   <- mzR::peaks(handle, raw_idx)   # matrix: col1=mz, col2=intensity

    mz_i  <- peaks[, 1]
    int_i <- peaks[, 2]

    mz_vals[[i]]         <- mz_i
    intensities[[i]]     <- int_i
    retention_times[i]   <- hdr_ms1$retentionTime[i] / 60  # convert s -> min
    n_peaks_vec[i]       <- length(int_i)
    tic_vec[i]           <- sum(int_i, na.rm = TRUE)

    skew_vec[i] <- if (length(int_i) > 2)
      tryCatch(e1071::skewness(int_i, na.rm = TRUE), error = function(e) NA)
      else NA

    kurt_vec[i] <- if (length(int_i) > 3)
      tryCatch(e1071::kurtosis(int_i, na.rm = TRUE), error = function(e) NA)
      else NA
  }

  df <- data.frame(
    spectrum_id        = seq_len(n),
    n_peaks            = n_peaks_vec,
    total_ion_current  = tic_vec,
    retention_time     = retention_times,
    group              = cultivar,
    intensity_skewness = skew_vec,
    intensity_kurtosis = kurt_vec,
    stringsAsFactors   = FALSE
  )

  df$intensities <- intensities
  df$mz          <- mz_vals

  gc()   # release any transient allocations
  return(df)
}



findRegions <- function(binned_spectrum,
                        window_size = 3,
                        n_peaks = 5) {
  valid_intensities = binned_spectrum$intensity[binned_spectrum$intensity > 0]
  noise_baseline = median(valid_intensities)
  noise_mad = mad(valid_intensities)
  min_peak_height = noise_baseline + (3 * noise_mad)

  peaks = pracma::findpeaks(
    binned_spectrum$intensity,
    minpeakheight = min_peak_height,
    minpeakdistance = 5,
    npeaks = n_peaks,
    sortstr = TRUE
  )

  if (is.null(peaks)) {
    return(NULL)
  }

  regions = list()
  for (i in 1:nrow(peaks)) {
    peak_mz    = binned_spectrum$mz[peaks[i, 2]]
    peak_height = peaks[i, 1]

    mz_start = peak_mz - window_size / 2
    mz_end   = peak_mz + window_size / 2

    region_idx = which(binned_spectrum$mz >= mz_start &
                       binned_spectrum$mz <= mz_end)
    region_int = binned_spectrum$intensity[region_idx]

    peak_to_noise = peak_height / median(region_int[region_int > 0])
    n_points      = sum(region_int > 0)

    regions[[i]] = list(
      peak_mz       = peak_mz,
      peak_height   = peak_height,
      mz_range      = c(mz_start, mz_end),
      n_points      = n_points,
      peak_to_noise = peak_to_noise,
      quality_score = peak_to_noise * sqrt(n_points)
    )
  }

  quality_scores = sapply(regions, function(r) r$quality_score)
  regions        = regions[order(quality_scores, decreasing = TRUE)]
  return(regions)
}



binData <- function(features, idx, bin_width = 0.015) {
  mz_vals       = features$mz[[idx]]
  intensity_vals = features$intensities[[idx]]

  keep_idx = intensity_vals > 0
  if (sum(keep_idx) == 0)
    return(data.frame(mz = numeric(0), intensity = numeric(0)))

  mz_vals       = mz_vals[keep_idx]
  intensity_vals = intensity_vals[keep_idx]

  mz_binned = round(mz_vals / bin_width) * bin_width

  temp_df = data.frame(mz = mz_binned, intensity = intensity_vals)
  binned_spectra = aggregate(intensity ~ mz, data = temp_df, FUN = max)
  return(binned_spectra)
}



groupPeaks <- function(mz_values, tolerance = 0.05) {
  ord       = order(mz_values)
  sorted_mz = mz_values[ord]

  gaps      = c(TRUE, diff(sorted_mz) > tolerance)
  cluster_ids = cumsum(gaps)

  final_ids        = integer(length(mz_values))
  final_ids[ord]   = cluster_ids
  return(final_ids)
}
