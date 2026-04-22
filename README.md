# Probabilistic Signal Extraction on Sparse Sensor Data

[![Live App](https://img.shields.io/badge/Shiny-Live%20App-blue?logo=r)](https://j-lamadrid.shinyapps.io/probabilistic-sensor-uq/)
[![R](https://img.shields.io/badge/R-4.5-276DC3?logo=r)](https://www.r-project.org/)
[![Stan](https://img.shields.io/badge/Stan-Variational%20Bayes-B2001D)](https://mc-stan.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Bayesian peak detection and noise estimation for LC-MS spectra, with an interactive 3D visualization app deployed on shinyapps.io.**

---

## Overview

Peak and noise estimation is a persistent challenge in scientific signal analysis. Confidently identifying and characterizing points of interest in highly complex, high-dimensional signals is essential for uncertainty quantification in sensor data — yet standard methods often lack rigorous probabilistic foundations or become computationally intractable at scale.

This project proposes a **Bayesian framework for peak and noise estimation** with temporal considerations, enabling comprehensive aggregate analysis and sensor validation across full LC-MS recordings. The approach draws on techniques from [Tokuda et al.](https://doi.org/10.1371/journal.pone.0197440), adapting their framework for estimating noise variance and peak count via Bayesian model selection and posterior density estimation to the practical constraints of large, complex spectra.

For spectra where full MCMC is computationally infeasible, the framework instead uses **Variational Bayes (VB)** with strongly informed priors, enabling tractable inference at scale without sacrificing the probabilistic interpretability of the results.

---

## Methods

### Model Selection via BIC

Each candidate peak region is evaluated under two competing models:

- **Peak model** — a Gaussian signal with 5 free parameters: amplitude $a$, position $\mu$, precision $\rho$, baseline, and noise $\varepsilon$
- **Noise model** — a flat baseline with 2 free parameters: baseline and noise $\varepsilon$

Model selection uses the **BIC difference** as an approximation to the log Bayes factor:

$$\text{Evidence} = \text{BIC}_{\text{noise}} - \text{BIC}_{\text{peak}}$$

A threshold of $> 2$ is used to accept a peak as statistically supported over noise.

### Variational Bayes Inference

Rather than full MCMC sampling, posterior densities are approximated using Stan's **Automatic Differentiation Variational Inference (ADVI)**. This reduces per-peak inference from minutes to milliseconds, making whole-run analysis of hundreds of spectra tractable.

### Temporal Aggregation

Detected peaks are grouped across spectra by m/z proximity (tolerance: 0.05 Da) and summarized as aggregate features with:

- Mean and 95% HDI for amplitude, position, precision, noise variance, and Bayes factor evidence
- Retention time range (first and last appearance)
- Spectrum count (used to filter sensor artifacts: < 1% or > 80% of total spectra)

---

## Peak Model

The generative model for a single peak region:

$$\hat{y}_i = a \cdot \exp\left(-\frac{\rho}{2}(x_i - \mu)^2\right) + \text{baseline}$$

$$y_i \sim \mathcal{N}(\hat{y}_i,\ \varepsilon^2)$$

With priors:

| Parameter | Prior | Description |
|---|---|---|
| $a$ | LogNormal | Peak amplitude |
| $\mu$ | Normal | Peak m/z position |
| $\rho$ | LogNormal | Peak precision (inverse width) |
| baseline | LogNormal | Spectral baseline |
| $\varepsilon$ | LogNormal | Noise standard deviation |

---

## Interactive App

A live Shiny app is deployed for interactive exploration of processed results:

**[Launch App →](https://j-lamadrid.shinyapps.io/probabilistic-sensor-uq/)**

Features:
- **2D fingerprint plot** — collapsed mass spectrum showing reconstructed peak profiles across m/z
- **3D temporal plot** — waterfall view of intensity across m/z and retention time
- **Aggregated features table** — per-m/z-group summary statistics matching the `_agg.rds` output
- Upload your own `.mzML` file (up to 1 GB) or select from pre-processed cultivar results

> A demo file (`data/raw/demo/demo_lcms.mzML`) with 6 simulated compounds across an 8-minute gradient is included for testing.

---

## Data

LC-MS recordings from [Metabolomics Workbench Study ST003061](https://www.metabolomicsworkbench.org/), comprising three sweet potato cultivars:

| File | Cultivar |
|---|---|
| `BTY-1.mzML` | Bantianyao |
| `FY6-1.mzML` | FY6 |
| `JFH-1.mzML` | Jifenghong |

Pre-processed aggregate results (`.rds` files) are included in `data/processed/`.

---

## Repository Structure

```
.
├── app.R                        # Shiny app entry point
├── R/
│   ├── utils.R                  # mzML loading, binning, peak finding
│   └── plot_3d.R                # 2D/3D plot generation
├── scripts/
│   ├── model.R                  # Offline Stan modeling pipeline
│   └── plot_utils.R             # Diagnostic plots (PPC, pairs, evidence)
├── models/
│   ├── peak_model.stan          # Gaussian peak Stan model
│   └── noise_model.stan         # Baseline noise Stan model
├── data/
│   ├── processed/               # Aggregated per-cultivar RDS results
│   └── raw/                     # Raw mzML files and demo data
└── figures/                     # Output figures
```

---

## Dependencies

| Package | Role |
|---|---|
| `rstan` | Variational Bayes inference (offline modeling) |
| `mzR` | mzML file parsing (onDisk, memory-efficient) |
| `pracma` | Peak finding (`findpeaks`) |
| `shiny` + `plotly` | Interactive app and 3D visualization |
| `dplyr` | Data wrangling |
| `HDInterval` | 95% HDI computation |
| `e1071` | Spectral skewness/kurtosis features |

---

## References

S. Tokuda, et al., "Simultaneous Estimation of Noise Variance and Number of Peaks in Bayesian Spectral Deconvolution," *Journal of the Physical Society of Japan*, vol. 86, no. 2, p. 024001, 2016. DOI: [10.7566/jpsj.86.024001](https://doi.org/10.7566/jpsj.86.024001)

Metabolomics Workbench, Project ID PR001907. Available at [metabolomicsworkbench.org](https://www.metabolomicsworkbench.org). DOI: [10.21228/M8714X](https://doi.org/10.21228/M8714X)

V. Mazet, et al., "Unsupervised Joint Decomposition of a Spectroscopic Signal Sequence," *Signal Processing*, vol. 109, pp. 193–205, 2014. DOI: [10.1016/j.sigpro.2014.10.032](https://doi.org/10.1016/j.sigpro.2014.10.032)

G. C. Allen and R. F. McMeeking, "Deconvolution of Spectra by Least-Squares Fitting," *Analytica Chimica Acta*, vol. 103, no. 1, pp. 73–108, 1978.

H. Sayama, "Mean-Field Approximation," in *Introduction to the Modeling and Analysis of Complex Systems*. Available at [math.libretexts.org](https://math.libretexts.org/)

P. Miketova, et al., "Tandem Mass Spectrometry Studies of Green Tea Catechins," *Journal of Mass Spectrometry*, vol. 35, no. 7, pp. 860–869, 2000.

W. Shao and H. Lam, "Denoising Peptide Tandem Mass Spectra for Spectral Libraries: A Bayesian Approach," *Journal of Proteome Research*, vol. 12, no. 7, pp. 3223–3232, 2013. DOI: [10.1021/pr400080b](https://doi.org/10.1021/pr400080b)

---

## Applications

While applied here to LC-MS metabolomics, the framework generalizes to any domain requiring principled peak detection in noisy, high-dimensional signals — including neuroscience (LFP/EEG), cosmology (spectral line detection), and engineering sensor validation.
