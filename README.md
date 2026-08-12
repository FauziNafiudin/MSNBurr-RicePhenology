# MSNBurr-RicePhenology

**On The Bayesian Neo-Normal Mixture Model for Unsupervised Phenological Zonation using Sentinel-2 Time Series**

Fauzi Nafi'udin, Nur Iriawan*, Tintrim Dwi Ary Widhianingsih
Department of Statistics, Faculty of Science and Data Analytics,
Institut Teknologi Sepuluh Nopember (ITS), Surabaya, Indonesia
\* Corresponding author: `nur.iriawan@its.ac.id`

---

## 📖 Overview

This repository provides the **complete, fully reproducible implementation** of the
companion study published in *IEEE Journal of Selected Topics in Applied Earth
Observations and Remote Sensing (JSTARS)*.

The study proposes a **Bayesian Neo-Normal Mixture Model** built on the
**Modified Skew-Normal Burr (MSNBurr)** distribution for unsupervised probabilistic
zonation of paddy phenology. Applied to Sentinel-2 NDVI time series over Lamongan
Regency, East Java, Indonesia (November 2022 – February 2025; **11,371–15,058
phenological cycles**), the framework:

- **Jointly captures right-skewness and heavy tails** through a closed-form algebraic likelihood — no numerical integration of the Student-t CDF, unlike Skew-t formulations;
- **Decisively outperforms a structurally matched GMM** across three segmentation scenarios (Strict / Default / Lenient) and both K = 2 and K = 3 configurations, with PSIS-LOO-CV |Δelpd|/SE ratios of **12.99–24.49**;
- Delivers an **≈100-fold computational speed-up** per effective sample over the Skew-t baseline (2,000 HMC iterations in ~2 h vs. 100 iterations in ~10 h);
- **Natively isolates a water-stressed rice sub-population** at K = 3 (truncated plateau length, suppressed peak NDVI), independently corroborated by CHIRPS precipitation data (Wilcoxon rank-sum, *p* < 2.2×10⁻¹⁶) during the 2023–2024 El Niño event;
- Yields **full posterior soft-membership probabilities** for pixel-level uncertainty quantification.

---

## 🔄 Pipeline

```
┌────────────────────────── STAGE 1 · Python ──────────────────────────┐
│ 01 notebook:                                                         │
│   Sentinel-2 L2A (GEE + SCL masking)                                 │
│   → 5% SRS over ESA WorldCover cropland frame (5,635 candidates)     │
│   → daily linear interpolation → Savitzky–Golay smoothing            │
│   → two-stage sensitivity analysis → peak-based segmentation         │
│   → 5 agronomic features under Strict / Default / Lenient scenarios  │
└───────────────────────────────┬──────────────────────────────────────┘
                                ▼
┌────────────────────────── STAGE 2 · R / Stan ────────────────────────┐
│ 02  GMM + MSNBurr mixture models (HMC-NUTS, RStan)   [production]    │
│ 03  Skew-t benchmark (loads 02's saved posteriors)   [optional]      │
│ 04  posterior tables, diagnostics & combined PPC                     │
│ 05  CHIRPS extraction within each cycle's biological window          │
│ 06  proxy validation (Wilcoxon / Spearman) + cropping-calendar plots │
└──────────────────────────────────────────────────────────────────────┘
```

---

## 📁 Repository Structure

```
MSNBurr-RicePhenology/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── requirements.txt                 # pinned Python dependencies
├── environment.yml                  # optional conda environment
│
├── data/
│   ├── README.md                    # provenance, Zenodo DOI, CHIRPS download script
│   ├── raw/                         # Lamongan_NDVI_Raw_5pct_2022_2024.csv   [Zenodo]
│   ├── processed/                   # 001–006 CSVs from the notebook         [Zenodo]
│   ├── chirps/                      # chirps-v2.0.{2022..2025}.days_p05.nc   [download only]
│   └── shapefiles/                  # BPS admin-2 boundaries (HDX)           [committed]
│
├── gee/
│   └── extract_ndvi_lamongan.js     # GEE: Sentinel-2 L2A, SCL masking, 5% SRS
│
├── notebooks/
│   └── 01_Sentinel2_NDVI_Preprocessing_and_Cycle_Extraction_Lamongan.ipynb
│
├── r_analysis/
│   ├── stan/
│   │   ├── normal.stan              # GMM baseline
│   │   ├── msn_burr.stan            # MSNBurr (proposed)
│   │   └── skewt.stan               # Skew-t (computational benchmark)
│   ├── scripts/
│   │   ├── 02_run_gmm_msnburr.R     # production fits: GMM + MSNBurr, pruning,
│   │   │                            #   traceplots, PPC, clustering CSV
│   │   ├── 03_run_skewt_benchmark.R # Skew-t only (100 iterations); reuses 02's RDS
│   │   ├── 04_post_analysis.R       # convergence/LOO/CI tables + combined PPC
│   │   ├── 05_extract_chirps.R      # offline NetCDF + terra extraction per cycle
│   │   └── 06_proxy_validation.R    # Wilcoxon/Spearman statistics + Figs. 6–7
│   └── outputs/
│       ├── fits/                    # fit_{gmm,msnburr}_pruned_*.rds  [gitignored → Zenodo]
│       ├── traceplots/              # traceplot_{gmm,msnburr}_*.png
│       ├── ppc/                     # ppc_{gmm,msnburr}_*.png, ppc_combined_*.png
│       └── clustering/              # clustering_results_*.csv,
│                                    #   clustering_results_with_skewt_*.csv,
│                                    #   clustering_with_chirps_validation_FINAL.csv
│
├── figures/
│   ├── main/                        # Fig. 1–7
│   └── supplementary/               # Fig. S1–S7
│
├── results/
│   └── tables/                      # Tables I–III and S1–S3 as CSV
│
└── paper/
    ├── manuscript.pdf
    └── supplementary.pdf
```

---

## 🗂️ Data

| Dataset | Role | Source |
|---|---|---|
| Sentinel-2 Level-2A (B4, B8) | NDVI time series, Nov 2022 – Feb 2025 | ESA / Copernicus via Google Earth Engine |
| ESA WorldCover 10 m 2021 v200 | Spatial sampling frame **only** | [Zenodo · 10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221) |
| CHIRPS v2.0 daily rainfall | Independent physical validation | [Climate Hazards Center, UC Santa Barbara](https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/) |
| BPS administrative boundaries (adm-2) | Study-area cartography | HDX / OCHA |

**Sampling design.** The regency was discretized into a 100 m × 100 m grid and a
Simple Random Sample with 5% sampling fraction (mirroring the official BPS
fraction) was drawn over the open-cropland frame, yielding **5,635 candidate
locations**, of which **5,311** produced at least one valid phenological cycle
(324 persistently bare/fallow pixels excluded).

**Archiving policy.** Large artefacts are **not committed to Git**: the raw and
processed CSVs and the pruned posterior RDS files (≈40 MB each) are archived on
Zenodo (see `data/README.md` for the download script), and the CHIRPS NetCDF
files are re-downloadable from the official server. Small derived artefacts
(clustering CSVs, tables, figures) are committed.

---

## ⚙️ Installation

**Python (Stage 1)**

```bash
conda env create -f environment.yml    # or: pip install -r requirements.txt
conda activate msnburr
```

Core dependencies: `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`,
`geopandas`, `contextily`, `matplotlib-scalebar`, `tqdm`.

**R (Stage 2)**

```r
install.packages(c("rstan", "loo", "posterior", "BH", "dplyr", "ggplot2",
                   "patchwork", "e1071", "bayesplot", "aricode", "digest",
                   "showtext", "terra", "tidyr", "data.table", "gridExtra",
                   "lubridate", "cowplot"))
```

Requires R ≥ 4.3 with a working C++ toolchain (RStan system requirements).

---

## 🚀 Running the Pipeline

Run from the repository root. Each R script reads its configuration
(`SCENARIO`, `K_CLUSTERS`) from the clearly marked config block at its top;
all outputs are written to `r_analysis/outputs/`.

```bash
# Stage 1 — NDVI preprocessing & feature extraction (Python)
jupyter notebook notebooks/01_Sentinel2_NDVI_Preprocessing_and_Cycle_Extraction_Lamongan.ipynb

# Stage 2 — production fits; repeat for all 6 configurations
#           (SCENARIO ∈ {strict, default, lenient} × K ∈ {2, 3})
Rscript r_analysis/scripts/02_run_gmm_msnburr.R

# Stage 3 — optional Skew-t computational benchmark (requires 02's RDS)
Rscript r_analysis/scripts/03_run_skewt_benchmark.R

# Stage 4 — posterior tables, diagnostics & combined PPC (Fig. 5 / Fig. S7)
Rscript r_analysis/scripts/04_post_analysis.R

# Stage 5 — CHIRPS extraction for the K=3 Lenient configuration
#           (requires data/chirps/*.nc; see data/README.md)
Rscript r_analysis/scripts/05_extract_chirps.R

# Stage 6 — proxy validation statistics + Figs. 6–7
Rscript r_analysis/scripts/06_proxy_validation.R
```

---

## 🖼️ Figure & Table Provenance

| Paper artefact | Produced by | Output file |
|---|---|---|
| Fig. 1 (study area) | `01` notebook | `lamongan_study_area.png` |
| Fig. 2 (sampled locations) | `01` notebook | `lamongan_points_only.png` |
| Fig. 3 (segmentation example) | `01` notebook | `fig_3_segmentation_example.png` |
| Fig. 4 (marginal distributions) | `01` notebook | `fig_4_marginal_distributions.png` |
| Fig. 5 (combined PPC, K=3) | `04_post_analysis.R` | `ppc_combined_*_K3_*.png` |
| Fig. 6 (proxy validation) | `06_proxy_validation.R` | `proxy_validation_IEEE.png` |
| Fig. 7 (SOS/EOS dynamics) | `06_proxy_validation.R` | `combined_trend_wide_IEEE.png` |
| Fig. S1–S6 | `01` notebook | `fig_S1…S6_*.png` |
| Fig. S7 (combined PPC, K=2) | `04_post_analysis.R` | `ppc_combined_*_K2_*.png` |
| Table I (HMC-NUTS config) | documented in scripts & this README | — |
| Table II (K=3 agronomic estimates) | `04_post_analysis.R` | `results/tables/` |
| Table III (rainfall comparison) | `06_proxy_validation.R` | `results/tables/` |
| Table S1 (extraction parameters) | `01` notebook config | `results/tables/` |
| Table S2 (full model comparison) | `04_post_analysis.R` (all configs) | `results/tables/` |
| Table S3 (K=2 estimates) | `04_post_analysis.R` | `results/tables/` |

---

## ✅ Reproducibility

| Item | Value |
|---|---|
| Python random seed | `42` |
| R / Stan global seed | `12345` |
| HMC-NUTS configuration | 4 chains × 2,000 iterations (1,000 warmup), storage thinning 5, `adapt_delta = 0.90`, `max_treedepth = 12` |
| Skew-t benchmark | 4 chains × 100 iterations (50 warmup), `adapt_delta = 0.85` |
| Initialization | K-means (n_start = 25), centroids sorted by anchor (cycle duration); MSNBurr α = 1.0 |
| Identifiability | Ordered constraint on the anchor feature only: μ₁ < μ₂ < ⋯ < μ_K |
| Convergence gate | R̂ < 1.01, sufficient ESS, zero divergent transitions |
| Hardware | 11th Gen Intel Core i7-11700F (2.50 GHz), 16 GB DDR4 |

**Expected Stage-1 outputs** (validated automatically by the notebook):

| Scenario | Cycles | Locations |
|---|---|---|
| Strict | 11,371 | 4,896 |
| Default | 13,448 | 5,125 |
| Lenient | 15,058 | 5,311 |

---

## 📊 Headline Results

| K | Scenario | \|Δelpd\|/SE (MSNBurr vs GMM) |
|---|---|---|
| 2 | Strict / Default / Lenient | 13.13 / 17.90 / 24.04 |
| 3 | Strict / Default / Lenient | 12.99 / 17.50 / 24.49 |

All ratios far exceed the decisive-evidence threshold of 10.
ARI(GMM, MSNBurr): 0.73–0.76 (K = 2) and 0.36–0.58 (K = 3), reflecting the
MSNBurr's unique capacity to resolve the intermediate, asymmetric sub-population.

---

## 📚 Citation

```bibtex
@article{nafiudin2026msnburr,
  title     = {On The Bayesian Neo-Normal Mixture Model for Unsupervised
               Phenological Zonation using Sentinel-2 Time Series},
  author    = {Nafi'udin, Fauzi and Iriawan, Nur and Widhianingsih, Tintrim Dwi Ary},
  journal   = {IEEE Journal of Selected Topics in Applied Earth Observations
               and Remote Sensing},
  year      = {2026},
  publisher = {IEEE}
}
```

---

## 🙏 Acknowledgments

Supported by the Ministry of Higher Education, Science, and Technology of the
Republic of Indonesia through BIMA under the Magister Thesis Research Grant
(Program Penelitian Baru) Tahun 2026, Master Contract No.
171/C3/DT.05.00/PL-BARU/2026 and Derivative Contract No. 2168/PKS/ITS/2026.
The authors thank ESA (Sentinel-2 & WorldCover), the Climate Hazards Center
(CHIRPS), and Google (Earth Engine) for data and platform access.

## 📄 License

MIT — see [`LICENSE`](LICENSE).

## 📮 Contact

- Fauzi Nafi'udin — `6003251024@student.its.ac.id`
- Nur Iriawan — `nur.iriawan@its.ac.id`
