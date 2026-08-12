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

## 🔄 Pipeline

```
┌────────────────────────── STAGE 1 · Python ──────────────────────────┐
│ Sentinel-2 L2A (GEE + SCL masking)                                   │
│   → 5% SRS over ESA WorldCover cropland frame (5,635 candidates)     │
│   → daily linear interpolation → Savitzky–Golay smoothing            │
│   → peak-based segmentation (3 scenarios)                            │
│   → 5 agronomic features: duration, peak NDVI, slope up,             │
│     slope down, plateau length (≥95% of peak)                        │
└───────────────────────────────┬──────────────────────────────────────┘
                                ▼
┌────────────────────────── STAGE 2 · R / Stan ────────────────────────┐
│ Bayesian mixture models: MSNBurr vs GMM vs Skew-t (HMC-NUTS, RStan)  │
│   → PSIS-LOO-CV · PPC · ARI · Shannon entropy                        │
│   → K=3 sub-population isolation                                     │
│   → independent CHIRPS rainfall validation (Wilcoxon rank-sum)       │
└──────────────────────────────────────────────────────────────────────┘
```

## 📁 Repository Structure

```
MSNBurr-RicePhenology/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── requirements.txt               # pinned Python dependencies
├── environment.yml                # optional conda environment
│
├── data/
│   ├── README.md                  # data provenance & download instructions
│   ├── raw/                       # Lamongan_NDVI_Raw_5pct_2022_2024.csv
│   ├── processed/                 # 001_*.csv … 006_final_features_{scenario}.csv
│   └── shapefiles/                # BPS admin-2 boundaries (HDX)
│
├── gee/
│   └── extract_ndvi_lamongan.js   # GEE script: Sentinel-2 L2A, SCL masking, 5% SRS
│
├── notebooks/
│   └── 01_Sentinel2_NDVI_Preprocessing_and_Cycle_Extraction_Lamongan.ipynb
│
├── r_analysis/
│   ├── 02_Bayesian_Mixture_Model.Rmd
│   ├── stan/
│   │   ├── msnburr_mixture.stan
│   │   ├── gmm_mixture.stan
│   │   └── skewt_mixture.stan
│   ├── R/
│   │   ├── utils_data.R
│   │   ├── utils_diagnostics.R    # R-hat, ESS, PPC, ARI, entropy
│   │   └── utils_loo.R            # PSIS-LOO-CV comparison
│   └── scripts/
│       ├── run_all_models.R       # HMC-NUTS configuration (Table I)
│       └── chirps_validation.R    # CHIRPS extraction + Wilcoxon (Table III, Figs. 6–7)
│
├── figures/
│   ├── main/                      # Fig. 1–7
│   └── supplementary/             # Fig. S1–S7
│
├── results/
│   ├── tables/                    # Tables I–III and S1–S3 as CSV
│   └── posteriors/                # thinned posterior draws
│
└── paper/
    ├── manuscript.pdf
    └── supplementary.pdf
```

## 🗂️ Data

| Dataset | Role | Source |
|---|---|---|
| Sentinel-2 Level-2A (B4, B8) | NDVI time series, Nov 2022 – Feb 2025 | ESA / Copernicus via Google Earth Engine |
| ESA WorldCover 10 m 2021 v200 | Spatial sampling frame **only** | [Zenodo · 10.5281/zenodo.7254221](https://doi.org/10.5281/zenodo.7254221) |
| CHIRPS v2.0 daily rainfall | Independent physical validation | Climate Hazards Center, UC Santa Barbara |
| BPS administrative boundaries (adm-2) | Study-area cartography | HDX / OCHA |

Sampling design: the regency was discretized into a 100 m × 100 m grid and a
**Simple Random Sample with 5% sampling fraction** (mirroring the official BPS
fraction) was drawn over the open-cropland frame, yielding **5,635 candidate
locations**, of which **5,311** produced at least one valid phenological cycle
(324 persistently bare/fallow pixels excluded).

> **Note:** raw and processed CSVs are archived on Zenodo due to size; see
> `data/README.md` for the download script.

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
install.packages(c("rstan", "loo", "posterior", "mclust", "ggplot2", "sf", "terra"))
```

Requires R ≥ 4.3 with a working C++ toolchain (RStan system requirements).

## 🚀 Running the Pipeline

```bash
# Stage 1 — NDVI preprocessing & phenological feature extraction (Python)
jupyter notebook notebooks/01_Sentinel2_NDVI_Preprocessing_and_Cycle_Extraction_Lamongan.ipynb

# Stage 2 — Bayesian mixture modelling & validation (R)
Rscript r_analysis/scripts/run_all_models.R
Rscript r_analysis/scripts/chirps_validation.R
```

## ✅ Reproducibility

| Item | Value |
|---|---|
| Python random seed | `42` |
| R / Stan global seed | `12345` |
| HMC-NUTS configuration | 4 chains × 2,000 iterations (1,000 warmup), storage thinning 5, `adapt_delta = 0.90`, `max_treedepth = 12` |
| Initialization | K-means (n_start = 25), centroids sorted by the anchor feature (cycle duration); MSNBurr α initialized at 1.0 |
| Identifiability | Ordered constraint on the anchor feature only: μ₁ < μ₂ < ⋯ < μ_K |
| Convergence gate | R̂ < 1.01, sufficient ESS, zero divergent transitions |
| Hardware | 11th Gen Intel Core i7-11700F (2.50 GHz), 16 GB DDR4 |

**Expected Stage-1 outputs** (the notebook validates these automatically):

| Scenario | Cycles | Locations |
|---|---|---|
| Strict | 11,371 | 4,896 |
| Default | 13,448 | 5,125 |
| Lenient | 15,058 | 5,311 |

## 📊 Headline Results

| K | Scenario | \|Δelpd\|/SE (MSNBurr vs GMM) |
|---|---|---|
| 2 | Strict / Default / Lenient | 13.13 / 17.90 / 24.04 |
| 3 | Strict / Default / Lenient | 12.99 / 17.50 / 24.49 |

All ratios far exceed the decisive-evidence threshold of 10.
ARI(GMM, MSNBurr): 0.73–0.76 (K = 2) and 0.36–0.58 (K = 3), reflecting the
MSNBurr's unique capacity to resolve the intermediate, asymmetric sub-population.

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
