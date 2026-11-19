# Odd Lindley-G Beta Prime Distribution: Properties and Applications

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language: R](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)
[![Language: Python](https://img.shields.io/badge/Language-Python-3776AB.svg)](https://www.python.org/)

Replication code and datasets for the manuscript:

**Odd Lindley-G Beta prime Distribution: Properties and Applications**  
*Authors:* C.S. Rajitha and Cliya Symon

---

## Project Summary

This repository contains the code and data used to evaluate the **Odd Lindley-G Beta Prime (OL-GBP)** distribution and compare it against competing models on three real-world datasets ($n = 33$ each) from Indian states/territories:

1. **Per Capita Income** (Section 8.1.1, Tables 4–6, Figure 9)
2. **GDP** (Section 8.1.2, Tables 7–9, Figure 10)
3. **Birth Rate** (Section 8.1.3, Tables 10–12, Figure 11)

### Models Included
* **OL-GBP**: Odd Lindley-G Beta Prime distribution ($a, \alpha, \beta$)
* **WBP**: Weibull-Beta Prime distribution ($a, b, \alpha, \beta$)
* **BP**: Beta Prime distribution ($\alpha, \beta$)
* **W**: Weibull distribution ($\alpha, \lambda$)
* **L**: Lindley distribution ($q$)

---

## Repository Structure

```text
OL-GBP-Distribution-India/
├── data/
│   └── indian_states_indicators_2024.csv   # 33 observations each for Per Capita Income, GDP, and Birth Rate
├── R/
│   ├── 01_models_definition.R              # PDF and CDF functions for OL-GBP and baseline models
│   ├── 02_real_data_analysis.R             # Model fitting and goodness-of-fit comparison (LogLik, D, AIC, BIC, HQIC, p-value)
│   ├── 03_simulation_study.R               # Monte Carlo simulation for MLEs (MLE, MSE, RMSE, Absolute Bias, Variance)
│   └── 04_plot_empirical_fits.R            # Histogram density overlays and Total Time on Test (TTT) plots
├── python/
│   └── moments_and_hazards.py              # Numerical moments, survival, and hazard rate functions
├── .gitignore
├── LICENSE
└── README.md
