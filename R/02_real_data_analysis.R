# ==============================================================================
# File: R/02_real_data_analysis.R
# Purpose: Model fitting and goodness-of-fit evaluation for Indian states data
# Replicates: Tables 4–12 (Per Capita Income, GDP, Birth Rate)
# ==============================================================================

if (!require("AdequacyModel")) install.packages("AdequacyModel")
library(AdequacyModel)

# Source distribution definitions
source("R/01_models_definition.R")

# Load socio-economic data
data_path <- "data/indian_states_data_2024.csv"
if (!file.exists(data_path)) {
  stop("Dataset not found. Please ensure data/indian_states_data_2024.csv exists.")
}
df <- read.csv(data_path)

# Generalized fitting function for AdequacyModel
fit_benchmark_models <- function(data_vector, dataset_name, init_olgbp) {
  cat("\n====================================================================\n")
  cat(paste("MODEL FITTING & EVALUATION:", toupper(dataset_name), "\n"))
  cat("====================================================================\n")
  
  n <- length(data_vector)
  
  # 1. Odd Lindley-G Beta Prime (OL-GBP)
  fit_ol <- goodness.fit(
    pdf = pdf_olgbp, cdf = cdf_olgbp, starts = init_olgbp,
    data = data_vector, method = "PSO", domain = c(0, Inf),
    lim_inf = c(1e-7, 1e-4, 1e-4), lim_sup = c(10, 10, 10), N = 50
  )
  
  # 2. Weibull-Beta Prime (WBP)
  fit_wbp <- goodness.fit(
    pdf = pdf_wbp, cdf = cdf_wbp, starts = c(0.1, 0.5, 0.5, 0.5),
    data = data_vector, method = "PSO", domain = c(0, Inf),
    lim_inf = c(1e-7, 1e-4, 1e-4, 1e-4), lim_sup = c(10, 10, 10, 10), N = 50
  )
  
  # 3. Beta Prime (BP)
  fit_bp <- goodness.fit(
    pdf = pdf_bp, cdf = cdf_bp, starts = c(0.5, 0.5),
    data = data_vector, method = "PSO", domain = c(0, Inf),
    lim_inf = c(1e-4, 1e-4), lim_sup = c(10, 10), N = 50
  )
  
  # 4. Weibull (W)
  fit_w <- goodness.fit(
    pdf = pdf_weibull, cdf = cdf_weibull, starts = c(1, 0.01),
    data = data_vector, method = "PSO", domain = c(0, Inf),
    lim_inf = c(1e-4, 1e-6), lim_sup = c(10, 5), N = 50
  )
  
  # 5. Lindley (L)
  fit_l <- goodness.fit(
    pdf = pdf_lindley, cdf = cdf_lindley, starts = c(0.01),
    data = data_vector, method = "PSO", domain = c(0, Inf),
    lim_inf = c(1e-6), lim_sup = c(5), N = 50
  )
  
  # Helper to summarize goodness-of-fit metrics
  summarize_model <- function(name, fit_obj, k) {
    loglik <- -fit_obj$Value
    ks_d <- fit_obj$KS$statistic
    p_val <- fit_obj$KS$p.value
    aic <- 2 * k - 2 * loglik
    bic <- k * log(n) - 2 * loglik
    hqic <- 2 * k * log(log(n)) - 2 * loglik
    
    data.frame(
      Model = name,
      LogLikelihood = round(loglik, 4),
      D_KS = round(ks_d, 5),
      AIC = round(aic, 4),
      BIC = round(bic, 4),
      HQIC = round(hqic, 4),
      p_value = format.pval(p_val, digits = 4)
    )
  }
  
  results_table <- rbind(
    summarize_model("OL-GBP", fit_ol, 3),
    summarize_model("WBP", fit_wbp, 4),
    summarize_model("BP", fit_bp, 2),
    summarize_model("Weibull", fit_w, 2),
    summarize_model("Lindley", fit_l, 1)
  )
  
  print(results_table)
  return(list(OL_GBP = fit_ol, WBP = fit_wbp, BP = fit_bp, Weibull = fit_w, Lindley = fit_l))
}

# --- Execution Across Empirical Datasets ---
# 1. Per Capita Income (Table 5 & 6)
cat("\nFitting Per Capita Income Data...\n")
res_income <- fit_benchmark_models(df$Per_Capita_Income, "Per Capita Income", c(0.001, 0.6, 0.5))

# 2. State GDP (Table 8 & 9)
cat("\nFitting State GDP Data...\n")
res_gdp <- fit_benchmark_models(df$GDP, "State GDP", c(0.007, 0.7, 0.3))

# 3. Birth Rate (Table 11 & 12)
cat("\nFitting Crude Birth Rate Data...\n")
res_birth <- fit_benchmark_models(df$Birth_Rate, "Birth Rate", c(0.9, 0.3, 0.9))
