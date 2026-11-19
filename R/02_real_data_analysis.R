# ==============================================================================
# File: R/02_real_data_analysis.R
# Purpose: Fit OL-GBP, WBP, BP, Weibull, and Lindley to the 3 real datasets
#          and compute Log-Likelihood, KS D, AIC, BIC, HQIC, and p-value
# ==============================================================================

if (!require("AdequacyModel")) install.packages("AdequacyModel")
library(AdequacyModel)

source("R/01_models_definition.R")

df <- read.csv("data/indian_states_indicators_2024.csv")

run_comparison <- function(x, label, ol_start) {
  cat("\n============================================================\n")
  cat("Dataset:", label, "(n =", length(x), ")\n")
  cat("============================================================\n")
  print(summary(x))

  n <- length(x)

  fit_ol  <- goodness.fit(pdf = pdf_olgbp, cdf = cdf_olgbp, starts = ol_start,
                          data = x, method = "PSO", domain = c(0, Inf),
                          lim_inf = c(1e-8, 1e-4, 1e-4), lim_sup = c(20, 20, 20), N = 50)
  fit_wbp <- goodness.fit(pdf = pdf_wbp, cdf = cdf_wbp, starts = c(0.1, 0.5, 0.5, 0.5),
                          data = x, method = "PSO", domain = c(0, Inf),
                          lim_inf = c(1e-8, 1e-4, 1e-4, 1e-4), lim_sup = c(20, 20, 20, 20), N = 50)
  fit_bp  <- goodness.fit(pdf = pdf_bp, cdf = cdf_bp, starts = c(0.5, 0.5),
                          data = x, method = "PSO", domain = c(0, Inf),
                          lim_inf = c(1e-4, 1e-4), lim_sup = c(20, 20), N = 50)
  fit_w   <- goodness.fit(pdf = pdf_weibull, cdf = cdf_weibull, starts = c(1, 0.01),
                          data = x, method = "PSO", domain = c(0, Inf),
                          lim_inf = c(1e-4, 1e-7), lim_sup = c(20, 10), N = 50)
  fit_l   <- goodness.fit(pdf = pdf_lindley, cdf = cdf_lindley, starts = c(0.01),
                          data = x, method = "PSO", domain = c(0, Inf),
                          lim_inf = c(1e-7), lim_sup = c(10), N = 50)

  metrics_row <- function(model_name, fit, k) {
    ll <- -fit$Value
    d_val <- fit$KS$statistic
    p_val <- fit$KS$p.value
    aic <- -2 * ll + 2 * k
    bic <- -2 * ll + k * log(n)
    hqic <- -2 * ll + 2 * k * log(log(n))

    data.frame(
      Distribution = model_name,
      Log_Likelihood = round(ll, 4),
      D = round(d_val, 5),
      AIC = round(aic, 4),
      BIC = round(bic, 4),
      HQIC = round(hqic, 4),
      p_value = format.pval(p_val, digits = 4),
      stringsAsFactors = FALSE
    )
  }

  gof_table <- rbind(
    metrics_row("OL-GBP", fit_ol, 3),
    metrics_row("WBP", fit_wbp, 4),
    metrics_row("BP", fit_bp, 2),
    metrics_row("W", fit_w, 2),
    metrics_row("L", fit_l, 1)
  )

  cat("\nParameter Estimates:\n")
  cat("OL-GBP:", paste(names(ol_start), round(fit_ol$mle, 7), sep = "=", collapse = ", "), "\n")
  cat("WBP:   ", paste(round(fit_wbp$mle, 7), collapse = ", "), "\n")
  cat("BP:    ", paste(round(fit_bp$mle, 7), collapse = ", "), "\n")
  cat("W:     ", paste(round(fit_w$mle, 7), collapse = ", "), "\n")
  cat("L:     ", paste(round(fit_l$mle, 7), collapse = ", "), "\n")

  cat("\nGoodness-of-Fit Summary:\n")
  print(gof_table, row.names = FALSE)
}

run_comparison(df$per_capita_income, "Per Capita Income", c(a = 0.001, alpha = 0.6, beta = 0.5))
run_comparison(df$gdp, "GDP", c(a = 0.007, alpha = 0.7, beta = 0.3))
run_comparison(df$birth_rate, "Birth Rate", c(a = 0.9, alpha = 0.3, beta = 0.9))
