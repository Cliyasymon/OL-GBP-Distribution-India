# ==============================================================================
# File: R/03_simulation_study.R
# Purpose: Monte Carlo simulation to evaluate MLEs (MLE, MSE, RMSE,
#          Absolute Bias, and Variance) for OL-GBP parameters
# ==============================================================================

if (!require("AdequacyModel")) install.packages("AdequacyModel")
if (!require("GoFKernel")) install.packages("GoFKernel")
library(AdequacyModel)
library(GoFKernel)

source("R/01_models_definition.R")

simulate_olgbp <- function(true_par, n_vals = c(10, 30, 50, 70), n_rep = 1000, seed = 1000) {
  set.seed(seed) # Set once outside replication loops for distinct Monte Carlo draws

  q_fn <- inverse(function(x) cdf_olgbp(true_par, x), lower = 1e-5, upper = 500)
  out_list <- list()

  for (n in n_vals) {
    est_mat <- matrix(NA_real_, nrow = n_rep, ncol = 3)
    colnames(est_mat) <- c("a", "alpha", "beta")

    for (r in seq_len(n_rep)) {
      u <- runif(n)
      x_sim <- sapply(u, q_fn)

      fit <- tryCatch(
        goodness.fit(
          pdf = pdf_olgbp, cdf = cdf_olgbp,
          starts = true_par, data = x_sim,
          method = "PSO", domain = c(0, Inf),
          lim_inf = c(1e-4, 1e-4, 1e-4), lim_sup = c(10, 10, 10),
          N = 30
        )$mle,
        error = function(e) rep(NA_real_, 3)
      )
      est_mat[r, ] <- fit
    }

    valid <- est_mat[complete.cases(est_mat), , drop = FALSE]
    mean_mle <- colMeans(valid)
    bias <- mean_mle - true_par
    abs_bias <- abs(bias)
    mse <- colMeans(sweep(valid, 2, true_par)^2)
    rmse <- sqrt(mse)
    variance <- apply(valid, 2, var)

    out_list[[as.character(n)]] <- data.frame(
      Sample_Size = n,
      Parameter = c("a", "alpha", "beta"),
      True = as.numeric(true_par),
      MLE = round(mean_mle, 6),
      MSE = round(mse, 6),
      RMSE = round(rmse, 6),
      Absolute_Bias = round(abs_bias, 6),
      Variance = round(variance, 6),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_list)
}

# Example parameter sets from the manuscript simulation tables
param_sets <- list(
  set_1 = c(a = 0.5, alpha = 0.40, beta = 0.40),
  set_2 = c(a = 0.9, alpha = 0.40, beta = 0.55),
  set_3 = c(a = 0.6, alpha = 0.65, beta = 0.75)
)

# Use n_rep = 1000 for full paper study (reduce to e.g. 50 for a fast local check)
res_set1 <- simulate_olgbp(param_sets$set_1, n_vals = c(10, 30, 50, 70), n_rep = 1000)
print(res_set1, row.names = FALSE)
