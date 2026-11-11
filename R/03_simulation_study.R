# ==============================================================================
# File: R/03_simulation_study.R
# Purpose: Monte Carlo simulation to evaluate MLE consistency (Tables 1, 2, 3)
# ==============================================================================

if (!require("AdequacyModel")) install.packages("AdequacyModel")
if (!require("GoFKernel")) install.packages("GoFKernel")
library(AdequacyModel)
library(GoFKernel)

source("R/01_models_definition.R")

run_mc_simulation <- function(true_par, sample_sizes = c(10, 30, 50, 70), n_rep = 1000, seed = 42) {
  set.seed(seed) # Global seed outside the replication loop
  
  a_true <- true_par[1]
  alpha_true <- true_par[2]
  beta_true <- true_par[3]
  
  # Numerical quantile function using GoFKernel inverse
  target_cdf <- function(x) cdf_olgbp(true_par, x)
  Q_numeric <- inverse(target_cdf, lower = 1e-5, upper = 500)
  
  all_results <- list()
  
  for (n in sample_sizes) {
    cat(sprintf("Running Monte Carlo: n = %d, reps = %d, par = (%.2f, %.2f, %.2f)...\n", 
                n, n_rep, a_true, alpha_true, beta_true))
    
    estimates <- matrix(NA, nrow = n_rep, ncol = 3)
    colnames(estimates) <- c("a", "alpha", "beta")
    
    for (i in 1:n_rep) {
      u <- runif(n)
      x_sim <- sapply(u, Q_numeric)
      
      fit_res <- tryCatch({
        goodness.fit(
          pdf = pdf_olgbp, cdf = cdf_olgbp,
          starts = true_par, data = x_sim,
          method = "PSO", domain = c(0, Inf),
          lim_inf = c(1e-4, 1e-4, 1e-4), lim_sup = c(10, 10, 10),
          N = 30
        )$mle
      }, error = function(e) rep(NA, 3))
      
      estimates[i, ] <- fit_res
    }
    
    valid_estimates <- estimates[complete.cases(estimates), , drop = FALSE]
    
    mean_mle <- colMeans(valid_estimates)
    bias <- mean_mle - true_par
    abs_bias <- abs(bias)
    mse <- colMeans(sweep(valid_estimates, 2, true_par)^2)
    rmse <- sqrt(mse)
    variance <- apply(valid_estimates, 2, var)
    
    res_df <- data.frame(
      SampleSize = n,
      Parameter = c("a", "alpha", "beta"),
      True = true_par,
      MLE = mean_mle,
      MSE = mse,
      RMSE = rmse,
      Absolute_Bias = abs_bias,
      Variance = variance
    )
    all_results[[as.character(n)]] <- res_df
  }
  
  final_table <- do.call(rbind, all_results)
  print(round(final_table[, 3:8], 6))
  return(final_table)
}

# Example run replicating Table 1 (Set 1: a=0.5, alpha=0.4, beta=0.4)
cat("\nReplicating Simulation: Parameter Set 1 (Table 1)\n")
sim_table1 <- run_mc_simulation(
  true_par = c(a = 0.5, alpha = 0.4, beta = 0.4),
  sample_sizes = c(10, 30, 50, 70),
  n_rep = 100 # Set to 1000 for full replication
)
