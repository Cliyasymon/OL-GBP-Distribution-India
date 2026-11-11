# ==============================================================================
# File: R/04_plot_empirical_fits.R
# Purpose: Generate fitted density histograms and Total Time on Test (TTT) plots
# Replicates: Figures 9, 10, and 11
# ==============================================================================

source("R/01_models_definition.R")

df <- read.csv("data/indian_states_data_2024.csv")

generate_diagnostics <- function(data_vector, title_label, fitted_ol_par, x_limit, out_prefix) {
  n <- length(data_vector)
  y_sorted <- sort(data_vector)
  i_vec <- 1:n
  
  # 1. Total Time on Test (TTT) Transform
  phi <- sapply(i_vec, function(r) {
    (sum(y_sorted[1:r]) + (n - r) * y_sorted[r]) / sum(y_sorted)
  })
  
  png(paste0(out_prefix, "_ttt.png"), width = 1600, height = 1400, res = 300)
  plot(i_vec / n, phi, type = "l", lwd = 2, col = "black",
       xlab = "i / n", ylab = "T(i/n)", main = paste("TTT Plot -", title_label))
  abline(0, 1, lty = 2, col = "gray50")
  dev.off()
  
  # 2. Fitted Histogram Overlay
  png(paste0(out_prefix, "_density.png"), width = 1800, height = 1400, res = 300)
  hist(data_vector, breaks = 12, probability = TRUE, col = "#76c893", border = "white",
       xlab = "Data", main = paste("Fitted Density -", title_label), xlim = x_limit)
  grid(col = "gray90")
  hist(data_vector, breaks = 12, probability = TRUE, col = "#76c893", border = "white", add = TRUE)
  
  x_seq <- seq(min(data_vector) * 0.8, max(data_vector) * 1.05, length.out = 500)
  lines(x_seq, pdf_olgbp(fitted_ol_par, x_seq), col = "#e63946", lwd = 2.5)
  legend("topright", legend = c("Histogram", "OL-GBP Fit"),
         col = c("#76c893", "#e63946"), lwd = c(8, 2.5), bty = "n")
  dev.off()
  
  cat(sprintf("Diagnostics saved: %s_ttt.png and %s_density.png\n", out_prefix, out_prefix))
}

# Generate plots for Birth Rate
generate_diagnostics(
  data_vector = df$Birth_Rate,
  title_label = "Birth Rate (2024)",
  fitted_ol_par = c(0.9037800, 0.3234193, 0.9760054),
  x_limit = c(0, 0.5),
  out_prefix = "birth_rate"
)
