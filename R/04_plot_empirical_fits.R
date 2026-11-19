# ==============================================================================
# File: R/04_plot_empirical_fits.R
# Purpose: Generate Total Time on Test (TTT) plots and fitted OL-GBP density
#          overlays for the 3 datasets
# ==============================================================================

source("R/01_models_definition.R")

df <- read.csv("data/indian_states_indicators_2024.csv")

plot_ttt <- function(x, main_title, file_name) {
  n <- length(x)
  xs <- sort(x)
  r <- seq_len(n)
  t_val <- sapply(r, function(k) (sum(xs[1:k]) + (n - k) * xs[k]) / sum(xs))

  png(file_name, width = 1500, height = 1300, res = 250)
  plot(r / n, t_val, type = "l", lwd = 2, col = "black",
       xlab = "i / n", ylab = "T(i/n)", main = main_title)
  abline(0, 1, lty = 2, col = "gray50")
  dev.off()
}

plot_density_fit <- function(x, mle_par, main_title, x_range, file_name) {
  png(file_name, width = 1600, height = 1300, res = 250)
  hist(x, breaks = 10, probability = TRUE, col = "#74c69d", border = "white",
       xlab = "Data", main = main_title, xlim = x_range)
  grid(col = "gray90")
  hist(x, breaks = 10, probability = TRUE, col = "#74c69d", border = "white", add = TRUE)

  x_grid <- seq(max(min(x) * 0.5, 1e-4), max(x) * 1.05, length.out = 400)
  lines(x_grid, pdf_olgbp(mle_par, x_grid), col = "#d62828", lwd = 2.5)
  legend("topright", legend = c("Data", "OL-GBP"),
         col = c("#74c69d", "#d62828"), lwd = c(8, 2.5), bty = "n")
  dev.off()
}

# Fitted parameters from Sections 8.1.1 - 8.1.3 of the manuscript
plot_density_fit(df$per_capita_income, c(0.001031915, 0.622944486, 0.581670882),
                 "Per Capita Income", c(0, max(df$per_capita_income) * 1.05), "pci_pdf.png")
plot_ttt(df$per_capita_income, "TTT Plot: Per Capita Income", "pci_ttt.png")

plot_density_fit(df$gdp, c(0.007859651, 0.751781857, 0.341639575),
                 "GDP", c(0, max(df$gdp) * 1.05), "gdp_pdf.png")
plot_ttt(df$gdp, "TTT Plot: GDP", "gdp_ttt.png")

plot_density_fit(df$birth_rate, c(0.9037800, 0.3234193, 0.9760054),
                 "Birth Rate", c(0, 0.5), "birth_rate_pdf.png")
plot_ttt(df$birth_rate, "TTT Plot: Birth Rate", "birth_rate_ttt.png")
