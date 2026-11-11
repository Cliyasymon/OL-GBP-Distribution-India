# ==============================================================================
# File: R/01_models_definition.R
# Purpose: Probability Density Functions (PDF) and Cumulative Distribution 
#          Functions (CDF) for OL-GBP and baseline distributions.
# Compatible with AdequacyModel::goodness.fit()
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Odd Lindley Generalized Beta Prime (OL-GBP) Distribution
# Parameters: par = c(a, alpha, beta)
# Support: x > 0, a > 0, alpha > 0, beta > 0
# ------------------------------------------------------------------------------
pdf_olgbp <- function(par, x) {
  a     <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  if (a <= 0 || alpha <= 0 || beta <= 0) {
    return(rep(1e-15, length(x)))
  }
  
  # Support check
  dens <- rep(1e-15, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(dens)
  
  x_pos <- x[pos_idx]
  z <- x_pos / (1 + x_pos)
  
  # Baseline Beta Prime CDF: G(x) = I_z(alpha, beta)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  diff_G <- pmax(1 - G, 1e-12)
  odds <- G / diff_G
  
  # Baseline Beta Prime PDF: g(x) = x^(alpha-1) * (1+x)^(-alpha-beta) / B(alpha, beta)
  log_g <- (alpha - 1) * log(x_pos) - (alpha + beta) * log(1 + x_pos) - lbeta(alpha, beta)
  g <- exp(log_g)
  
  # Compounded OL-GBP PDF
  dens[pos_idx] <- (a^2 / (1 + a)) * (g / (diff_G^3)) * exp(-a * odds)
  return(pmax(dens, 1e-15))
}

cdf_olgbp <- function(par, x) {
  a     <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  if (a <= 0 || alpha <= 0 || beta <= 0) {
    return(rep(0, length(x)))
  }
  
  prob <- rep(0, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(prob)
  
  x_pos <- x[pos_idx]
  z <- x_pos / (1 + x_pos)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  diff_G <- pmax(1 - G, 1e-12)
  odds <- G / diff_G
  
  term <- (a + 1 - G) / ((1 + a) * diff_G)
  prob[pos_idx] <- 1 - term * exp(-a * odds)
  return(pmin(pmax(prob, 0), 1))
}


# ------------------------------------------------------------------------------
# 2. Weibull-Beta Prime (WBP) Distribution
# Parameters: par = c(a, b, alpha, beta)
# Support: x > 0, a > 0, b > 0, alpha > 0, beta > 0
# ------------------------------------------------------------------------------
pdf_wbp <- function(par, x) {
  a     <- par[1]
  b     <- par[2]
  alpha <- par[3]
  beta  <- par[4]
  
  if (any(par <= 0)) {
    return(rep(1e-15, length(x)))
  }
  
  dens <- rep(1e-15, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(dens)
  
  x_pos <- x[pos_idx]
  z <- x_pos / (1 + x_pos)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  diff_G <- pmax(1 - G, 1e-12)
  odds <- G / diff_G
  
  log_g <- (alpha - 1) * log(x_pos) - (alpha + beta) * log(1 + x_pos) - lbeta(alpha, beta)
  g <- exp(log_g)
  
  dens[pos_idx] <- (a * b * g * (G^(b - 1))) / (diff_G^(b + 1)) * exp(-a * (odds^b))
  return(pmax(dens, 1e-15))
}

cdf_wbp <- function(par, x) {
  a     <- par[1]
  b     <- par[2]
  alpha <- par[3]
  beta  <- par[4]
  
  if (any(par <= 0)) {
    return(rep(0, length(x)))
  }
  
  prob <- rep(0, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(prob)
  
  x_pos <- x[pos_idx]
  z <- x_pos / (1 + x_pos)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  diff_G <- pmax(1 - G, 1e-12)
  odds <- G / diff_G
  
  prob[pos_idx] <- 1 - exp(-a * (odds^b))
  return(pmin(pmax(prob, 0), 1))
}


# ------------------------------------------------------------------------------
# 3. Standard Beta Prime (BP) Distribution
# Parameters: par = c(alpha, beta)
# Support: x > 0, alpha > 0, beta > 0
# ------------------------------------------------------------------------------
pdf_bp <- function(par, x) {
  alpha <- par[1]
  beta  <- par[2]
  
  if (alpha <= 0 || beta <= 0) {
    return(rep(1e-15, length(x)))
  }
  
  dens <- rep(1e-15, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(dens)
  
  x_pos <- x[pos_idx]
  log_g <- (alpha - 1) * log(x_pos) - (alpha + beta) * log(1 + x_pos) - lbeta(alpha, beta)
  dens[pos_idx] <- exp(log_g)
  return(pmax(dens, 1e-15))
}

cdf_bp <- function(par, x) {
  alpha <- par[1]
  beta  <- par[2]
  
  if (alpha <= 0 || beta <= 0) {
    return(rep(0, length(x)))
  }
  
  prob <- rep(0, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(prob)
  
  x_pos <- x[pos_idx]
  z <- x_pos / (1 + x_pos)
  prob[pos_idx] <- pbeta(z, shape1 = alpha, shape2 = beta)
  return(pmin(pmax(prob, 0), 1))
}


# ------------------------------------------------------------------------------
# 4. Weibull (W) Distribution
# Parameters: par = c(alpha, lambda)
# Support: x > 0, alpha > 0, lambda > 0
# Formulation: f(x) = alpha * lambda * x^(alpha - 1) * exp(-lambda * x^alpha)
# ------------------------------------------------------------------------------
pdf_weibull <- function(par, x) {
  alpha  <- par[1]
  lambda <- par[2]
  
  if (alpha <= 0 || lambda <= 0) {
    return(rep(1e-15, length(x)))
  }
  
  dens <- rep(1e-15, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(dens)
  
  x_pos <- x[pos_idx]
  dens[pos_idx] <- alpha * lambda * (x_pos^(alpha - 1)) * exp(-lambda * (x_pos^alpha))
  return(pmax(dens, 1e-15))
}

cdf_weibull <- function(par, x) {
  alpha  <- par[1]
  lambda <- par[2]
  
  if (alpha <= 0 || lambda <= 0) {
    return(rep(0, length(x)))
  }
  
  prob <- rep(0, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(prob)
  
  x_pos <- x[pos_idx]
  prob[pos_idx] <- 1 - exp(-lambda * (x_pos^alpha))
  return(pmin(pmax(prob, 0), 1))
}


# ------------------------------------------------------------------------------
# 5. Lindley (L) Distribution
# Parameters: par = c(q)
# Support: x > 0, q > 0
# Formulation: f(x) = (q^2 / (1 + q)) * (1 + x) * exp(-q * x)
# ------------------------------------------------------------------------------
pdf_lindley <- function(par, x) {
  q <- par[1]
  
  if (q <= 0) {
    return(rep(1e-15, length(x)))
  }
  
  dens <- rep(1e-15, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(dens)
  
  x_pos <- x[pos_idx]
  dens[pos_idx] <- (q^2 / (1 + q)) * (1 + x_pos) * exp(-q * x_pos)
  return(pmax(dens, 1e-15))
}

cdf_lindley <- function(par, x) {
  q <- par[1]
  
  if (q <= 0) {
    return(rep(0, length(x)))
  }
  
  prob <- rep(0, length(x))
  pos_idx <- which(x > 0)
  if (length(pos_idx) == 0) return(prob)
  
  x_pos <- x[pos_idx]
  prob[pos_idx] <- 1 - (1 + (q * x_pos) / (1 + q)) * exp(-q * x_pos)
  return(pmin(pmax(prob, 0), 1))
}
