# ==============================================================================
# File: R/01_models_definition.R
# Purpose: PDF and CDF functions for OL-GBP and comparison distributions
#          (WBP, BP, Weibull, Lindley)
# ==============================================================================

# 1. Odd Lindley-G Beta Prime (OL-GBP): par = c(a, alpha, beta)
pdf_olgbp <- function(par, x) {
  a <- par[1]; alpha <- par[2]; beta <- par[3]
  if (a <= 0 || alpha <= 0 || beta <= 0) return(rep(1e-15, length(x)))

  out <- rep(1e-15, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  z <- xp / (1 + xp)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  one_minus_G <- pmax(1 - G, 1e-12)
  log_g <- (alpha - 1) * log(xp) - (alpha + beta) * log(1 + xp) - lbeta(alpha, beta)
  g <- exp(log_g)

  out[idx] <- (a^2 / (1 + a)) * (g / (one_minus_G^3)) * exp(-a * (G / one_minus_G))
  return(pmax(out, 1e-15))
}

cdf_olgbp <- function(par, x) {
  a <- par[1]; alpha <- par[2]; beta <- par[3]
  if (a <= 0 || alpha <= 0 || beta <= 0) return(rep(0, length(x)))

  out <- rep(0, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  z <- xp / (1 + xp)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  one_minus_G <- pmax(1 - G, 1e-12)

  surv <- ((a + 1 - G) / ((1 + a) * one_minus_G)) * exp(-a * (G / one_minus_G))
  out[idx] <- 1 - surv
  return(pmin(pmax(out, 0), 1))
}

# 2. Weibull-Beta Prime (WBP): par = c(a, b, alpha, beta)
pdf_wbp <- function(par, x) {
  a <- par[1]; b <- par[2]; alpha <- par[3]; beta <- par[4]
  if (any(par <= 0)) return(rep(1e-15, length(x)))

  out <- rep(1e-15, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  z <- xp / (1 + xp)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  one_minus_G <- pmax(1 - G, 1e-12)
  odds <- G / one_minus_G
  log_g <- (alpha - 1) * log(xp) - (alpha + beta) * log(1 + xp) - lbeta(alpha, beta)
  g <- exp(log_g)

  out[idx] <- (a * b * g * (G^(b - 1))) / (one_minus_G^(b + 1)) * exp(-a * (odds^b))
  return(pmax(out, 1e-15))
}

cdf_wbp <- function(par, x) {
  a <- par[1]; b <- par[2]; alpha <- par[3]; beta <- par[4]
  if (any(par <= 0)) return(rep(0, length(x)))

  out <- rep(0, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  z <- xp / (1 + xp)
  G <- pbeta(z, shape1 = alpha, shape2 = beta)
  one_minus_G <- pmax(1 - G, 1e-12)
  odds <- G / one_minus_G

  out[idx] <- 1 - exp(-a * (odds^b))
  return(pmin(pmax(out, 0), 1))
}

# 3. Beta Prime (BP): par = c(alpha, beta)
pdf_bp <- function(par, x) {
  alpha <- par[1]; beta <- par[2]
  if (alpha <= 0 || beta <= 0) return(rep(1e-15, length(x)))

  out <- rep(1e-15, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  log_g <- (alpha - 1) * log(xp) - (alpha + beta) * log(1 + xp) - lbeta(alpha, beta)
  out[idx] <- exp(log_g)
  return(pmax(out, 1e-15))
}

cdf_bp <- function(par, x) {
  alpha <- par[1]; beta <- par[2]
  if (alpha <= 0 || beta <= 0) return(rep(0, length(x)))

  out <- rep(0, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  z <- xp / (1 + xp)
  out[idx] <- pbeta(z, shape1 = alpha, shape2 = beta)
  return(pmin(pmax(out, 0), 1))
}

# 4. Weibull (W): par = c(alpha, lambda)
pdf_weibull <- function(par, x) {
  alpha <- par[1]; lambda <- par[2]
  if (alpha <= 0 || lambda <= 0) return(rep(1e-15, length(x)))

  out <- rep(1e-15, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  out[idx] <- alpha * lambda * (xp^(alpha - 1)) * exp(-lambda * (xp^alpha))
  return(pmax(out, 1e-15))
}

cdf_weibull <- function(par, x) {
  alpha <- par[1]; lambda <- par[2]
  if (alpha <= 0 || lambda <= 0) return(rep(0, length(x)))

  out <- rep(0, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  out[idx] <- 1 - exp(-lambda * (xp^alpha))
  return(pmin(pmax(out, 0), 1))
}

# 5. Lindley (L): par = c(q)
pdf_lindley <- function(par, x) {
  q <- par[1]
  if (q <= 0) return(rep(1e-15, length(x)))

  out <- rep(1e-15, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  out[idx] <- (q^2 / (1 + q)) * (1 + xp) * exp(-q * xp)
  return(pmax(out, 1e-15))
}

cdf_lindley <- function(par, x) {
  q <- par[1]
  if (q <= 0) return(rep(0, length(x)))

  out <- rep(0, length(x))
  idx <- which(x > 0)
  if (length(idx) == 0) return(out)

  xp <- x[idx]
  out[idx] <- 1 - (1 + (q * xp) / (1 + q)) * exp(-q * xp)
  return(pmin(pmax(out, 0), 1))
}
