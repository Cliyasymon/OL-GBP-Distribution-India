"""
Analytical properties and numerical integration for Odd Lindley-G Beta Prime (OL-GBP)
Computes theoretical PDF, CDF, survival rate, hazard rate, and numerical moments.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import beta as B_func, betainc
from scipy.integrate import quad

def olgbp_pdf(x, a, alpha, beta):
    """Probability density function of OL-GBP distribution."""
    if x <= 0:
        return 0.0
    z = x / (1.0 + x)
    G = betainc(alpha, beta, z)
    diff = max(1.0 - G, 1e-15)
    B_val = B_func(alpha, beta)
    g = (x**(alpha - 1.0) * (1.0 + x)**(-alpha - beta)) / B_val
    return (a**2 / (1.0 + a)) * (g / (diff**3)) * np.exp(-a * (G / diff))

def olgbp_cdf(x, a, alpha, beta):
    """Cumulative distribution function of OL-GBP distribution."""
    if x <= 0:
        return 0.0
    z = x / (1.0 + x)
    G = betainc(alpha, beta, z)
    diff = max(1.0 - G, 1e-15)
    term = (a + 1.0 - G) / ((1.0 + a) * diff)
    return np.clip(1.0 - term * np.exp(-a * (G / diff)), 0.0, 1.0)

def olgbp_survival(x, a, alpha, beta):
    """Survival function S(x) = 1 - F(x)."""
    return 1.0 - olgbp_cdf(x, a, alpha, beta)

def olgbp_hazard(x, a, alpha, beta):
    """Hazard rate function h(x) = f(x) / S(x)."""
    S_val = olgbp_survival(x, a, alpha, beta)
    if S_val <= 1e-15:
        return 0.0
    return olgbp_pdf(x, a, alpha, beta) / S_val

def compute_moments(a, alpha, beta):
    """Computes Mean, Variance, Skewness, and Kurtosis numerically."""
    mean, _ = quad(lambda x: x * olgbp_pdf(x, a, alpha, beta), 0, np.inf)
    m2, _ = quad(lambda x: (x**2) * olgbp_pdf(x, a, alpha, beta), 0, np.inf)
    m3, _ = quad(lambda x: (x**3) * olgbp_pdf(x, a, alpha, beta), 0, np.inf)
    m4, _ = quad(lambda x: (x**4) * olgbp_pdf(x, a, alpha, beta), 0, np.inf)
    
    variance = m2 - mean**2
    sigma = np.sqrt(variance)
    skewness = (m3 - 3 * mean * m2 + 2 * mean**3) / (sigma**3)
    kurtosis = (m4 - 4 * mean * m3 + 6 * (mean**2) * m2 - 3 * mean**4) / (variance**2)
    
    return {
        "Mean": mean,
        "Variance": variance,
        "Skewness": skewness,
        "Kurtosis": kurtosis
    }

if __name__ == "__main__":
    test_params = {"a": 0.5, "alpha": 2.0, "beta": 3.0}
    stats = compute_moments(**test_params)
    print("Numerical Properties for OL-GBP (a=0.5, alpha=2.0, beta=3.0):")
    for k, v in stats.items():
        print(f"  {k}: {v:.6f}")
