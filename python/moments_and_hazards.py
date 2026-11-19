"""
Numerical functions for OL-GBP PDF, CDF, Survival, Hazard Rate,
and first four moments (Mean, Variance, Skewness, Kurtosis).
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import betainc, gammaln


def _log_beta_prime_pdf(x, alpha, beta):
    log_b = gammaln(alpha) + gammaln(beta) - gammaln(alpha + beta)
    return (alpha - 1.0) * np.log(x) - (alpha + beta) * np.log(1.0 + x) - log_b


def olgbp_pdf(x, a, alpha, beta):
    if x <= 0.0:
        return 0.0
    z = x / (1.0 + x)
    g_cdf = float(betainc(alpha, beta, z))
    bar_g = max(1.0 - g_cdf, 1e-15)
    g_pdf = np.exp(_log_beta_prime_pdf(x, alpha, beta))
    return (a**2 / (1.0 + a)) * (g_pdf / (bar_g**3)) * np.exp(-a * (g_cdf / bar_g))


def olgbp_cdf(x, a, alpha, beta):
    if x <= 0.0:
        return 0.0
    z = x / (1.0 + x)
    g_cdf = float(betainc(alpha, beta, z))
    bar_g = max(1.0 - g_cdf, 1e-15)
    survival = ((a + 1.0 - g_cdf) / ((1.0 + a) * bar_g)) * np.exp(-a * (g_cdf / bar_g))
    return float(np.clip(1.0 - survival, 0.0, 1.0))


def olgbp_survival(x, a, alpha, beta):
    return 1.0 - olgbp_cdf(x, a, alpha, beta)


def olgbp_hazard(x, a, alpha, beta):
    s_val = olgbp_survival(x, a, alpha, beta)
    if s_val <= 1e-15:
        return 0.0
    return olgbp_pdf(x, a, alpha, beta) / s_val


def olgbp_moments(a, alpha, beta):
    mean, _ = quad(lambda t: t * olgbp_pdf(t, a, alpha, beta), 0.0, np.inf)
    m2, _ = quad(lambda t: (t**2) * olgbp_pdf(t, a, alpha, beta), 0.0, np.inf)
    m3, _ = quad(lambda t: (t**3) * olgbp_pdf(t, a, alpha, beta), 0.0, np.inf)
    m4, _ = quad(lambda t: (t**4) * olgbp_pdf(t, a, alpha, beta), 0.0, np.inf)

    var = m2 - mean**2
    sd = np.sqrt(var)
    skew = (m3 - 3.0 * mean * m2 + 2.0 * (mean**3)) / (sd**3)
    kurt = (m4 - 4.0 * mean * m3 + 6.0 * (mean**2) * m2 - 3.0 * (mean**4)) / (var**2)
    return {"Mean": mean, "Variance": var, "Skewness": skew, "Kurtosis": kurt}


if __name__ == "__main__":
    summary = olgbp_moments(a=0.5, alpha=2.0, beta=3.0)
    for k, v in summary.items():
        print(f"{k}: {v:.6f}")
