"""
xray_uplim.statistics
---------------------
Count-rate estimation methods shared across all observatory modules.

Three methods are implemented:

net_count_rate()
    Background-subtracted point estimate.  Not an upper limit.
    Useful as a sanity check and for comparison with the proper ULs.

marginalized_upper_limit()
    Bayesian upper limit marginalising over the unknown background rate B.
    Treats n~Poisson(S·t + α·B·t) and m~Poisson(B·t) and integrates over
    B analytically (weighted mixture of Gamma distributions).
    Supersedes kraft_upper_limit(), which fixes B at its point estimate.
    As m → ∞ the two methods converge.  Used across all telescopes.

gehrels_upper_limit()
    Gehrels 1986, ApJ 303, 336.
    Closed-form Poisson approximation.  Slightly overestimates at low N.
    Printed as a cross-check alongside the marginalized UL.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm as sp_norm
from scipy.special import gammaincc, gammaln, gammainc


def net_count_rate(N_src, B_scaled, t_eff, N_bkg_raw, area_ratio):
    """
    Background-subtracted net count rate.

        CR_net = (N_src - B_scaled) / t_eff

    This is a POINT ESTIMATE, not an upper limit.  A negative value is
    perfectly valid and indicates a clean non-detection (the source
    aperture collected fewer counts than the expected background).

    Uncertainty
    -----------
    Correct propagation through the area scaling:

        Var(B_scaled) = (A_src / A_bkg)^2 * N_bkg_raw

        sigma_counts  = sqrt(N_src + N_bkg_raw * area_ratio^2)
        sigma_CR      = sigma_counts / t_eff

    Note: using sqrt(N_src + B_scaled) is only correct when area_ratio = 1,
    which is never the case for a source circle + background annulus.

    Parameters
    ----------
    N_src      : int    — counts in source aperture
    B_scaled   : float  — expected background in source aperture
                          (= N_bkg_raw * area_ratio)
    t_eff      : float  — effective exposure time in seconds
    N_bkg_raw  : int    — raw background counts (before area scaling)
    area_ratio : float  — A_src / A_bkg

    Returns
    -------
    CR    : float  — net count rate in cts/s
    sigma : float  — 1-sigma Poisson uncertainty in cts/s
    """
    CR    = (N_src - B_scaled) / t_eff
    sigma = np.sqrt(N_src + N_bkg_raw * area_ratio**2) / t_eff
    return CR, sigma


def kraft_upper_limit(N_obs, B_scaled, confidence):
    """
    Kraft, Burrows & Nousek 1991 Bayesian upper limit on source signal S.

    EXACT SOLUTION via the regularised incomplete Gamma function
    ------------------------------------------------------------
    Prior      : uniform on S >= 0
    Likelihood : Poisson(N; S + B)
    Posterior  : p(S | N, B)  proportional to  (S+B)^N * exp(-(S+B))

    With the substitution lambda = S + B the normalisation integral is:

        norm  =  integral_B^inf  lambda^N * exp(-lambda)  d(lambda)
              =  Gamma(N+1, B)                [upper incomplete gamma]

    The posterior CDF therefore has the closed form:

        P(S <= s_up | N, B)
            = 1 - gammaincc(N+1, s_up+B) / gammaincc(N+1, B)

    where gammaincc(a, x) = Gamma(a, x) / Gamma(a) is scipy's regularised
    upper incomplete gamma function, evaluated in log-space internally and
    numerically stable for any N.

    Parameters
    ----------
    N_obs      : int    — counts observed in source aperture
    B_scaled   : float  — expected background in source aperture
    confidence : float  — one-sided confidence level (e.g. 0.9973)

    Returns
    -------
    S_upper : float  — upper limit on net source counts
    """
    a = float(N_obs) + 1.0
    B = float(B_scaled)

    norm_reg = gammaincc(a, B)
    if norm_reg == 0.0:
        return 0.0   # B >> N: entire posterior mass is at S ~ 0

    target = (1.0 - confidence) * norm_reg

    def equation(s_up):
        return gammaincc(a, s_up + B) - target

    # Build a bracket: s=0 gives norm_reg > target (since confidence < 1).
    # Expand s_hi until gammaincc falls below target.
    s_hi = max(float(N_obs), 1.0) + 10.0 * np.sqrt(max(float(N_obs), 1.0)) + 50.0
    for _ in range(40):
        if gammaincc(a, s_hi + B) < target:
            break
        s_hi *= 2.0

    try:
        S_upper = brentq(equation, 0.0, s_hi, xtol=1e-5, maxiter=500)
    except ValueError:
        S_upper = s_hi   # fallback; should not occur with the bracket above

    return S_upper


def marginalized_upper_limit(n_obs, n_bkg_raw, area_ratio, t_eff, confidence):
    """
    Bayesian upper limit that marginalises over the unknown background rate.

    This is what CIAO aprates implements for Chandra, generalised to work
    for any telescope.  It supersedes kraft_upper_limit() by treating the
    background as an unknown Poisson parameter rather than a fixed number.

    Statistical model
    -----------------
    Source aperture  :  n ~ Poisson(S·t  +  α·B·t)
    Background aperture :  m ~ Poisson(B·t)

    where S = source rate (cts/s), B = background rate (cts/s),
    α = A_src / A_bkg (area ratio), t = t_eff (equal for both apertures
    in the same observation; for co-added obs use the summed exposure).

    Flat (improper) priors on both S ≥ 0 and B ≥ 0.

    The marginal posterior on S (analytically integrating over all B > 0)
    is a weighted mixture of Gamma distributions:

        p(S | n, m) ∝ e^{-S·t}  Σ_{j=0}^{n}  w_j · (S·t)^{n-j}

    where:
        w_j = C(n,j) · α^j · Γ(m+j+1) / (1+α)^{m+j+1}

    Kraft et al. (1991) corresponds to keeping only the j=0 term after
    replacing B with its point estimate α·m/t (i.e. treating B as known).

    As m → ∞ (background perfectly measured) the two methods converge.

    Parameters
    ----------
    n_obs      : int    — counts in source aperture
    n_bkg_raw  : int    — raw counts in background aperture (not scaled)
    area_ratio : float  — A_src / A_bkg
    t_eff      : float  — effective exposure time in seconds
    confidence : float  — one-sided CL (e.g. 0.9973 for 3σ)

    Returns
    -------
    S_upper : float  — upper limit on source count RATE in cts/s
    """
    n     = int(n_obs)
    m     = int(n_bkg_raw)
    alpha = float(area_ratio)
    t     = float(t_eff)

    j = np.arange(n + 1, dtype=float)

    # ---- log-weights (computed in log-space to handle large n, m) -----------
    log_binom  = gammaln(n + 1) - gammaln(j + 1) - gammaln(n - j + 1)
    log_alphaj = np.where(j == 0,
                          0.0,
                          j * np.log(alpha) if alpha > 0 else -np.inf)
    log_w = (log_binom
             + log_alphaj
             + gammaln(m + j + 1)
             - (m + j + 1) * np.log(1.0 + alpha))

    # ---- normalisation terms: w_j × Γ(n-j+1) -------------------------------
    log_nt    = log_w + gammaln(n - j + 1)
    log_shift = log_nt.max()                   # subtract max for numerical safety
    nt        = np.exp(log_nt - log_shift)     # relative normalisation terms
    norm      = nt.sum()

    # ---- CDF: P(S·t ≤ L) = Σ_j nt_j · gammainc(n-j+1, L) / norm -----------
    #  gammainc(a, x) = regularised lower incomplete gamma = P(a, x)
    a_arr = (n - j + 1).astype(int)            # a values (all ≥ 1)

    def cdf(L):
        if L <= 0.0:
            return 0.0
        vals = np.array([gammainc(int(a), L) for a in a_arr])
        return float(np.dot(nt, vals) / norm)

    if cdf(1e-300) >= confidence:
        return 0.0

    # ---- bracket: expand L_hi until cdf(L_hi) ≥ confidence -----------------
    L_hi = max(n, 1.0) + 10.0 * np.sqrt(max(n, 1.0)) + 50.0
    for _ in range(60):
        if cdf(L_hi) >= confidence:
            break
        L_hi *= 2.0

    L_ul = brentq(lambda L: cdf(L) - confidence, 0.0, L_hi,
                  xtol=1e-6, maxiter=500)

    return L_ul / t          # convert counts → cts/s


def gehrels_upper_limit(N_obs, B_scaled, confidence):
    """
    Gehrels 1986 closed-form Poisson upper limit.

        S_upper = (N_obs + 1 + sqrt(N_obs + 0.75) * z) - B_scaled

    where z is the one-sided Gaussian quantile for `confidence`.
    Slightly overestimates at very low N.  Useful as a cross-check.

    Parameters
    ----------
    N_obs      : int
    B_scaled   : float
    confidence : float

    Returns
    -------
    S_upper : float
    """
    z = sp_norm.ppf(confidence)
    return max(N_obs + 1.0 + np.sqrt(N_obs + 0.75) * z - B_scaled, 0.0)
