Estimates the log marginal likelihood (model evidence) from
    posterior samples using the Ellipsoidal Covering Marginal Likelihood
    Estimator (ECMLE). The method constructs a non-overlapping ellipsoidal
    cover of the Highest Posterior Density (HPD) region and uses the second
    half of posterior samples to form a ratio estimator of the marginal
    likelihood. The approach is applicable to multi-dimensional posterior
    distributions and requires only the posterior samples and their
    log-posterior values, without additional likelihood evaluations.
    See Naderi et al. (2025) <arXiv:2510.20617> for details.
