# ECMLE

`ECMLE` is an R package for the **Elliptical Covering Marginal Likelihood Estimator**.
It implements a geometric marginal-likelihood estimator based on posterior draws,
log-posterior evaluations, and ellipsoidal coverings of a high-posterior-density region.

The method is described in:

> Naderi et al. (2025). *Approximating Evidence via Bounded Harmonic Means.*
> <https://doi.org/10.48550/arXiv.2510.20617>

## Installation

```r
install.packages("ECMLE")
```

## Main function

```r
ecmle(post_samples, lps, log_post_fn, hpd_level = 0.75)
```

### Arguments

- `post_samples`: matrix of posterior draws, one row per draw.
- `lps`: numeric vector of log-posterior values at `post_samples`.
- `log_post_fn`: function returning the log unnormalized posterior density (log prior + log likelihood) at a parameter vector.
- `hpd_level`: HPD fraction used to define the packing region.

## Minimal example

```r
set.seed(1)
post_samples <- cbind(rnorm(400), rnorm(400))
lps <- apply(post_samples, 1, function(z) sum(dnorm(z, log = TRUE)))
log_post_fn <- function(theta) sum(dnorm(theta, log = TRUE))

fit <- ecmle(
  post_samples = post_samples,
  lps          = lps,
  log_post_fn  = log_post_fn,
  hpd_level    = 0.75
)
print(fit)
plot(fit)
```

## Package contents

- `ecmle()`: estimates the log marginal likelihood.
- `summary()`: returns a compact summary of an `ecmle` fit.
- `plot()`: shows the running estimate across the evaluation half-sample.
- `plot_ecmle_2d()`, `draw_ellipse_2d()`: 2-D visualisation of fitted ellipsoids.
- `pair_plot()`: lower-triangular pair plot for posterior sample matrices.
- `rosenbrock_generate_data()`, `rosenbrock_exact_posterior()`,
  `rosenbrock_log_post()`, `rosenbrock_log_post_vec()`:
  Rosenbrock (banana) posterior benchmark functions.

## License

GPL-3