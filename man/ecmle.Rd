\name{ecmle}
\alias{ecmle}
\title{Estimate the log marginal likelihood with ECMLE}
\usage{
ecmle(
  post_samples,
  lps,
  log_post_fn,
  hpd_level = 0.75,
  subsample_frac = 0.05,
  r_max = NULL,
  bisect_tol = 1e-04,
  seed = NULL,
  verbose = FALSE
)
}
\arguments{
  \item{post_samples}{Numeric matrix of posterior draws; each row is one draw.}
  \item{lps}{Numeric vector of log unnormalized posterior values (log prior +
  log likelihood, up to an additive constant) evaluated at the rows of
  \code{post_samples}. Must be consistent with \code{log_post_fn};
  see the note below.}
  \item{log_post_fn}{A function that accepts a numeric vector of length \eqn{d}
  and returns a single finite numeric scalar equal to the log unnormalized
  posterior density at that point. Called repeatedly during ellipsoid
  construction, so it should be reasonably fast. Must be consistent
  with \code{lps}; see the note below.}
  \item{hpd_level}{Fraction in \code{(0, 1)} used to define the HPD region.
  Default \code{0.75}.}
  \item{subsample_frac}{Fraction of HPD points retained for candidate ellipsoid
  centers. Default \code{0.05}.}
  \item{r_max}{Upper bracketing radius used by the directional bisection solver.
  If \code{NULL}, estimated from the HPD points.}
  \item{bisect_tol}{Tolerance for the directional radius solver. Default
  \code{1e-4}.}
  \item{seed}{Optional integer seed used for HPD subsampling.}
  \item{verbose}{Logical; if \code{TRUE}, emits basic progress information.}
}
\note{
\strong{Consistency of \code{lps} and \code{log_post_fn} is critical.}
Both must evaluate the same log unnormalized posterior (log prior +
log likelihood) up to the same additive constant. If they differ by more
than a constant the marginal likelihood estimate will be incorrect.
}
\value{
An object of class \code{"ecmle"}, a list with components:
\describe{
  \item{\code{log_marginal_likelihood}}{Final log marginal likelihood estimate.}
  \item{\code{log_marginal_likelihood_iter}}{Running estimate over the second half-sample.}
  \item{\code{ellipsoids}}{List of fitted ellipsoids.}
  \item{\code{total_volume}}{Total ellipsoid volume.}
  \item{\code{n_ellipsoids}}{Number of ellipsoids.}
  \item{\code{points_in_ellipsoids}}{Number of evaluation points inside the ellipsoids.}
  \item{\code{n_samples}}{Number of evaluation points in the second half-sample.}
  \item{\code{coverage_rate}}{Fraction of evaluation points covered by the ellipsoids.}
  \item{\code{hpd_level}}{HPD level used in the fit.}
}
}
\description{
Computes an ECMLE estimate of the log marginal likelihood from posterior
draws, corresponding log unnormalized posterior evaluations, and a function
that evaluates the log unnormalized posterior density at arbitrary points.
}
\details{
The algorithm splits posterior draws into two halves. The first half is used
to define and pack ellipsoids inside a high-posterior-density region. The
second half is used to estimate the coverage-adjusted harmonic-mean quantity
that yields the final estimate.
}
\examples{
set.seed(1)
d <- 2
b <- rep(-1, d - 1)
a <- rep(0,  d - 1)

Y_bar <- rosenbrock_generate_data(
  theta_true = seq(1, 2.8, length.out = d),
  n = 20, b = b, a = a, sigma = 8
)
samps <- rosenbrock_exact_posterior(Y_bar, n = 20,
  b = b, a = a, sigma = 8, n_samples = 500L)
lps <- rosenbrock_log_post_vec(samps, Y_bar,
  n = 20, b = b, a = a, sigma = 8)
log_post_fn <- function(theta)
  rosenbrock_log_post(theta, Y_bar, n = 20, b = b, a = a, sigma = 8)

fit <- ecmle(samps, lps, log_post_fn, hpd_level = 0.75, seed = 1L)
fit
summary(fit)
plot_ecmle_2d(fit, post_samples = samps)
}