\name{rosenbrock}
\alias{rosenbrock_generate_data}
\alias{rosenbrock_exact_posterior}
\alias{rosenbrock_log_post}
\alias{rosenbrock_log_post_vec}
\title{Rosenbrock (banana) posterior}
\usage{
rosenbrock_generate_data(theta_true, n, b, a, sigma)

rosenbrock_exact_posterior(Y_bar, n, b, a, sigma, n_samples = 10000L)

rosenbrock_log_post(theta, Y_bar, n, b, a, sigma)

rosenbrock_log_post_vec(theta_matrix, Y_bar, n, b, a, sigma)
}
\arguments{
  \item{theta_true}{Numeric vector of length \eqn{d}; true parameter values.}
  \item{n}{Positive integer; number of observations used to form the sufficient statistic.}
  \item{b}{Numeric vector of length \eqn{d-1}; curvature (banana) parameters.}
  \item{a}{Numeric vector of length \eqn{d-1}; shift parameters.}
  \item{sigma}{Positive scalar; observation standard deviation.}
  \item{Y_bar}{Numeric vector of length \eqn{d}; sufficient statistic, typically
    the output of \code{rosenbrock_generate_data()}.}
  \item{n_samples}{Positive integer; number of posterior draws to generate.
    Default \code{10000}.}
  \item{theta}{Numeric vector of length \eqn{d}; a single parameter value.}
  \item{theta_matrix}{Numeric matrix with \eqn{d} columns; one parameter
    vector per row.}
}
\value{
\describe{
  \item{\code{rosenbrock_generate_data}}{Numeric vector of length \eqn{d}.}
  \item{\code{rosenbrock_exact_posterior}}{Numeric matrix of dimensions
    \code{n_samples} by \eqn{d}, with column names
    \code{theta1}, ..., \code{thetad}.}
  \item{\code{rosenbrock_log_post}}{A single numeric value.}
  \item{\code{rosenbrock_log_post_vec}}{Numeric vector of length
    \code{nrow(theta_matrix)}.}
}
}
\description{
Functions for the Rosenbrock (banana-shaped) posterior, providing data
generation, exact posterior sampling, and log-posterior evaluation.

\code{rosenbrock_generate_data} draws one sufficient statistic vector.
\code{rosenbrock_exact_posterior} draws exact posterior samples.
\code{rosenbrock_log_post} evaluates the log-posterior at a single parameter
vector (suitable as the \code{log_post_fn} argument to \code{\link{ecmle}()}).
\code{rosenbrock_log_post_vec} is the vectorised version over a sample matrix.
}
\examples{
set.seed(1)
d <- 2
b <- rep(-1, d - 1)
a <- rep(0,  d - 1)

# 1. Generate one dataset
Y_bar <- rosenbrock_generate_data(
  theta_true = seq(1, 2.8, length.out = d),
  n = 20, b = b, a = a, sigma = 8
)

# 2. Draw exact posterior samples
samps <- rosenbrock_exact_posterior(Y_bar, n = 20,
  b = b, a = a, sigma = 8, n_samples = 500L)

# 3. Evaluate the log-posterior
lps <- rosenbrock_log_post_vec(samps, Y_bar,
  n = 20, b = b, a = a, sigma = 8)

# 4. Build a log_post_fn closure for ecmle()
log_post_fn <- function(theta)
  rosenbrock_log_post(theta, Y_bar, n = 20, b = b, a = a, sigma = 8)

fit <- ecmle(samps, lps, log_post_fn, hpd_level = 0.75, seed = 1L)
fit
}