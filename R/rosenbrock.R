# ============================================================
#  Rosenbrock (banana) posterior
# ============================================================

#' Generate a single Rosenbrock sufficient statistic
#'
#' Draws one realisation of the sufficient statistic \eqn{\bar{Y}} from the
#' Rosenbrock (banana-shaped) observation model with flat prior.
#'
#' The observation model is:
#' \deqn{
#'   \bar{Y}_1 \sim N(\theta_1,\; \sigma^2/n), \quad
#'   \bar{Y}_j \sim N\!\bigl(\theta_j + b_{j-1}(\theta_{j-1}^2 - a_{j-1}),\;
#'                            \sigma^2/n\bigr), \; j = 2,\ldots,d.
#' }
#'
#' @param theta_true Numeric vector of length \eqn{d}; true parameter values.
#' @param n          Integer; number of observations used to form \eqn{\bar{Y}}.
#' @param b          Numeric vector of length \eqn{d-1}; curvature parameters.
#' @param a          Numeric vector of length \eqn{d-1}; shift parameters.
#' @param sigma      Positive scalar; observation standard deviation.
#'
#' @return Numeric vector of length \eqn{d} containing \eqn{\bar{Y}}.
#'
#' @examples
#' set.seed(1)
#' d <- 2
#' Y_bar <- rosenbrock_generate_data(
#'   theta_true = seq(1, 2.8, length.out = d),
#'   n = 20, b = rep(-1, d - 1), a = rep(0, d - 1), sigma = 8
#' )
#' @export
rosenbrock_generate_data <- function(theta_true, n, b, a, sigma) {
  d <- length(theta_true)
  if (length(b) != d - 1 || length(a) != d - 1)
    stop("`b` and `a` must each have length `length(theta_true) - 1`.",
         call. = FALSE)
  Y_bar    <- numeric(d)
  sd_j     <- sigma / sqrt(n)
  Y_bar[1] <- stats::rnorm(1, mean = theta_true[1], sd = sd_j)
  if (d > 1) {
    for (j in 2:d) {
      mu_j     <- theta_true[j] + b[j - 1] * (theta_true[j - 1]^2 - a[j - 1])
      Y_bar[j] <- stats::rnorm(1, mean = mu_j, sd = sd_j)
    }
  }
  Y_bar
}


#' Draw exact posterior samples from the Rosenbrock model
#'
#' Given the sufficient statistic \eqn{\bar{Y}}, draws \code{n_samples} exact
#' posterior samples using sequential Gaussian sampling (the posterior is a
#' product of univariate Gaussians given the preceding component).
#'
#' @param Y_bar     Numeric vector of length \eqn{d}; sufficient statistic
#'   produced by [rosenbrock_generate_data()].
#' @param n         Integer; number of observations.
#' @param b         Numeric vector of length \eqn{d-1}; curvature parameters.
#' @param a         Numeric vector of length \eqn{d-1}; shift parameters.
#' @param sigma     Positive scalar; observation standard deviation.
#' @param n_samples Integer; number of posterior draws to generate.
#'
#' @return Numeric matrix of dimensions \code{n_samples} \eqn{\times} \eqn{d},
#'   with column names \code{theta1}, \ldots, \code{thetad}.
#'
#' @examples
#' set.seed(1)
#' d <- 2
#' Y_bar <- rosenbrock_generate_data(
#'   theta_true = seq(1, 2.8, length.out = d),
#'   n = 20, b = rep(-1, d - 1), a = rep(0, d - 1), sigma = 8
#' )
#' samps <- rosenbrock_exact_posterior(Y_bar, n = 20,
#'   b = rep(-1, d - 1), a = rep(0, d - 1),
#'   sigma = 8, n_samples = 500L)
#' head(samps)
#' @export
rosenbrock_exact_posterior <- function(Y_bar, n, b, a, sigma,
                                       n_samples = 1e4L) {
  d    <- length(Y_bar)
  sd_j <- sigma / sqrt(n)
  samples <- matrix(NA_real_, nrow = n_samples, ncol = d,
                    dimnames = list(NULL, paste0("theta", seq_len(d))))
  samples[, 1] <- stats::rnorm(n_samples, mean = Y_bar[1], sd = sd_j)
  if (d > 1) {
    for (j in 2:d) {
      mu_j_vec     <- Y_bar[j] - b[j - 1] * (samples[, j - 1]^2 - a[j - 1])
      samples[, j] <- stats::rnorm(n_samples, mean = mu_j_vec, sd = sd_j)
    }
  }
  samples
}


#' Log-posterior of the Rosenbrock model (flat prior)
#'
#' Evaluates \eqn{\log p(\theta \mid \bar{Y})} (up to an additive constant)
#' for a single parameter vector \eqn{\theta}.
#'
#' @param theta Numeric vector of length \eqn{d}.
#' @param Y_bar Numeric vector of length \eqn{d}; sufficient statistic.
#' @param n     Integer; number of observations.
#' @param b     Numeric vector of length \eqn{d-1}; curvature parameters.
#' @param a     Numeric vector of length \eqn{d-1}; shift parameters.
#' @param sigma Positive scalar; observation standard deviation.
#'
#' @return A single numeric value.
#'
#' @examples
#' d <- 2
#' rosenbrock_log_post(
#'   theta = c(1.5, 2.0),
#'   Y_bar = c(1.4, 2.1),
#'   n = 20, b = rep(-1, d - 1), a = rep(0, d - 1), sigma = 8
#' )
#' @export
rosenbrock_log_post <- function(theta, Y_bar, n, b, a, sigma) {
  d     <- length(theta)
  sd_j  <- sigma / sqrt(n)
  mu    <- numeric(d)
  mu[1] <- theta[1]
  if (d > 1) {
    for (j in 2:d)
      mu[j] <- theta[j] + b[j - 1] * (theta[j - 1]^2 - a[j - 1])
  }
  sum(stats::dnorm(Y_bar, mean = mu, sd = sd_j, log = TRUE))
}


#' Vectorised log-posterior for Rosenbrock model
#'
#' Same as [rosenbrock_log_post()] but operates on a matrix of parameter
#' vectors (one row per parameter vector) for efficiency.
#'
#' @param theta_matrix Numeric matrix; each row is one parameter vector.
#' @param Y_bar        Numeric vector of length \eqn{d}.
#' @param n            Integer; number of observations.
#' @param b            Numeric vector of length \eqn{d-1}.
#' @param a            Numeric vector of length \eqn{d-1}.
#' @param sigma        Positive scalar.
#'
#' @return Numeric vector of log-posterior values, one per row of
#'   \code{theta_matrix}.
#'
#' @examples
#' set.seed(1)
#' d <- 2
#' samps <- matrix(rnorm(200), ncol = d)
#' Y_bar <- c(1.4, 2.1)
#' lp    <- rosenbrock_log_post_vec(samps, Y_bar,
#'   n = 20, b = rep(-1, d - 1), a = rep(0, d - 1), sigma = 8)
#' head(lp)
#' @export
rosenbrock_log_post_vec <- function(theta_matrix, Y_bar, n, b, a, sigma) {
  n_samp      <- nrow(theta_matrix)
  d           <- ncol(theta_matrix)
  sd_j        <- sigma / sqrt(n)
  mu_mat      <- matrix(NA_real_, nrow = n_samp, ncol = d)
  mu_mat[, 1] <- theta_matrix[, 1]
  if (d > 1) {
    for (j in 2:d)
      mu_mat[, j] <- theta_matrix[, j] +
        b[j - 1] * (theta_matrix[, j - 1]^2 - a[j - 1])
  }
  rowSums(stats::dnorm(
    matrix(Y_bar, nrow = n_samp, ncol = d, byrow = TRUE),
    mean = mu_mat, sd = sd_j, log = TRUE
  ))
}