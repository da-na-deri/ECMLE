#' @export
print.ecmle <- function(x, ...) {
  cat("<ecmle>\n")
  cat("  log marginal likelihood:", format(x$log_marginal_likelihood, digits = 6), "\n")
  cat("  ellipsoids:", x$n_ellipsoids, "\n")
  cat("  coverage rate:", format(x$coverage_rate, digits = 4), "\n")
  invisible(x)
}

#' Summarize an ECMLE fit
#'
#' Returns a compact summary of an `ecmle` fit.
#'
#' @param object An object returned by [ecmle()].
#' @param ... Unused.
#'
#' @return An object of class `"summary.ecmle"`.
#' @export
summary.ecmle <- function(object, ...) {
  out <- list(
    log_marginal_likelihood = object$log_marginal_likelihood,
    n_ellipsoids = object$n_ellipsoids,
    total_volume = object$total_volume,
    points_in_ellipsoids = object$points_in_ellipsoids,
    n_samples = object$n_samples,
    coverage_rate = object$coverage_rate,
    hpd_level = object$hpd_level
  )
  class(out) <- "summary.ecmle"
  out
}

#' @export
print.summary.ecmle <- function(x, ...) {
  cat("ECMLE summary\n")
  cat("  log marginal likelihood:", format(x$log_marginal_likelihood, digits = 6), "\n")
  cat("  ellipsoids:", x$n_ellipsoids, "\n")
  cat("  total volume:", format(x$total_volume, digits = 6), "\n")
  cat("  HPD level:", format(x$hpd_level, digits = 4), "\n")
  invisible(x)
}

#' Plot an ECMLE fit
#'
#' Plots the running log marginal likelihood estimate produced by [ecmle()].
#'
#' @param x An object returned by [ecmle()].
#' @param y Unused.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Plot title.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @export
plot.ecmle <- function(x,
                       y,
                       xlab = "Iteration",
                       ylab = "Running log marginal likelihood",
                       main = "ECMLE running estimate",
                       ...) {
  graphics::plot(
    seq_along(x$log_marginal_likelihood_iter),
    x$log_marginal_likelihood_iter,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  invisible(x)
}