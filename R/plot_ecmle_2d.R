# ============================================================
#  2-D visualisation of ECMLE ellipsoids and posterior samples
# ============================================================

#' Draw a 2-D ellipse on an existing plot
#'
#' Low-level helper that traces a single ellipse (boundary of an ellipsoid in
#' 2-D) on the current graphics device via [graphics::lines()].  The ellipse
#' is defined by its center vector and covariance matrix \eqn{\Sigma}:
#' \deqn{\{x : (x - \mu)^\top \Sigma^{-1} (x - \mu) \le 1\}.}
#'
#' @param center  Numeric vector of length 2; ellipse center.
#' @param Sigma   \eqn{2 \times 2} symmetric positive-definite matrix.
#' @param npoints Integer; number of points on the ellipse boundary.
#'   Default 100.
#' @param col     Line colour.  Default \code{"red"}.
#' @param lwd     Line width.  Default \code{1.5}.
#' @param ...     Additional arguments passed to [graphics::lines()].
#'
#' @return Invisibly returns \code{NULL}.  Called for its side-effect.
#'
#' @examples
#' \donttest{
#' plot(0, 0, xlim = c(-3, 3), ylim = c(-3, 3), type = "n", asp = 1)
#' draw_ellipse_2d(center = c(0, 0), Sigma = diag(2))
#' }
#' @export
draw_ellipse_2d <- function(center, Sigma, npoints = 100L,
                             col = "red", lwd = 1.5, ...) {
  theta  <- seq(0, 2 * pi, length.out = npoints)
  circle <- cbind(cos(theta), sin(theta))
  # Cholesky of Sigma: ellipse = circle %*% chol(Sigma) + center
  chol_S <- chol(Sigma)
  pts    <- circle %*% chol_S
  pts[, 1] <- pts[, 1] + center[1]
  pts[, 2] <- pts[, 2] + center[2]
  graphics::lines(pts[, 1], pts[, 2], col = col, lwd = lwd, ...)
  invisible(NULL)
}


#' Plot ECMLE results in 2-D
#'
#' Visualises the output of [ecmle()] for a 2-dimensional posterior:
#' posterior samples are colour-coded according to whether they fall inside
#' any fitted ellipsoid, and all ellipsoid boundaries are drawn in red.
#'
#' For \eqn{d > 2} the function emits an informative message and returns
#' invisibly without plotting.
#'
#' @param fit           An object of class \code{"ecmle"} returned by [ecmle()].
#' @param post_samples  Numeric matrix of posterior draws (\eqn{N \times d});
#'   the same matrix that was passed to [ecmle()].  Only the **second**
#'   half-sample (rows \code{(floor(N/2)+1):N}) is used for the scatter, which
#'   matches the evaluation set used internally by ECMLE.
#' @param col_inside    Colour for samples inside the ellipsoids.
#'   Default \code{"navy"}.
#' @param col_outside   Colour for samples outside the ellipsoids.
#'   Default \code{"grey"}.
#' @param col_ellipse   Colour for the ellipsoid boundaries.
#'   Default \code{"red"}.
#' @param pch           Point character.  Default \code{19} (solid circle).
#' @param cex           Point size.  Default \code{0.3}.
#' @param lwd_ellipse   Line width for ellipse boundaries.  Default \code{1.5}.
#' @param npoints       Number of boundary points per ellipse.  Default
#'   \code{100}.
#' @param xlab          X-axis label.  Default \code{expression(theta[1])}.
#' @param ylab          Y-axis label.  Default \code{expression(theta[2])}.
#' @param main          Plot title.  When \code{NULL} (default) a title of
#'   the form \emph{"ECMLE: k ellipsoids"} is constructed automatically.
#' @param legend_pos    Position of the legend; passed to [graphics::legend()].
#'   Set to \code{NULL} to suppress the legend.  Default \code{"topright"}.
#' @param ...           Additional graphical parameters passed to
#'   [graphics::plot()].
#'
#' @return Invisibly returns \code{fit}.  Called for its side-effect.
#'
#' @examples
#' set.seed(42)
#' # Bivariate Gaussian toy example
#' post_samps <- cbind(rnorm(400), rnorm(400))
#' lps        <- apply(post_samps, 1, function(z) sum(stats::dnorm(z, log = TRUE)))
#' log_post   <- function(theta) sum(stats::dnorm(theta, log = TRUE))
#' fit <- ecmle(post_samps, lps, log_post, hpd_level = 0.75)
#' plot_ecmle_2d(fit, post_samples = post_samps)
#' @export
plot_ecmle_2d <- function(fit,
                           post_samples,
                           col_inside   = "navy",
                           col_outside  = "grey",
                           col_ellipse  = "red",
                           pch          = 19,
                           cex          = 0.3,
                           lwd_ellipse  = 1.5,
                           npoints      = 100L,
                           xlab         = expression(theta[1]),
                           ylab         = expression(theta[2]),
                           main         = NULL,
                           legend_pos   = "topright",
                           ...) {

  if (!inherits(fit, "ecmle")) {
    stop("`fit` must be an object of class 'ecmle' returned by ecmle().",
         call. = FALSE)
  }

  post_samples <- as.matrix(post_samples)
  d <- ncol(post_samples)

  if (d != 2L) {
    message("plot_ecmle_2d() is only available for 2-D posteriors (d = 2). ",
            "Skipping plot for d = ", d, ".")
    return(invisible(fit))
  }

  # Extract the second half-sample — the evaluation set used by ECMLE
  N    <- nrow(post_samples)
  half <- floor(N / 2)
  eval_pts <- post_samples[seq.int(half + 1L, N), , drop = FALSE]

  # Determine which evaluation points fall inside any ellipsoid
  n_ell <- fit$n_ellipsoids
  if (n_ell > 0L) {
    centers_mat  <- do.call(rbind, lapply(fit$ellipsoids, `[[`, "center"))
    Sigmas_inv   <- lapply(fit$ellipsoids, function(e) solve(e$Sigma))

    inside <- .is_point_in_any_ellipsoid(eval_pts, centers_mat, Sigmas_inv)
  } else {
    inside <- rep(FALSE, nrow(eval_pts))
  }

  pt_col <- ifelse(inside, col_inside, col_outside)

  # ---- base plot ---------------------------------------------------------
  if (is.null(main))
    main <- paste0("ECMLE: ", n_ell, " ellipsoid", if (n_ell != 1) "s")

  graphics::plot(
    eval_pts[, 1], eval_pts[, 2],
    col  = pt_col,
    pch  = pch,
    cex  = cex,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )

  # ---- ellipsoid boundaries ---------------------------------------------
  if (n_ell > 0L) {
    for (i in seq_len(n_ell)) {
      draw_ellipse_2d(
        center  = fit$ellipsoids[[i]]$center,
        Sigma   = fit$ellipsoids[[i]]$Sigma,
        npoints = npoints,
        col     = col_ellipse,
        lwd     = lwd_ellipse
      )
    }
  }

  # ---- legend ------------------------------------------------------------
  if (!is.null(legend_pos)) {
    graphics::legend(
      legend_pos,
      legend = c("Inside ellipsoids", "Outside ellipsoids",
                 "Ellipsoid boundaries"),
      col    = c(col_inside, col_outside, col_ellipse),
      pch    = c(pch, pch, NA),
      lty    = c(NA,  NA,  1),
      pt.cex = 0.8,
      cex    = 0.7,
      bty    = "n"
    )
  }

  invisible(fit)
}
