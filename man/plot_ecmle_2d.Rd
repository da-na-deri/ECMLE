\name{plot_ecmle_2d}
\alias{plot_ecmle_2d}
\alias{draw_ellipse_2d}
\title{2-D visualisation of ECMLE ellipsoids and posterior samples}
\usage{
plot_ecmle_2d(
  fit,
  post_samples,
  col_inside = "navy",
  col_outside = "grey",
  col_ellipse = "red",
  pch = 19,
  cex = 0.3,
  lwd_ellipse = 1.5,
  npoints = 100L,
  xlab = expression(theta[1]),
  ylab = expression(theta[2]),
  main = NULL,
  legend_pos = "topright",
  ...
)

draw_ellipse_2d(center, Sigma, npoints = 100L, col = "red", lwd = 1.5, ...)
}
\arguments{
  \item{fit}{An object of class \code{"ecmle"} returned by
    \code{\link{ecmle}()}.}
  \item{post_samples}{Numeric matrix of posterior draws (N by d); the same
    matrix passed to \code{\link{ecmle}()}. Only the second half-sample is
    used for the scatter.}
  \item{col_inside}{Colour for samples inside at least one ellipsoid.
    Default \code{"navy"}.}
  \item{col_outside}{Colour for samples outside all ellipsoids.
    Default \code{"grey"}.}
  \item{col_ellipse}{Colour for ellipsoid boundaries. Default \code{"red"}.}
  \item{pch}{Point character; default \code{19} (solid circle).}
  \item{cex}{Point size; default \code{0.3}.}
  \item{lwd_ellipse}{Line width for ellipse boundaries; default \code{1.5}.}
  \item{npoints}{Integer; number of boundary points per ellipse.
    Default \code{100}.}
  \item{xlab}{X-axis label; default \code{expression(theta[1])}.}
  \item{ylab}{Y-axis label; default \code{expression(theta[2])}.}
  \item{main}{Plot title. When \code{NULL} (default) a title of the form
    "ECMLE: k ellipsoids" is constructed automatically.}
  \item{legend_pos}{Legend position string passed to
    \code{\link[graphics]{legend}()}, or \code{NULL} to suppress the legend.
    Default \code{"topright"}.}
  \item{center}{Numeric vector of length 2; ellipse center.}
  \item{Sigma}{2 by 2 symmetric positive-definite covariance matrix defining
    the ellipse boundary.}
  \item{col}{Line colour for \code{draw_ellipse_2d}; default \code{"red"}.}
  \item{lwd}{Line width for \code{draw_ellipse_2d}; default \code{1.5}.}
  \item{...}{Additional graphical parameters passed to
    \code{\link[graphics]{plot}()} or \code{\link[graphics]{lines}()}.}
}
\value{
\code{plot_ecmle_2d} returns \code{fit} invisibly.
\code{draw_ellipse_2d} returns \code{NULL} invisibly.
Both are called for their side-effect (drawing on the active graphics device).
}
\description{
\code{plot_ecmle_2d} visualises the output of \code{\link{ecmle}()} for a
two-dimensional posterior. The second-half evaluation samples are colour-coded
by ellipsoid membership (navy = inside, grey = outside) and the fitted
ellipsoid boundaries are overlaid in red. For \eqn{d > 2} the function emits
an informative message and returns \code{fit} invisibly without plotting.

\code{draw_ellipse_2d} is the low-level helper that traces a single ellipse
boundary on the current graphics device. It can also be called directly to
add ellipses to an existing plot.
}
\examples{
set.seed(42)
post_samps <- cbind(rnorm(400), rnorm(400))
lps <- apply(post_samps, 1, function(z) sum(dnorm(z, log = TRUE)))
log_post_fn <- function(theta) sum(dnorm(theta, log = TRUE))

fit <- ecmle(post_samps, lps, log_post_fn, hpd_level = 0.75, seed = 1L)
plot_ecmle_2d(fit, post_samples = post_samps)

# draw_ellipse_2d standalone usage
plot(0, 0, xlim = c(-3, 3), ylim = c(-3, 3), type = "n", asp = 1,
     xlab = "", ylab = "")
draw_ellipse_2d(center = c(0, 0),
                Sigma  = matrix(c(1, 0.4, 0.4, 1), 2, 2))
}
}
