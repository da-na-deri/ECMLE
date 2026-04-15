# ============================================================
#  Pair-plot visualisation for posterior samples
# ============================================================

#' Pair plot of posterior samples
#'
#' Produces a lower-triangular scatter / density pair plot of a posterior
#' sample matrix.  Diagonal panels show marginal histograms; off-diagonal
#' lower panels show a 2-D pixel-density image using a blue colour ramp.
#' Upper panels are left blank.
#'
#' @param samples  Numeric matrix with one row per posterior draw and one
#'   column per parameter.  Column names are used as axis labels when
#'   \code{labels} is \code{NULL}.
#' @param pixs     Positive scalar controlling the pixel size of the 2-D
#'   density image.  Smaller values give finer resolution but require more
#'   computation.  Default is \code{1}.
#' @param labels   Optional character vector of axis labels (one per
#'   column).  When \code{NULL} (default) column names of \code{samples}
#'   are used; if \code{samples} has no column names, labels are set to
#'   \code{theta1}, \ldots, \code{thetad}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect of
#'   drawing a plot.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- 2
#' b <- rep(-1, d - 1)
#' a <- rep(0, d - 1)
#' Y_bar <- rosenbrock_generate_data(seq(1, 2.8, length.out = d),
#'   n = 20, b = b, a = a, sigma = 8)
#' samps <- rosenbrock_exact_posterior(Y_bar, n = 20, b = b, a = a,
#'   sigma = 8, n_samples = 2000L)
#' pair_plot(samps, pixs = 0.5)
#' }
#' @export
pair_plot <- function(samples, pixs = 1, labels = NULL) {
  
  samples <- as.matrix(samples)
  d <- ncol(samples)
  
  if (is.null(labels)) {
    labels <- if (!is.null(colnames(samples))) colnames(samples) else
      paste0("theta", seq_len(d))
  }
  
  .panel_hist <- function(x, ...) {
    usr <- graphics::par("usr")
    on.exit(graphics::par(usr = usr))
    graphics::par(usr = c(usr[1:2], 0, 1.5))
    h      <- graphics::hist(x, plot = FALSE)
    breaks <- h$breaks
    h2     <- graphics::hist(x, breaks = length(breaks) * 4, plot = FALSE)
    y      <- h2$counts
    y_norm <- if (max(y) > 0) y / max(y) else y
    graphics::rect(
      h2$breaks[-length(h2$breaks)], 0,
      h2$breaks[-1],                 y_norm,
      col = "blue4", border = NA, ...
    )
  }
  
  .panel_image <- function(x, y, pixs_val, ...) {
    xy      <- IDPmisc::NaRV.omit(data.frame(x = x, y = y))
    xy      <- sapply(xy, as.numeric)
    pixs_in <- (pixs_val / 10) / 2.54
    usr     <- graphics::par("usr")
    
    bx <- seq(usr[1], usr[2],
              length.out = round(graphics::par("pin")[1] / pixs_in) + 1)
    by <- seq(usr[3], usr[4],
              length.out = round(graphics::par("pin")[2] / pixs_in) + 1)
    
    zz   <- table(cut(xy[, 1], breaks = bx), cut(xy[, 2], breaks = by))
    zmax <- max(ceiling(max(zz)), 2)
    xx   <- 0.5 * (bx[-1] + bx[-length(bx)])
    yy   <- 0.5 * (by[-1] + by[-length(by)])
    
    graphics::image(
      x      = xx, y = yy, z = zz,
      col    = IDPmisc::IDPcolorRamp(zmax),
      breaks = seq(0.5, zmax + 1, 1),
      xaxs   = "r", yaxs = "r", add = TRUE
    )
    graphics::box()
    invisible(NULL)
  }
  
  graphics::pairs(
    samples,
    labels      = labels,
    lower.panel = function(...) {
      graphics::par(new = TRUE)
      args <- list(...)
      .panel_image(args[[1]], args[[2]], pixs_val = pixs)
    },
    diag.panel  = .panel_hist,
    upper.panel = NULL
  )
  
  invisible(NULL)
}