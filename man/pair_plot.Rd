\name{pair_plot}
\alias{pair_plot}
\title{Pair plot of posterior samples}
\usage{
pair_plot(samples, pixs = 1, labels = NULL)
}
\arguments{
  \item{samples}{Numeric matrix with one row per posterior draw and one column
    per parameter. Column names are used as axis labels when \code{labels} is
    \code{NULL}.}
  \item{pixs}{Positive scalar controlling the pixel size of the 2-D density
    image. Smaller values produce finer resolution at higher computational
    cost. Default \code{1}.}
  \item{labels}{Optional character vector of length \code{ncol(samples)}.
    When \code{NULL} (default), column names of \code{samples} are used; if
    there are no column names, labels default to
    \code{theta1}, \code{theta2}, etc.}
}
\value{
Invisibly returns \code{NULL}. Called for its side-effect (drawing a plot).
}
\description{
Produces a lower-triangular scatter/density pair plot of a posterior sample
matrix. Diagonal panels display scaled marginal histograms. Lower off-diagonal
panels display a 2-D pixel-density image using a blue colour ramp. Upper
panels are left blank.
}
\examples{
\donttest{
  set.seed(1)
  d <- 2
  b <- rep(-1, d - 1)
  a <- rep(0, d - 1)
  Y_bar <- rosenbrock_generate_data(seq(1, 2.8, length.out = d),
                                    n = 20, b = b, a = a, sigma = 8)
  samps <- rosenbrock_exact_posterior(Y_bar, n = 20, b = b, a = a,
                                      sigma = 8, n_samples = 2000L)
  pair_plot(samps, pixs = 0.5)
}
}

