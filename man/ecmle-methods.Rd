\name{ecmle-methods}
\alias{print.ecmle}
\alias{summary.ecmle}
\alias{print.summary.ecmle}
\alias{plot.ecmle}
\title{Methods for ECMLE objects}
\usage{
\method{print}{ecmle}(x, ...)

\method{summary}{ecmle}(object, ...)

\method{print}{summary.ecmle}(x, ...)

\method{plot}{ecmle}(
  x,
  y,
  xlab = "Iteration",
  ylab = "Running log marginal likelihood",
  main = "ECMLE running estimate",
  ...
)
}
\arguments{
  \item{x}{An object of class \code{"ecmle"} or \code{"summary.ecmle"}.}
  \item{object}{An object of class \code{"ecmle"}.}
  \item{y}{Unused.}
  \item{xlab}{X-axis label for the plot.}
  \item{ylab}{Y-axis label for the plot.}
  \item{main}{Plot title.}
  \item{...}{Unused or additional graphical parameters.}
}
\value{
\code{summary.ecmle()} returns an object of class \code{"summary.ecmle"}.
The print and plot methods return their input invisibly.
}
\description{
S3 methods for inspecting ECMLE fits.
}
