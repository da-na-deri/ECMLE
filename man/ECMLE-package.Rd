\name{ECMLE-package}
\alias{ECMLE-package}
\docType{package}
\title{ECMLE: Elliptical Covering Marginal Likelihood Estimator}
\description{
Implements Elliptical Covering Marginal Likelihood Estimator
    (ECMLE), a geometric method for approximating marginal likelihood from
    posterior draws and log-posterior evaluations. The approach constructs a
    collection of non-overlapping ellipsoids within a high-posterior-density
    region, computes the covered volume, and combines this with posterior sample
    coverage to estimate model evidence. The method is designed to stabilize
    harmonic-mean-based evidence approximation and can be applied in multimodal
    settings. The methodology is described in Naderi et al. (2025)
    <doi:10.48550/arXiv.2510.20617>.
}
\details{
The package provides the main estimator \code{ecmle()} with print, summary,
and plot methods, visualisation functions \code{plot_ecmle_2d()} and
\code{draw_ellipse_2d()}, a pair plot utility \code{pair_plot()}, and
Rosenbrock (banana) posterior benchmark functions for testing and reproducible
experiments.
}
\author{
Dana Naderi [aut, cre],
Christian P Robert [aut],
Kaniav Kamary [aut],
Darren Wraith [aut]
}
\seealso{\code{\link{ecmle}}}
\keyword{package}
