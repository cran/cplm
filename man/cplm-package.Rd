\name{cplm-package}
\alias{cplm-package}
\alias{cplm}
\docType{package}
\title{
  Tweedie compound Poisson linear models 
}
\description{This R package provides various methods for fitting Tweedie compound Poisson linear models. The Tweedie compound Poisson distribution is a subclass of the exponential dispersion family with a power variance function, where the value of the power index lies in the interval (1,2). The distribution is composed of a probability mass at the origin accompanied by a skewed continuous distribution on the positive values, and has been widely used in applied statistical fields for modeling continuous data with exact zeros. However, it is well known that the density function is not analytically tractable, and in order to carry out maximum likelihood approximations, we can either rely on numerical approximations of the density (available in the \code{tweedie} package), or introduce latent variables and work on the complete data to avoid direct density approximations.  Both methods are available in the package.  

The following features of the package may be of special interest to the users:

\enumerate{
\item{All methods available in the package enable the index parameter to be estimated from the data.}
\item{The compound Poisson generalized linear model handles large data set using the bounded memory regression facility in \code{biglm}.}
\item{For mixed-effect models, we provide Laplace and adaptive Gauss-Hermit quadrature methods. These likelihood approximation methods generally have better asymptotic properties than the Penalized Quasi-likelihood approach.}
\item{An convenient interface is offered to fit compound Poisson additive models using the mixed-effect model structure.}
\item{Self-tuned Markov Chain Monte Carlo procedures are available for both the compound Poisson generalized linear model and the mixed-effect model.}
}

More information is available on the hosting web site of the project \url{http://code.google.com/p/cplm/}
}

\author{
Wayne Zhang <actuary_zhang@hotmail.com>
}
\references{
\cite{ Dunn, P.K. and Smyth, G.K. (2005). Series evaluation of Tweedie exponential dispersion models densities. \emph{Statistics and Computing}, 15, 267-280.}

}


\keyword{ package }


