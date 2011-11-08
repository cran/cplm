\name{prior-struct}
\alias{igamma}
\alias{iwish}

\title{Constructor of the prior distribution
}
\description{Utilities functions to construct the prior distributions in the \code{cplm} family. 

}
\usage{
igamma(scale=0.01,shape=0.01)

iwish(df=3, scale=diag(1,df))
}

\arguments{
  \item{scale}{scale parameter in the inverse-Gamma distribution, or the scale matrix in the inverse-Wishart distribution.
}
  \item{shape}{shape parameter in the inverse-Gamma distribution.
}
  \item{df}{degree of freedom in the inverse-Wishart distribution.
}

}


\value{
  a list object. 
}


\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
See also \code{\link{bcpglmm}}.
}


\examples{
igamma(shape=0.001, scale=0.001)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

