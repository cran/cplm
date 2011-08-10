\name{cpglm}
\alias{cpglm}
\title{Compound Poisson Generalized Linear Model
}
\description{This function implements the Monte Carlo EM algorithm for fitting generalized linear models with Tweedie compound Poisson distribution. 
}
\usage{
cpglm(formula, link = "log", data, weights, subset, 
  na.action, betastart = NULL,phistart = NULL, pstart = NULL, 
  offset,  contrasts = NULL, control = list(),
  ...)
}

\arguments{
  \item{formula}{an object of class \code{formula}. See also in \code{glm}.
}
  \item{link}{a specification for the model link function. This can only be a literal character string, and must be either "log", "identity", "sqrt" or "inverse". 
}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. 
}
  \item{weights}{an optional vector of weights. Should be \code{NULL} or a numeric vector.
}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of options, and is \code{na.fail} if that is unset.  Another possible value is \code{NULL}, no action. Value \code{na.exclude} can be useful.
}
  \item{betastart}{starting values for the vector of mean parameters.
}
  \item{phistart}{starting value for the scale parameter.
}
  \item{pstart}{starting value for the index parameter. Must be between 1 and 2. 
}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. 
}
  \item{contrasts}{an optional list. See the \code{contrasts.arg}.
}
  \item{control}{a list of parameters for controlling the fitting process. See 'Details' below. 
} 
  \item{\dots}{ not used now. 
}

}

\details{
The Tweedie compound Poisson distribution is a subclass of the exponential dispersion family with a power variance function, where the value of the power index lies in the interval (1,2). The distribution is composed of a probability mass at the origin accompanied by a skewed continuous distribution on the positive values. Despite its ability to handle continuous zero-inflated data, the density of the compound Poisson distribution is not analytically tractable. As a result, full maximum likelihood inference remains challenging when the power index is unknown. The \code{cpglm} function implements a Monte Carlo EM [MCEM] algorithm that yields maximum likelihood estimates for all parameters in the model, in particular,  the index parameter. This approach is in light of the fact that the compound Poisson distribution involves an unobserved Poisson process. The observed data is thus augmented by the latent Poisson variable so that one could work on the joint likelihood of the complete data and invoke the MCEM algorithm to locate the maximum likelihood estimations.  

Specifically, a Monte Carlo method is implemented in the E-step, where samples of the unobserved Poisson process is simulated from its posterior via rejection sampling (using a zero-truncated Poisson as a proposal), and the required expectation in the E-step is approximated by a Monte Carlo estimate. The M-step involves conditional component-wise maximization: scoring algorithm is used for the mean parameters, constrained optimization for the index parameter and close-form solution exists for the scale parameter. The E-step and the M-step are iterated until the following stopping rule reached:
\deqn{\max_i \frac{|\theta_i^{(r+1)}-\theta_i^{(r)}|}{\theta_i^{(r)}+\epsilon_1}<\epsilon_2,
}
where \eqn{\theta_i} is the i_th parameter, and \eqn{\epsilon_1} and \eqn{\epsilon_2} are predefined constants. See the description for the \code{control} parameters below. 

The MCEM algorithm enables one to work one the more tractable joint likelihood, thus simplifying the estimation problem down to posterior simulation and block/univariate optimizations. However, it also has some overheads:
\itemize{
  \item The simulation of the latent variable could be computationally demanding when a large number of samples is needed;
  \item The inherent Monte Carlo error must be accounted for;
  \item The EM algorithm could have slow convergence rate.
}

The computational demand can be largely reduced by supplementing the algorithm with importance sampling after it converges to the neighborhood of the MLEs. The importance sampling retains the latest simulated samples for use in all later iterations, thus eliminating the need for simulating new samples. 

To account for the Monte Carlo error, many authors have suggests increasing the sample size as the algorithm proceeds. Booth and Hobert (1999) derives a Normal approximation for the Monte Carlo error, which is further used to construct a (1-\code{alpha}) confidence ellipsoid at each iteration. If the old parameters lie in this ellipsoid, then the sample size is increase as \eqn{m + m/k}, where \eqn{m} is the current sample size, and \eqn{k} is a predefined constant. This method is implemented when \code{fixed.size} in the control argument is \code{FALSE}. See the description for the \code{control} parameters below.    



The \code{control} argument is a list that can supply any of the following components:

\describe{
  \item{\code{init.size}}{initial number of samples used to generate the latent Poisson variable. Default is 100.}
\item{\code{sample.iter}}{the iteration number after which importance sampling is used where the samples in iteration \code{sample.iter}-1 will be recycled in all later iterations. The could improve the speed in the E-step, but should only be used when the algorithm has converged to the neighborhood of the MLE. The default value is 50. To suppress importance sampling, set this number to be large enough.}
\item{\code{maxit}}{maximum number of iterations allowed in the algorithm. The default value is 200.}
\item{\code{epsilon1}}{a constant used in the stopping rule. The default value is 0.001.}
\item{\code{epsilon2}}{a constant used in the stopping rule. The default value is 0.0001.}
\item{\code{fixed.size}}{a logical value indicating whether the sample size in each iteration is to be increased. If \code{FALSE}, approximate variance is computed in each iteration and the approximate (1-\code{alpha}) confidence ellipsoid is constructed. If the previously simulated value falls in this ellipsoid, then the sample size is increased by 1/\code{k} of the current sample size.  }
\item{\code{alpha}}{the confidence level used in the above calculation of confidence ellipsoid. The default value is 0.25.}
\item{\code{k}}{a constant used in the calculation of new sample size. The default value is 5.}
\item{\code{max.size}}{the maximum sample size used in each iteration. When \code{fixed.size} is \code{FALSE}, the sample size in later iterations could be increasing without bound and as a result, the algorithm could be impractically slow. \code{max.size} sets the limit of the required sample size to avoid such slow computation. The default value is 10000. }
\item{\code{bound.p}}{a vector of lower and upper bound for the index parameter p. The default is \code{c(1.01,1.99)}.}
\item{\code{verbose}}{if \code{TRUE}, tracing information on the progress of the fitting is produced. }
}

}

\value{
  \code{cpglm} returns an object of class \code{cpglm}. See \code{\link{cpglm-class}}. 
}

\references{
\cite{Booth, J. G., and Hobert, J. P. (1999). Maximizing Generalized Linear Mixed Model Likelihoods with an Automated
Monte Carlo EM Algorithm.  \emph{Journal of the Royal Statistical Society Series B}, 61, 265-285.}

\cite{McCulloch, C. E. (1997). Maximum likelihood algorithms for generalized linear mixed
models. \emph{Journal of the American Statistical Association}, 92, 162-170.}

\cite{Zhang Y. (2011). A Latent Variable Approach for Statistical Inference  in Tweedie Compound Poisson Linear Models, \emph{preprint}, \url{http://actuaryzhang.com/publication/bayesianTweedie.pdf}}.

}

\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
See also in \code{\link{cpglm-class}}.
}


\examples{

\dontrun{
set.seed(10)
fit <- cpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
            control=list(init.size=10,fixed.size=FALSE,
                max.size=3000,sample.iter=60))
summary(fit)                
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models}

