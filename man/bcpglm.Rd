\name{bcpglm}
\alias{bcpglm}

\title{
Bayesian Compound Poisson Generalized Linear Model
}
\description{
This function handles Bayesian Tweedie compound Poisson generalized linear models and generates posterior simulations using Markov Chain Monte Carlo methods.
}
\usage{
bcpglm(formula, link = "log", data, inits = NULL, 
  weights, offset, subset, na.action, contrasts = NULL, 
  n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2), 
  n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), 
  n.sims = 1000, n.report = 1000, prior.beta.mean, prior.beta.var, 
  bound.phi=100, bound.p = c(1.01, 1.99),method="dtweedie", 
  tune.iter=4000, n.tune=10, tune.weight=0.25,...)
}

\arguments{
 \item{formula}{an object of class \code{formula}. See also in \code{\link[stats]{glm}}.
}
  \item{link}{a specification for the model link function. This can be either a literal character string or a numeric number. If it is a character string, it must be one of "log", "identity", "sqrt" or "inverse". If it is numeric, it is the same as the \code{link.power} argument in the \code{\link[statmod]{tweedie}} function. The default is \code{link="log"}.
}
 \item{inits}{a list of initial values to be used for each chain. It must be of length \code{n.chains}. For each element, it is a named list with three components 'beta', 'phi' and 'p'. If not supplied, the function will generate initial values automatically, which are based on a GLM with the supplied model structure. 
}
  \item{weights}{an optional vector of weights. Should be \code{NULL} or a numeric vector. When it is numeric, it must be positive. Zero weights are not allowed in \code{bcpglm}. 
}
 
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. 
}
 
  \item{data, subset, na.action, contrasts}{same as those in \code{\link[stats]{glm}} and \code{\link{cpglm}}. 
}
  \item{n.chains}{an integer indicating the number of Markov chains (default: \code{3}).  
}
  \item{n.iter}{
number of total iterations per chain (including burn in; default: \code{2000})
}
  \item{n.burnin}{
length of burn in, i.e. number of iterations to discard at the beginning. Default
is \code{n.iter/2}, that is, discarding the first half of the simulations.
}
  \item{n.thin}{
thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and
computation time if n.iter is large. Default is \code{max(1, floor(n.chains*(n.iter-n.burnin) / 1000))} which will only thin if there are at
least 2000 simulations.
}
  \item{n.sims}{
The approximate number of simulations to keep after thinning.
}
  \item{n.report}{if greater than zero, fitting information will be printed out after each \code{n.report} iterations. 
}
  \item{prior.beta.mean}{a vector of prior means for the model coefficients. Default is a vector of zeros. 
}
  \item{prior.beta.var}{
a vector of prior variance for the model coefficients. Default is a vector of \code{10000}'s.
}

\item{bound.phi}{upper bound of the uniform prior for the dispersion parameter phi. The default is \code{100}. 
}

  \item{bound.p}{a vector of lower and upper bound for the index parameter \eqn{p}. The default is \code{c(1.01,1.99)}.
}

\item{method}{determines how the MCMC is implemented. If \code{method="dtweedie"} (the default), then full conditionals are computed using numerical methods to approximate the tweedie density. If \code{method="latent"}, then the latent Poisson variables are introduced into the MCMC so that the need for direct evaluation of the tweedie density is eliminated. The latter method is also much faster, but may fail in certain circumstances.  
}

\item{tune.iter}{number of iterations used for tuning the variance for the Normal proposal distribution. These iterations will not be used in 
the final output. Default is \code{4000}, and set it to be zero if the tuning process is not desired.  
}
\item{n.tune}{if positive, the \code{tune.iter} iterations will be divided into \code{n.tune} loops. Default is \code{10}.
}
\item{tune.weight}{the weight to be given to the old covariance matrix as opposed to the covariance matrix estimated from simulated samples. Default is \code{0.25}. }

  \item{\dots}{ not used.
}
}
\details{
The function fits the Tweedie compound Poisson generalized linear model within the Bayesian framework, where posterior distributions of model parameters are estimated using Monte Carlo Markov Chains [MCMC]. Two methods are provided here: one relies on the capability to approximate the Tweedie density (i.e., the normalizing constant involving \eqn{\phi} and \eqn{p}), the other relies on the use of latent Poisson variables to avoid direct evaluation of the Tweedie density (\code{method="latent"}). The first method is more general while the second method is much faster when it works.   

Since this is a Bayesian model, prior distributions have to be specified for all parameters in the model. Here, Normal distributions are used for the model coefficients (\eqn{\beta}), a Uniform distribution for the  dispersion parameter (\eqn{\phi}) and a Uniform distribution for the index parameter (\eqn{p}). Prior means and variances for the model coefficients can be supplied using the argument \code{prior.beta.mean} and \code{prior.beta.var}, respectively. The prior distribution for \eqn{\phi} is uniform on (0, \code{bound.phi}). And the bounds of the Uniform for \eqn{p} can be specified in the argument \code{bound.p}. See details in Section 'Arguments'. 

In implementing the MCMC, a Gibbs sampler is constructed as follows:
\enumerate{
  \item{If \code{method="latent"}, simulation of the latent Poisson variables relies on rejection sampling, where zero-truncated Poisson distributions are used as proposals. 
}
  \item{Model coefficients are updated block-wise using the Metropolis-Hastings algorithm, where a random-walk multivariate Normal proposal distribution is used. The variance-covariance matrix used in the multivariate Normal is the variance-covariance matrix estimated from a preliminary \code{cpglm} fit with the specified model structure.}
  \item{Updates of the index parameter and the dispersion parameter rely on univariate Metropolis-Hastings, where a truncated Normal distribution is used as a proposal. }
  }

Before the MCMC, there is a tuning process where the proposal variance from the (truncated) Normal distribution is updated according to the sample variance computed from the previous simulations. The goal is to make the acceptance rate at about 25\% for \eqn{\beta}, and about 45\% for \eqn{p} and \eqn{\phi}. The argument \code{tune.iter} determines how many iterations are used for the tuning process, and \code{n.tune} determines how many loops the \code{tune.iter} iterations should be divided into. These iterations will not be used in the final output. 

The simulated values of all model parameters are stored in the slot \code{sims.list} of the returned \code{bcpglm} object. It is a list of \code{n.chains} matrices and each matrix has approximately \code{n.sims} rows and (No. coefficients + 2) columns. \code{sims.list} is further coerced to be class \code{"mcmc.list"} so that various methods from the package \code{coda} can be directly applied to \code{sims.list} to get Markov Chain diagnostics, posterior summary and plots. See \code{coda} for available methods.    
}

\value{
  \code{bcpglm} returns an object of class \code{"bcpglm"}. See \code{\link{bcpglm-class}} for details of the return values as well as various methods available for this class.
}



\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{bcpglm-class}}, \code{\link{cpglm}}, \code{\link[coda]{mcmc}}, and \code{\link[statmod]{tweedie}} for related information.
}

\examples{

\dontrun{

# fit the fineroot data with Bayesian models
# first use the latent variable approach
set.seed(10)
fit1 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000, method="latent")
gelman.diag(fit1$sims.list)
# diagnostic plots                             
acfplot(fit1$sims.list,lag.max=20)
xyplot(fit1$sims.list)                              
densityplot(fit1$sims.list)               

# summary results
summary(fit1)                             

# now try direct density evaluation (This is much slower)
set.seed(10)
fit2 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000)
gelman.diag(fit2$sims.list)
summary(fit2)                             


# now fit the Bayesian model to an insurance loss triangle 
# (see Peters et al. 2009)
fit3 <- bcpglm(increLoss~ factor(year)+factor(lag), data=insLoss,
               tune.iter=5000, n.tune = 10,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000, bound.p=c(1.1,1.95))
gelman.diag(fit3$sims.list)
summary(fit3)                             

}

}

\keyword{models}

