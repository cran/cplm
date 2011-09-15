\name{bcpglm}
\alias{bcpglm}

\title{
Bayesian Compound Poisson Generalized Linear Model
}
\description{
This function fits generalized linear models with the Tweedie compound Poisson distribution using Markov Chain Monte Carlo methods.
}
\usage{
bcpglm(formula, link = "log", data, inits = NULL, 
  weights, offset, subset, na.action, contrasts = NULL, 
  n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2), 
  n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), 
  n.sims = 1000, n.report = 1000, prior.beta.mean, prior.beta.var, 
  phi.shape = 0.001, phi.scale = 0.001, bound.p = c(1.01, 1.99), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{formula}{an object of class \code{formula}. See also in \code{\link[stats]{glm}}.
}
  \item{link}{a specification for the model link function. This can be either a literal character string or a numeric number. If it is a character string, it must be one of "log", "identity", "sqrt" or "inverse". If it is numeric, it is the same as the \code{link.power} argument in the \code{\link[statmod]{tweedie}} function. The default is \code{link="log"}.
}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. 
}
 \item{inits}{a list of initial values to be used for each chain. It must be of length \code{n.chains}. For each element, the last two values are for the dispersion and the index parameters, respectively. If not supplied, the function will generate initial values automatically, which are based on a GLM with the supplied model structure. 
}
  \item{weights}{an optional vector of weights. Should be \code{NULL} or a numeric vector. When it is numeric, it must be positive. Zero weights are not allowed in \code{cpglm}. 
}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of options, and is \code{na.fail} if that is unset.  Another possible value is \code{NULL}, no action. Value \code{na.exclude} can be useful.
}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}.
}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. 
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
  \item{phi.shape}{shape parameter in the Inverse-Gamma distribution, the prior of the dispersion parameter. Default is \code{0.001}.
}
  \item{phi.scale}{scale parameter in the Inverse-Gamma distribution, the prior of the dispersion parameter. Default is \code{0.001}.
}
  \item{bound.p}{a vector of lower and upper bound for the index parameter \eqn{p}. The default is \code{c(1.01,1.99)}.
}
  \item{\dots}{ not used.
}
}
\details{
The function fits the Tweedie compound Poisson generalized linear model within the Bayesian framework, where posterior distributions of model parameters are estimated using Monte Carlo Markov Chains [MCMC]. Similar to the Monte Carlo EM [MCEM] algorithm in \code{\link{cpglm}}, posterior simulations also rely on the use of latent Poisson variables to avoid evaluating the Tweedie density and to improve the mixing of the chains. 

Since this is a Bayesian model, prior distributions have to be specified for all parameters in the model. Here, Normal distributions are used for the model coefficients (\eqn{\beta}), Invese-Gamma distribution for the  dispersion parameter (\eqn{\phi}) and Uniform distribution for the index parameter (\eqn{p}). Prior means and variances for the model coefficients can be supplied using the argument \code{prior.beta.mean} and \code{prior.beta.var}, respectively. The prior shape and scale parameter for the dispersion parameter can be supplied using the  argument \code{phi.shape} and \code{phi.scale}, respectively. And the bounds of the Uniform for \eqn{p} can be specified in the argument \code{bound.p}. See details in Section 'Arguments'. 

In implementing the MCMC, a Gibbs sampler is constructed by modifying the MCEM algorithm in \code{\link{cpglm}} in the following ways:
\enumerate{
  \item{The posterior distribution is derived by incorporating the prior distributions of the parameters into the joint loglikelihood in MCEM;}
  \item{In the E-step, draw only one simulation of the latent Poisson variable in each iteration; }
  \item{Replace the M-step with posterior simulations of the parameters in each iteration.}
  }

In simulating the latent variables, rejection sampling with zero-truncated Poisson proposals are still used as in MCEM. Model coefficients are updated block-wise using the Metropolis-Hastings algorithm, where a random-walk multivariate Normal proposal distribution is used. The variance-covariance matrix used in the multivariate Normal is the variance-covariance matrix estimated from a preliminary GLM fit with the specified model structure. Update of the index parameter relies on the Adaptive Rejection Sampling (ARS), where a Metropolis step is also used if the resulting posterior is not log-concave. The posterior distribution of the dispersion parameter is also Inverse-Gamma, and is sampled directly. 

The simulated values of all model parameters are stored in the slot \code{sims.list} of the returned \code{bcpglm} object, which is a list of \code{n.chains} matrices. Each matrix has approximately \code{n.sims} rows and (# coefficients + 2) columns. \code{sims.list} is further coerced to be class \code{"mcmc.list"} so that various methods from the package \code{coda} can be directly applied to \code{sims.list} to get Markov Chain diagnostics, posterior summary and plots. See \code{coda} for details.    
}

\value{
  \code{bcpglm} returns an object of class \code{bcpglm}. See \code{\link{bcplm-class}} for details of the return values as well as various methods available for this class.
}

\note{
The Adaptive Rejection Metropolis Sampling is implemented using Wally Gilks' original C code, available at 
\url{http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html}
}

\references{
\cite{Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995) Adaptive rejection Metropolis sampling. \emph{Applied Statistics}, 44, 455-472.
}

\cite{Zhang Y. (2011). A Latent Variable Approach for Statistical Inference  in Tweedie Compound Poisson Linear Models, \emph{preprint}, \url{http://actuaryzhang.com/publication/bayesianTweedie.pdf}}.
}

\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{bcpglm-class}}, \code{\link{cpglm}}, \code{\link[coda]{mcmc}}, and \code{\link[statmod]{tweedie}} for related information.
}

\examples{

\dontrun{
# default fit, 3 chains, no initial values 
set.seed(10)
fit1 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot)
gelman.diag(fit1$sims.list)
# looks good but has serious autocorrelations                             
acfplot(fit1$sims.list,lag.max=20)                             
                             
# supply initial values
# run the profiled likelihood approach to
# get rough estimates                             
M <- cpglm(RLD~ factor(Zone)*factor(Stock),
  data=fineroot,method="profile", 
  control=list(decimal=1,trace=FALSE))
n.chains <- 3
# inits must be a list of length n.chains           
init.vals <- vector("list",n.chains)
# each element is a vector (beta, phi, p)           
init.vals <- lapply(init.vals, function(x) 
                      c(coef(M)+rnorm(length(coef(M))),
                      M$phi+runif(1), 
                      runif(1,(M$p+1)/2,(M$p+2)/2)))
fit2 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               inits=init.vals, n.chains=n.chains)
               
# run longer chains to get better mixing
fit3 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.chains=3, n.iter=11000, n.burnin=1000,
               n.thin=10)
# since the slot sims.list is a "mcmc.list" object
# we can apply various methods  or functions defined from 
# the package "coda"
# use methods(class="mcmc.list") to see available methods
# and see the documentation of coda for more functions

sims <- fit3$sims.list
gelman.diag(sims)
xyplot(sims)               
acfplot(sims,lag.max=20)                  
densityplot(sims)               
                 
# I also defined plot and summary directly for a "bcpglm" object                 
summary(fit3)
plot(fit3)
}

}

\keyword{models}

