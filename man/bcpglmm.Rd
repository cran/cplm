\name{bcpglmm}
\alias{bcpglmm}

\title{
Bayesian Compound Poisson Generalized Linear Mixed Models
}
\description{
This function fits Bayesian compound Poisson generalized linear mixed models using MCMC methods.
}
\usage{
bcpglmm(formula, link = "log", data, inits = NULL, weights, 
    offset, subset, na.action, contrasts = NULL, n.chains = 3, 
    n.iter = 2000, n.burnin = floor(n.iter/2), 
    n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims)), 
    n.sims = 1000, n.report = 1000, prior.beta.mean = NULL, 
    prior.beta.var = NULL, bound.phi = 100, bound.p = c(1.01, 1.99), 
    prior.Sigma = NULL,  tune.iter=4000, n.tune=10, tune.weight=0.25, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a two-sided linear formula object describing the fixed-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. The vertical bar character "|" separates an expression for a model matrix and a grouping factor. See \code{\link[lme4]{glmer}}.
}
  \item{link}{a specification for the model link function. This can be either a literal character string or a numeric number. If it is a character string, it must be one of "log", "identity", "sqrt" or "inverse". If it is numeric, it is the same as the \code{link.power} argument in the \code{\link[statmod]{tweedie}} function. The default is \code{link="log"}.
}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model.
}
 
  \item{subset, weights, na.action, offset, contrasts}{further model specification arguments as in \code{\link[stats]{glm}}; see there for details.}
 
 \item{inits}{a list of initial values to be used for each chain. It must be of length \code{n.chains}. For each element, it is a named list with five components 'beta' (fixed effects), 'phi' (dispersion), 'p' (index parameter), 'u' (random effects) and 'Sigma' (variance components). 'Sigma' must be a list of the same format as the 'ST' slot in \code{cpglmm}. If not supplied, the function will generate initial values automatically from a call of \code{cpglmm} with the supplied model structure. 
}

  \item{n.chains, n.iter, n.burnin, n.thin, n.sims, n.report}{parameters that control the number of chains, iterations, burnins, thinning, simulations to keep and report intervals. See \code{\link{bcpglm}} for details. 
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


\item{tune.iter}{number of iterations used for tuning the variance for the Normal proposal distribution. These iterations will not be used in 
the final output. Default is \code{4000}, and set it to be zero if the tuning process is not desired.  
}
\item{n.tune}{if positive, the \code{tune.iter} iterations will be divided into \code{n.tune} loops. Default is \code{10}.
}
\item{tune.weight}{the weight to be given to the old covariance matrix as opposed to the covariance matrix estimated from simulated samples. Default is \code{0.25}. }

  \item{prior.Sigma}{ a list that specifies the prior for each of the variance components. 
}
  \item{\dots}{
not used }
}
\details{
This function implements the MCMC algorithm for the Bayesian compound Poisson generalized linear mixed model. It is similar to the \code{\link{bcpglm}} function, but only the direct density evaluation method is available currently. It employs the block-wise Metropolis-Hasting for the random effects (\eqn{u}) and direct simulation for the variance components (\eqn{\Sigma}) due to conjugacy. The defualt prior of the variance component for one group is either inverse-Gamma (\code{igamma(scale=0.001, shape=0.001)}) if there is one variance component parameter, or inverse-Wishart (\code{iwish(df=3,scale=diag(1,df))}) if there are more than one parameters. These can be changed in the \code{prior.Sigma} argument.  
}

\value{
\code{bcpglmm} returns an object of class \code{"bcpglmm"}. See \code{\link{bcpglmm-class}} for details of the return values as well as various methods available for this class.
}
\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{bcpglmm-class}}, \code{\link{bcpglm}}, \code{\link{cpglmm}}, \code{\link[coda]{mcmc}}, and \code{\link[statmod]{tweedie}} for related information.
}

\examples{
\dontrun{
fit1 <- bcpglmm(RLD~Zone*Stock+(1|Plant), data = fineroot,  
                   n.chains=3, n.iter=25000, n.burnin=5000,
                   n.sims=2000, n.report=5000)
gelman.diag(fit1$sims.list)                   
summary(fit1)

# now change the prior distribution for the random component
fit2 <- bcpglmm(RLD~Zone*Stock+(1|Plant), data = fineroot,  
                   n.chains=3, n.iter=25000, n.burnin=5000,
                   n.sims=2000, n.report=5000, tune.iter=6000, n.tune=20,                  
                   prior.Sigma = list(igamma(scale=0.01,shape=0.01)))
summary(fit2)                   
}
}
\keyword{models}