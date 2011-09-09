\name{cpglm}
\alias{cpglm}
\alias{cpglm_em}
\alias{cpglm_profile}
\title{Compound Poisson Generalized Linear Model
}
\description{This function fits generalized linear models with the Tweedie compound Poisson distribution. 
}
\usage{
cpglm(formula, link = "log", data, weights, offset, subset, 
  na.action, betastart = NULL,phistart = NULL, pstart = NULL, 
  contrasts = NULL, control = list(), method="MCEM", ...)
  
cpglm_em(X, Y, weights = NULL, offset = NULL, link.power = 0,
  betastart, phistart, pstart, intercept = TRUE, control = list())  
  
cpglm_profile(X, Y, weights = NULL, offset = NULL, 
  link.power = 0, intercept = TRUE, contrasts, control = list())                      
}

\arguments{
  \item{formula}{an object of class \code{formula}. See also in \code{\link[stats]{glm}}.
}
  \item{link}{a specification for the model link function. This can be either a literal character string or a numeric number. If it is a character string, it must be one of "log", "identity", "sqrt" or "inverse". If it is numeric, it is the same as the \code{link.power} argument in the \code{\link[statmod]{tweedie}} function. The default is \code{link="log"}.
}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. 
}
  \item{weights}{an optional vector of weights. Should be \code{NULL} or a numeric vector. When it is numeric, it must be positive. Zero weights are not allowed in \code{cpglm}. 
}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of options, and is \code{na.fail} if that is unset.  Another possible value is \code{NULL}, no action. Value \code{na.exclude} can be useful.
}
  \item{betastart}{starting values for the vector of mean parameters.
}
  \item{phistart}{starting value for the dispersion parameter.
}
  \item{pstart}{starting value for the index parameter. Must be between 1 and 2. When \code{bound.p} is specified, \code{pstart} must be between the bounds set there.  
}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. 
}
  \item{contrasts}{an optional list. See the \code{contrasts.arg}.
}
  \item{control}{a list of parameters for controlling the fitting process. See 'Details' below. 
} 
  \item{method}{a character string, either \code{"MCEM"} or \code{"profile"}, indicating whether the Monte Carlo EM algorithm or the profiled log-likelihood method should be used.  See 'Details' below.
}
  \item{X, Y}{Used in \code{cpglm_em} and \code{cpglm_profile}: \code{X} is a design matrix and \code{Y} is a vector of observations.  
}
  \item{intercept}{logical. Should an intercept be included in the \code{null} model?
}
 \item{link.power}{see  \code{\link[statmod]{tweedie}}.
}
  \item{\dots}{ not used now. 
}

}

\details{
The Tweedie compound Poisson distribution is a subclass of the exponential dispersion family with a power variance function, where the value of the power index lies in the interval (1,2). The distribution is composed of a probability mass at the origin accompanied by a skewed continuous distribution on the positive values. Despite its ability to handle continuous zero-inflated data, the density of the compound Poisson distribution is not analytically tractable. As a result, full maximum likelihood inference remains challenging when the power index is unknown. The \code{cpglm} function implements two methods to fit generalized linear models with the compound Poisson distribution, namely, the  Monte Carlo EM [MCEM]  algorithm and the profile likelihood approach, both of which yield maximum likelihood estimates for all parameters in the model, in particular, for the dispersion and the index parameter. 
 
The MCEM approach is in light of the fact that the compound Poisson distribution involves an unobserved Poisson process. The observed data is thus augmented by the latent Poisson variable so that one could work on the joint likelihood of the complete data and invoke the MCEM algorithm to locate the maximum likelihood estimations. Specifically, a Monte Carlo method is implemented in the E-step, where samples of the unobserved Poisson process is simulated from its posterior via rejection sampling (using a zero-truncated Poisson as a proposal), and the required expectation in the E-step is approximated by a Monte Carlo estimate. The M-step involves conditional component-wise maximization: scoring algorithm is used for the mean parameters, constrained optimization for the index parameter and close-form solution exists for the scale parameter. The E-step and the M-step are iterated until the following stopping rule is reached:
\deqn{\max_i \frac{|\theta_i^{(r+1)}-\theta_i^{(r)}|}{\theta_i^{(r)}+\epsilon_1}<\epsilon_2,
}
where \eqn{\theta_i} is the \eqn{i_{th}} element of the vector of parameters, \eqn{r} represents the \eqn{r_{th}} iteration, and \eqn{\epsilon_1} and \eqn{\epsilon_2} are predefined constants. See the description for the \code{control} parameters below. 

Alternatively, one could use the profiled likelihood approach, where the index parameter is estimated based on the profiled likelihood first and then a traditional GLM is fitted. This is because the mean parameters and the index parameter (\eqn{p}) are orthogonal, meaning that the mean parameters vary slowly as \eqn{p} changes. Specifically, to get a profiled maximum likelihood estimate of \eqn{p}, one needs to specify a grid of possible values of \eqn{p}, estimates the corresponding mean and dispersion parameters, and computes the associated likelihood. Then these likelihoods are compared to locate the value of \eqn{p} corresponding to the maximum likelihood.  To compute the likelihood, one has to resort to numerical methods provided in the  \code{tweedie} package that approximate the density of the compound Poisson distribution. Indeed, the function  \code{\link[tweedie]{tweedie.profile}} in that package makes available the profiled likelihood approach. What is included here is simply a wrapper that provides automatic estimation of \eqn{p} to any pre-specified decimal places (see \code{decimal} in the following). Only MLE estimate for \eqn{\phi} is included here, while \code{\link[tweedie]{tweedie.profile}} provides several other possibilities.  


These two methods are specified in the \code{method} argument in \code{cpglm}. The MCEM algorithm is invoked by specifying \code{method="MCEM"} in \code{cplm}, while the profile likelihood approach is used when \code{method="profile"}. The two methods can also be invoked by calling the function \code{cpglm_em} and \code{cpglm_profile}, however, the users are encouraged to use \code{cpglm}, which uses \code{cpglm_em} and \code{cpglm_profile} internally.


It is to be noted that the MCEM algorithm enables one to work on the more tractable joint likelihood, thus simplifying the estimation problem down to posterior simulation and block/univariate optimizations. However, it also has some overheads:
\itemize{
  \item The simulation of the latent variable could be computationally demanding when a large number of samples is needed;
  \item The inherent Monte Carlo error must be accounted for;
  \item The EM algorithm could have slow convergence rate.
}
The computational demand can be largely reduced by supplementing the algorithm with importance sampling after it converges to the neighborhood of the MLEs. The importance sampling retains the latest simulated samples for use in all later iterations, thus eliminating the need for simulating new samples. 

To account for the Monte Carlo error, many authors have suggests increasing the sample size as the algorithm proceeds. Booth and Hobert (1999) derives a Normal approximation for the Monte Carlo error, which is further used to construct a (1-\code{alpha}) confidence ellipsoid at each iteration. If the old parameters lie in this ellipsoid, then the sample size is increase as \eqn{m + m/k}, where \eqn{m} is the current sample size, and \eqn{k} is a predefined constant. This method is implemented when \code{fixed.size} in the control argument is \code{FALSE}. See the description for the \code{control} parameters below.    



The \code{control} argument is a list that can supply various controlling elements used in the fitting process. Most of these elements are used in the MCEM algorithm, and those used in the profile likelihood approach are explicitly identified in the following:

\describe{
  \item{\code{init.size}}{initial number of samples used to generate the latent Poisson variable. Default is 100.}
\item{\code{sample.iter}}{the iteration number after which importance sampling is used where the samples in iteration \code{sample.iter}-1 will be recycled in all later iterations. This could improve the speed in the E-step, but should only be used when the algorithm has converged to the neighborhood of the MLE. The default value is 50. To suppress importance sampling, set this number to be large enough.}
\item{\code{max.iter}}{maximum number of iterations allowed in the algorithm. The default value is 200.}
\item{\code{epsilon1}}{a constant used in the stopping rule. The default value is 0.001.}
\item{\code{epsilon2}}{a constant used in the stopping rule. The default value is 0.0001.}
\item{\code{fixed.size}}{a logical value indicating whether the sample size in each iteration is to be increased. If \code{FALSE}, approximate variance is computed in each iteration and the approximate (1-\code{alpha}) confidence ellipsoid is constructed. If the previously simulated value falls in this ellipsoid, then the sample size is increased by 1/\code{k} of the current sample size.  }
\item{\code{alpha}}{the confidence level used in the above calculation of confidence ellipsoid. The default value is 0.25.}
\item{\code{k}}{a constant used in the calculation of new sample size. The default value is 5.}
\item{\code{max.size}}{the maximum sample size used in each iteration. When \code{fixed.size} is \code{FALSE}, the sample size in later iterations could be increasing without bound and as a result, the algorithm could be impractically slow. \code{max.size} sets the limit of the required sample size to avoid such slow computation. The default value is 10000. }
\item{\code{bound.p}}{a vector of lower and upper bound for the index parameter \eqn{p}. The default is \code{c(1.01,1.99)}. Used in both \code{cpglm_em} and \code{cpglm_profile}.}
\item{\code{beta.step}}{an integer. Since the mean parameters are orthogonal to \eqn{p} and not dependent on the latent variable, it is not necessary to update them in each step. By specifying \code{beta.step}, update of the mean parameters will only be performed every \code{beta.step} steps. It defaults to 10.}
\item{\code{trace}}{if \code{TRUE} (the default), tracing information on the progress of the fitting is produced. Used in both \code{cpglm_em} and \code{cpglm_profile}.}
\item{\code{profile.method}}{a character string indicating the method used in the Tweedie density evaluation. The default is \code{"inversion"}. If the default method fails, try other methods such as \code{"series"}. See \code{\link[tweedie]{dtweedie}} and \code{\link[tweedie]{tweedie.profile}}. This is only used the profile likelihood approach.}
\item{\code{decimal}}{an integer indicating the accuracy of the profile likelihood approach in terms of the number of decimal places of \eqn{p}. For example, if \code{decimal=1}, then \eqn{p=1.1, 1.2, \dots, 1.9} (if they are within \code{bound.p}) will be used as the grid of possible values to estimate \eqn{p}. If \code{decimal=2} and suppose the estimated \eqn{p} from \code{decimal=1} is 1.4, then the grid is updated to be \code{seq(1.31,1.49, by=0.01)} and a new estimate of \eqn{p} is obtained. This procedure continues for higher accuracy. This is only used the profiled likelihood approach. }

}

}

\value{
  \code{cpglm} returns an object of class \code{cpglm}. See \code{\link{cpglm-class}} for details of the return values as well as various method available for this class. 
}

\references{
\cite{Booth, J. G., and Hobert, J. P. (1999). Maximizing Generalized Linear Mixed Model Likelihoods with an Automated
Monte Carlo EM Algorithm.  \emph{Journal of the Royal Statistical Society Series B}, 61, 265-285.}

\cite{ Dunn, P.K. and Smyth, G.K. (2005). Series evaluation of Tweedie exponential dispersion models densities. \emph{Statistics and Computing}, 15, 267-280.}

\cite{McCulloch, C. E. (1997). Maximum likelihood algorithms for generalized linear mixed
models. \emph{Journal of the American Statistical Association}, 92, 162-170.}

\cite{Zhang Y. (2011). A Latent Variable Approach for Statistical Inference  in Tweedie Compound Poisson Linear Models, \emph{preprint}, \url{http://actuaryzhang.com/publication/bayesianTweedie.pdf}}.

}

\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see the documentation for \code{\link{cpglm-class}}, \code{\link[stats]{glm}}, \code{\link[statmod]{tweedie}}, and \code{\link[tweedie]{tweedie.profile}} for related information.
}


\examples{

\dontrun{
# MCEM fit
set.seed(10)
fit1 <- cpglm(RLD~ factor(Zone)*factor(Stock),
  data=fineroot,
  control=list(init.size=5,sample.iter=50,
              max.size=2000,fixed.size=FALSE),
  pstart=1.6)

# show the iteration history of p
plot(fit1$theta.all[,ncol(fit1$theta.all)],
     type="o", cex=0.5)
     
# residual and qq plot
parold <- par(mfrow=c(2,2), mar=c(5,5,2,1))
# 1. regular plot
r1 <- resid(fit1)/sqrt(fit1$phi)
plot(r1~fitted(fit1), cex=0.5)
qqnorm(r1, cex=0.5)
# 2. quantile residual plot to avoid overlappling
u <- ptweedie(fit1$y, fit1$p, fitted(fit1), fit1$phi)
u[fit1$y == 0] <- runif(sum(fit1$y == 0), 0, u[fit1$y == 0])
r2 <- qnorm(u)
plot(r2~fitted(fit1), cex=0.5)
qqnorm(r2, cex=0.5)

par(parold)

# profile likelihood         
fit2 <- cpglm(RLD~ factor(Zone)*factor(Stock),
  data=fineroot,method="profile", 
  control=list(decimal=3))      

# compare the two
summary(fit1)
summary(fit2)              
            
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models}

