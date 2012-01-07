\name{cpglmm}
\alias{cpglmm}
\title{
Compound Poisson Generalized Linear Mixed Models
}
\description{
Laplace approximation and adaptive Gauss-Hermite quadrature methods for compound Poisson generalized linear mixed models and additive models. 
}
\usage{
cpglmm(formula, link = "log", data, weights, offset, subset, 
    na.action, inits = NULL,  contrasts = NULL, 
    control = list(), basisGenerators = c("tp","tpU","bsp","sp2d"),
    optimizer = "nlminb", doFit = TRUE, nAGQ = 1)
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
  
\item{inits}{a named list with three components 'beta', 'phi', 'p', 'Sigma' that supply the initial values used in the optimization. If not supplied, the function will generate initial values automatically, which are based on a GLM with the supplied model structure.
}
  \item{control}{
a list of parameters for controlling the fitting process. See 'Details' below. 
}
  \item{basisGenerators}{
a character vector of names of functions that generate spline bases. This is used when smoothing effects are to be included in the model. See \code{\link[amer]{amer}} and \code{\link[amer]{tp}} for details. 
}
  \item{optimizer}{
a character string that determines which optimization routine is to be used. Possible choices are \code{"nlminb"} (the default, see \code{\link[stats]{nlminb}}), \code{"bobyqa"} (\code{\link[minqa]{bobyqa}}) and \code{"L-BFGS-B"} (\code{\link[stats]{optim}}).   
}
  \item{doFit}{
if \code{FALSE}, no model will be fitted and the constructed \code{"cpglmm"} object is returned. 
}
  \item{nAGQ}{
  a positive integer - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. This defaults to 1, corresponding to the Laplacian approximation. Values greater than 1 produce greater accuracy in the evaluation of the log-likelihood at the expense of speed.
  }
}

\details{
The estimation of the mixed-effect models in the compound Poisson distribution in  existing software has been limited to the PQL approach (Penalized Quasi-Likelihood, see \code{\link[MASS]{glmmPQL}}). While straightforward and fast, this method, because of its quasi-likelihood nature, does not offer a neat solution to the problem of estimating the index parameter. In contrast, the function \code{cpglmm} implements true likelihood-based inferential procedures, i.e., the Laplace approximation and the Adaptive Gauss-Hermite Quadrature (for single grouping factor), so that all parameters in the model can be estimated using maximum likelihood estimations. These methods are often considered to have better asymptotic properties.  
 
This implementation is based on the code for  \code{\link[lme4]{glmer}} in the \code{lme4} package, with changes made on the updating of the mean, the variance function and the marginal loglikelihood. For the Laplace methods, one big difference to \code{glmer} is that the contribution of the dispersion parameter to the approximated loglikelihood is explicitly accounted for, which I think is more accurate and more consistent with the quadrature estimate. As a result, both the dispersion parameter and  the index parameter are included as a part of the optimization process. In computing the marginal loglikelihood, the density of the compound Poisson distribution is approximated using numerical methods provided in the \code{tweedie} package. For details of the Laplace approximation and the Gauss-Hermite quadrature method for generalized linear mixed models, see the documentations associated with \code{lme4}. 

In addition, similar to the package \code{amer}, we provide convenient interfaces for fitting additive models with penalized splines using the mixed-effect model facility.  All the spline constructors available in \code{amer} or those defined by the users following the \code{amer} requirement can be directly used here. See the 'example' section for one such application.  

The \code{control} argument is a list that can supply various controlling elements used in the optimization process:

\describe{
\item{\code{max.iter}}{maximum number of iterations allowed in the optimizer. The default value is 300.}
\item{\code{max.fun}}{maximum number of function evaluations allowed in the optimizer. }
\item{\code{bound.p}}{a vector of lower and upper bound for the index parameter \eqn{p}. The default is \code{c(1.01,1.99)}. }
\item{\code{trace}}{if greater than 0, tracing information on the progress of the fitting is produced. For \code{optimizer == "nlminb"} or \code{optimizer == "L-BFGS-B"} , this is the same as the \code{trace} control parameter, and for \code{optimizer == "bobyqa"}, this is the same as the \code{iprint} control parameter. See the corresponding documentation for details. 
}

}

}
\value{
  \code{cpglmm} returns an object of class \code{cpglmm}. See \code{\link{cpglmm-class}} for details of the return values as well as various method available for this class.
}
\author{
Wayne (Yanwei) Zhang \email{actuary_zhang@hotmail.com}
}

\seealso{
The users are recommended to see \code{\link{cpglm}} for a general introduction to the compound Poisson distribution, \code{\link[lme4]{glmer}} for syntax and usage of mixed-effect models and \code{\link{cpglmm-class}} for detailed explanation of the return value.
}
\examples{

# use Stock and Spacing as main effects and Plant as random effect
(f1 <- cpglmm(RLD ~ Stock + Spacing +  (1|Plant), data = fineroot))
            
# most of the methods defined in lme4 are directly applicable
coef(f1); fixef(f1); ranef(f1)  #coefficients
VarCorr(f1)  #variance components

# add another random effect
(f2 <- update(f1, . ~ . + (1|Zone)))
# test the additional random effect
anova(f1,f2)

# try a different optimizer 
(f3 <- cpglmm(RLD ~  Stock + Spacing +  (1|Plant), 
            data = fineroot, optimizer = "bobyqa", 
            control = list(trace = 2)))

# adaptive G-H quadrature  
(f4 <- cpglmm(RLD ~  Stock + Spacing +  (1|Plant), 
            data = fineroot, nAGQ = 3))

# a model with smoothing effects
(f5 <- cpglmm(increLoss ~ tp(lag, k = 4) + (1|year) , 
            data = insLoss))
            
}

\keyword{ models}
