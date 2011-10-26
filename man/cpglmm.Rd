\name{cpglmm}
\alias{cpglmm}
\title{
Compound Poisson Generalized Linear Mixed Models
}
\description{
This function fits compound Poisson generalized linear mixed models using Laplacian approximation.
}
\usage{
cpglmm(formula, link = "log", data, weights, offset, subset, 
    na.action, betastart = NULL, phistart = NULL, pstart = NULL, 
    contrasts = NULL, control = list())
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
  
\item{betastart}{starting values for the vector of mean parameters.
}
  \item{phistart}{starting value for the dispersion (scale) parameter.
}
  \item{pstart}{starting value for the index parameter. Must be between 1 and 2. When \code{bound.p} is specified, \code{pstart} must be between the bounds set there.
}  
  \item{control}{
a list of parameters for controlling the fitting process. See 'Details' below. 
}
}

\details{
The Tweedie compound Poisson distribution is composed of a probability mass at the origin accompanied by a skewed continuous distribution on the positive values. It has been widely used to accommodate zero-inflated continuous data, but primarily in the form of generalized linear models (see \code{\link{cpglm}} for more details). The function \code{cpglmm} now enables users to apply mixed-effect models when working with the Tweedie compound Poisson distribution, thus enriching the set of tools to handle zero-inflated data.
 
Specifically, \code{cpglmm} fits compound Poisson generalized linear mixed models [CPGLMM] using Laplacian approximation. The implementation is based on the code for  \code{\link[lme4]{glmer}} in the \code{lme4} package, with changes made on the updating of the mean, the variance function and the marginal loglikelihood. In computing the marginal logliklihood, the density of the compound Poisson distribution is approximated using numerical methods provided in the \code{tweedie} package, and the parameters in the optimization include the dispersion parameter and the index parameter in additional to the mean parameters and the variance components. For details of the implementation, see the documentations in \code{lme4}. 



The \code{control} argument is a list that can supply various controlling elements used in the fitting process. Most of these elements are used in the MCEM algorithm, and those used in the profile likelihood approach are explicitly identified in the following:

\describe{
\item{\code{max.iter}}{maximum number of iterations allowed in the algorithm. The default value is 300.}
\item{\code{max.fun}}{maximum number of function evaluations allowed in the algorithm. The default value is 900.}
\item{\code{bound.p}}{a vector of lower and upper bound for the index parameter \eqn{p}. The default is \code{c(1.01,1.99)}. }
\item{\code{trace}}{if \code{TRUE}, tracing information on the progress of the fitting is produced. }

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
f1 <- cpglmm(RLD~  Stock + Spacing +  (1|Plant) , 
            link="log", data = fineroot)
summary(f1)
# most of the methods defined in lme4 are directly applicable
coef(f1); fixef(f1); ranef(f1)  #coefficients
VarCorr(f1)  #variance components

# add another random effect
f2 <- cpglmm(RLD~  Stock + Spacing +  (1|Plant) + (1|Zone), 
            link="log", data = fineroot)
# test the additional random effect
anova(f1,f2)

}

\keyword{ models}
