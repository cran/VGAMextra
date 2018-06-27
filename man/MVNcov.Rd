\name{MVNcov}
\alias{MVNcov}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multivariate Normal Distribution Family Function}
\description{

  Maximum likelihood estimation of the Multivariate Normal
  distribution. The covariances (not correlation coefficients)
  are included in the parameter vector.

}
\usage{
      MVNcov(zero = c("var", "cov"),
             lmean = "identitylink",
             lvar  = "loge",
             lcov  = "identitylink")
 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zero}{
   Integer or character--string vector.
   Which linear predictors are intercept--only.
   Details at \code{\link[VGAM:zero]{zero}} or
   \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
   
  
  }
  \item{lmean, lvar, lcov}{
  VGLM--link functions applied to the means, variances and covariances.
  
   }
  }
  \details{
  
  For the \eqn{K}--dimensional normal distribution, 
  this fits a linear model to the \eqn{K} means 
  \eqn{\mu_j}{\mu[j]} \eqn{j = 1, \ldots, K}, which are the 
  first entries in the linear predictor. The variances
  \eqn{\sigma^2_j}{\sigma^2[j]} \eqn{j = 1, \ldots, K} 
  and then the covariances
  \eqn{cov_{ij}}{cov[ij]} \eqn{i = 1, \ldots, K,
  j = i + 1, \ldots, K}, are next aligned.
  The fitted means are returned as the fitted values.
  
  The log--likelihood is computed via
  \code{\link[VGAMextra:dmultinorm]{dmultinorm}},
  an implementation of the multivariate Normal density.
  
  The score and 
  expected information matrices are internally computed
  at each Fisher scoring step, using its vectorized form.
  
  The response should be an \eqn{K}--column matrix.
  The covariances may be any real number so that
  the \code{\link[VGAM:identitylink]{identitylink}}
  is a reasonable choice.
  For further details on VGLM/VGAM--link functions, see
  \code{\link[VGAM:Links]{Links}}.

}
\value{

 An object of class \code{"vglmff"} 
(see \code{\link[VGAM:vglmff-class]{vglmff-class}}) to be
used by VGLM/VGAM modelling functions, e.g.,
\code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.
  

}

\note{


Unlike other implementations, e.g., 
\code{\link[VGAM:binormal]{binormal}} from 
\pkg{VGAM} in terms of \eqn{\rho} and standard deviations,
\code{MVNcov} estimates
the variances and covariances, modeled as intercept--only.
See argument \code{zero}, whose default is \code{c("var", "cov")},
to change this.

% To get rid of this, incorporate the check at @last or implement @validparams
Thus far, there is no guarantee that the estimated var--cov matrix will be 
positive--definite. Proper procedures to validate this will
be incorporated shortly, such as the \code{@validparams} slot.

Although the function has been tested on \eqn{K \leq 5}{K \le 5}
data sets, it is recommended that \eqn{K \leq 3}{K \le 3},
unless the data are \emph{nice} and \eqn{n} is sufficiently large.


}

\author{

 Victor Miranda.
 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{

  \code{\link[VGAMextra:dmultinorm]{dmultinorm}},
  \code{\link[VGAM:binormal]{binormal}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[VGAM:vglm]{vglm}}.
  
}

\examples{
# K = 3.
set.seed(180227)
nn  <- 85
mvndata <- data.frame(x2 = runif(nn), x3 = runif(nn))
mvndata <- transform(mvndata, 
                     y = rbinorm(nn, mean1 = 2 - 2 * x2 + 1 * x3,
                          mean2 = 2 - 1.5 * x3,
                          var1 = exp(1.0), var2 = exp(-0.75),
                          cov12 = 0.5 * exp(1.0) * exp(-0.75)))
mvndata <- transform(mvndata, y3 = rnorm(nn, mean = 2 + x2, sd = exp(1.5)))
colnames(mvndata) <- c("x2", "x3", "y1", "y2", "y3")

mvnfit <- vglm(cbind(y1, y2, y3) ~ x2 + x3, MVNcov, data = mvndata, trace = TRUE)
(mvncoef  <- coef(mvnfit, mat = TRUE))

## Check variances and covariances: var1, var2 and var3.
exp(mvncoef[c(10, 13, 16)])  # True are var1 = exp(1.0) = 2.718, 
                             # var2 = exp(-0.75) = 0.472
                             # and var3 = exp(1.5)^2 = 20.08554
vcov(mvnfit)
constraints(mvnfit)
summary(mvnfit)


}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
