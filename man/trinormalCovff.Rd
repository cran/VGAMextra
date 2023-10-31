\name{trinormalCovff}
\alias{trinormalCovff}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trivariate Normal Distribution Family Function}
\description{
  Estimates the means and the upper-half of the (symmetric) covariance matrix of 
  a trivariate normal distribution by maximum likelihood.

}
\usage{
      trinormalCovff(zero = c("var", "cov"),
                     lmean = "identitylink",
                     lvar  = "loglink",
                     lcov  = "identitylink")
                    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zero}{
    The linear predictors modelled as intercept--only.
   See \code{\link[VGAM:zero]{zero}} for further details.
  
  }
  \item{lmean, lvar, lcov}{
  Link functions applied to the means, variances (diagonal elements of the 
  covariance matrix), and covariances (off-diagonal elements).
  %Same as \code{\link[VGAMextra:MVNcov]{MVNcov}}.
  See \code{\link[VGAM:Links]{Links}} for more choices.
  }

}
\details{

This family function is similar to
\code{\link[VGAM:trinormal]{trinormal}}. The only difference is
that the variances and covariances, instead of the 
standard deviations and correlation coefficients,
are directly modelled and estimated.
Similarly, \code{\link[VGAMextra:trinormalCovff]{trinormalCovff}}
also fits linear models to the means of a trivariate normal distribution.


The fitted means are returned as the fitted values in the form of a
three--column matrix. By default, the variances and covariances 
are modelled as intercept--only, where a
\code{\link[VGAM:loglink]{loglink}} link is applied to the variances
and an \code{\link[VGAM:identitylink]{identitylink}} over
the covariances.

}
\value{

An object of class \code{"vglmff"} 
(see \code{\link[VGAM:vglmff-class]{vglmff-class}}) to be
used by VGLM/VGAM modelling functions, e.g.,
\code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.
  

}


\author{

 Victor Miranda and Thomas Yee.
 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{

  \code{\link[VGAM:trinormal]{trinormal}},
  \code{\link[VGAM:zero]{zero}},
  \code{\link[VGAM:Links]{Links}},
  \code{\link[VGAM:vglm]{vglm}}.
  
  
}

\examples{
set.seed(123); nn <- 350
var1 <- exp(1.5); var2 <- exp(0.75); var3 <- exp(1.0)

### Artificial data, with two covariates.
tdata <- data.frame(x2 = runif(nn), x3 = runif(nn))
tdata <- transform(tdata,
                   y1 = rnorm(nn, 1 + 2 * x2, sd = sqrt(var1)),
                   y2 = rnorm(nn, 3 + 1 * x2, sd = sqrt(var2)),
                   y3 = rnorm(nn, 3 - 1 * x3, sd = sqrt(var2 * var3)))

### Fit the model using VGAMextra::trinormalCovff().
fit.trinormCovff <- vglm(cbind(y1, y2, y3) ~ x2 + x3,
                         trinormalCovff,
                         data = tdata, trace = TRUE)

summary(fit.trinormCovff)
vcov(fit.trinormCovff)
                         
### Fitting the model using VGAM::trinormal()
fit.trinormVGAM <- vglm(cbind(y1, y2, y3) ~ x2 + x3,
                        trinormal,
                        data = tdata, trace = TRUE)
                         
summary(fit.trinormVGAM)
vcov(fit.trinormVGAM)

                         
#### Compare the estimated coefficients. Note that 
#### trinormal() estimates the sd's and correlation coeffs.
coef(fit.trinormCovff, matrix = TRUE)
coef(fit.trinormVGAM, matrix = TRUE)


}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
