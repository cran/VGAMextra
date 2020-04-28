\name{VARff}
\alias{VARff}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ VGLTSM family function for the
 Order--\eqn{p} Vector Auto(R)egressive Model}
\description{
  Estimates an Order(\eqn{p}) Vector Autoregressive Models (VAR(p)) with
  white noise random errors
  by maximum likelihood estimation using Fisher scoring.

}
\usage{
            VARff(VAR.order = 1,
                  zero = c("var", "cov"),
                  lmean = "identitylink",
                  lvar  = "loglink",
                  lcov  = "identitylink")
                    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VAR.order}{
   Length--1 (positive) integer vector. The order of the VAR to be 
   fitted.
  
  }
  \item{zero}{
  Integer or character - string vector.
  Same as \code{\link[VGAMextra:MVNcov]{MVNcov}}.
  Details at \code{\link[VGAM:zero]{zero}}.
  
  }
  \item{lmean, lvar, lcov}{
  Same as \code{\link[VGAMextra:MVNcov]{MVNcov}}.
  
  }

}
\details{

Let \eqn{\boldsymbol{x}_t = (x_{1, t}, \ldots, x_{K, t})^T}{
x[t] = (x[1, t], \ldots, x[K, t])^T} be a time dependent
vector of responses, with index \eqn{t = 1, \ldots, T},
and \eqn{\boldsymbol{\varepsilon}_t = (\varepsilon_{1, t},
\ldots, \varepsilon_{K, t})}{e[t] = (e[1, t], \ldots, e[K, t])^T}
white noise with covariance matrix 
\eqn{\boldsymbol{\textrm{V}}}{V}.

\code{VARff} fits a linear model to the means of a 
\eqn{K}--variate normal distribution, where
each variable, \eqn{x_{i, t}}{x[i, t]}, \eqn{i = 1, \ldots, K},
is a linear function of \eqn{p}--past
lags of itself and past \eqn{p}--lags of the other variables.
The model has the form

\deqn{\boldsymbol{x}_t = \boldsymbol{\Phi_1} \boldsymbol{x}_{t - 1} +
\cdots + \boldsymbol{\Phi_p} \boldsymbol{x}_{t - p} +
\boldsymbol{\varepsilon}_t,}{
x[t] = \Phi[1] * x[t - 1] + \ldots + \Phi[p] * x[t - p] + e[t]
}
where \eqn{\boldsymbol{\Phi_j}}{\Phi[j]} are 
\eqn{K \times K}{K x K} matrices of coefficients, \eqn{j = 1, \ldots, K},
to be estimated.

The elements of the covariance matrix are intercept--only
by default.

}
\value{

 An object of class \code{"vglmff"} 
(see \code{\link[VGAM:vglmff-class]{vglmff-class}}) to be
used by VGLM/VGAM modelling functions, e.g.,
\code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.
  

}


\author{

 Victor Miranda.
 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{

  \code{\link[VGAMextra:MVNcov]{MVNcov}},
  \code{\link[VGAM:zero]{zero}},
  \code{\link[VGAM:Links]{Links}},
  \code{\link[VGAMextra:ECM.EngleGran]{ECM.EngleGran}},
  \code{\link[VGAM:vglm]{vglm}}.
  
  
}

\examples{
set.seed(20170227)
nn <- 60
var.data <- data.frame(x2 = runif(nn, -2.5, 2.5))
var.data <- transform(var.data, y1 = rnorm(nn, 1.5 - 2 * x2, sqrt(exp(1.5))),
                                y2 = rnorm(nn, 1.0 - 1 * x2, sqrt(exp(0.75))),
                                y3 = rnorm(nn, 0.5 + 1 * x2, sqrt(exp(1.0))))

fit.var <- vglm(cbind(y1, y2, y3) ~ x2, VARff(VAR.order = 2),
                trace = TRUE, data = var.data)
coef(fit.var, matrix = TRUE)

summary(fit.var)
vcov(fit.var)
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
