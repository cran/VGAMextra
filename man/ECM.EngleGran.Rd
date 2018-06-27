\name{ECM.EngleGran}
\alias{ECM.EngleGran}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ VGLTSM family function for the
Two--dimensional Error--Correction Model 
(Engle and Granger, 1987) for \eqn{I(1)}--variables}
\description{


  Estimates a bidimensional error-correction model of
  order--(\eqn{K}, \eqn{L}), as proposed by Engle--Granger
  (Two step--approach; 1987), with bivariate normal errors
  by maximum likelihood estimation using Fisher scoring.
  
  

}
\usage{
      ECM.EngleGran(ecm.order = c(1, 1),
                    zero = c("var", "cov"),
                    resids.pattern = c("intercept", "trend",
                                       "neither", "both")[1],
                    lag.res = 1, 
                    lmean = "identitylink",
                    lvar  = "loge",
                    lcov  = "identitylink",
                    ordtsDyn = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ecm.order}{
   Length--2 (positive) integer vector. The order of the ECM model.
  
  }
  \item{zero}{
  Integer or character--string vector.
  %Same as \code{\link[VGAMextra:MVNcov]{MVNcov}}.
  Details at \code{\link[VGAM:zero]{zero}}.
  
  }
  \item{resids.pattern}{
  Character. How the static linear regression
  \eqn{y_{2, t} \sim y_{1, t}}{y[2, t] ~ y[1, t]} must be 
  settle to estimate the residuals \eqn{\widehat{z_t}}{zhat[t]}.
  The default is a linear model with intercept, and no 
  trend term.
  See below for details. 
  
  
  }
  \item{lag.res}{
  Numeric, single positive integer. The error
  term for the long--run equilibrium path is
  lagged up to order \code{lag.res}. See below
  for further details.
  
  }
  \item{lmean, lvar, lcov}{
  Same as \code{\link[VGAMextra:MVNcov]{MVNcov}}.
  
  }
  \item{ordtsDyn}{
  Positive integer. Allows to compare the estimated coefficients with
  those provided by the package 'tsDyn'. See below for further details.
  
  }

}
\details{

This is an implementation of the two--step approach as
proposed by Engle--Granger [1987] to estimate an
order--(\eqn{K}, \eqn{L}) bidimensional error correction
model (ECM) with bivariate normal errors. 


This ECM class models the dynamic behaviour of two 
cointegrated \eqn{I(1)}-variables,
say \eqn{y_{1, t}}{y[1, t]} and \eqn{y_{2, t}}{y[2, t]} with,
probably, \eqn{y_{2, t}}{y[2, t]} a function of 
\eqn{y_{1, t}}{y[1, t]}.
Note, the response must be a two--column matrix,
 where the first entry is the regressor,
 i.e, \eqn{y_{1, t}}{y[1, t]} above, and the regressand in the
 second colum. See Example 2 below.

The general specification of the ECM class described by 
this family function is

\deqn{ \Delta y_{1, t} |\Phi_{t - 1} =
~\phi_{0, 1} + \gamma_1 \widehat{z}_{t - k} +
 \sum_{i = 1}^K \phi_{1, i} \Delta y_{2, t - i} +
 \sum_{j = 1}^L \phi_{2, j} \Delta y_{1, t - j} +
 \varepsilon_{1, t},}{
 \Deltay[1, t] |\Phi[t - 1] = \phi[0, 1] + \gamma[1] zhat[t - k] + 
 \sum_i^K \phi[1, i] \Delta y[2, t - i] +
 \sum_j^L \phi[2, j] \Delta y[1, t - j] + e[1, t],
 }
 
 
\deqn{ \Delta y_{2, t} |\Phi_{t - 1}=
~\psi_{0, 1} + \gamma_2 \widehat{z}_{t - k} +
 \sum_{i = 1}^K \psi_{1, i} \Delta y_{1, t - i} +
 \sum_{j = 1}^L \psi_{2, j} \Delta y_{2, t - j} +
 \varepsilon_{2, t}.}{
 \Deltay[2, t] |\Phi [t - 1] = \psi[0, 1] + \gamma[2] zhat[t - k] + 
 \sum_i^K \psi[1, i] \Delta y[1, t - i] +
 \sum_j^L \psi[2, j] \Delta y[2, t - j] + e[1, t].
 }

Under the binormality assumption on the errors
\eqn{(\varepsilon_{1, t}, \varepsilon_{2, t})^T}{
(e[1, t], e[2, t])^T
} with covariance matrix \eqn{\boldsymbol{\textrm{V}}}{V},
model above can be seen as a
VGLM fitting linear models over the conditional means,
\eqn{\mu_{\Delta y_{1, t} } = E(\Delta y_{1, t} | \Phi_{t - 1})}{
\mu\Deltay[1, t] = E(\Deltay[1, t] | \Phi[t - 1])} and
\eqn{\mu_{\Delta y_{2, t} } = E(\Delta y_{2, t} | \Phi_{t - 1} )}{
\mu\Deltay[2, t] = E(\Deltay[2, t] | \Phi[t - 1])},
producing

\deqn{ (\Delta y_{1, t} |\Phi_{t - 1} , 
\Delta y_{2, t} |\Phi_{t - 1} )^T \sim N_{2}(\mu_{\Delta y_{1, t}}  , 
\mu_{\Delta y_{2, t}}, \boldsymbol{\textrm{V}})}{
(\Deltay[1, t] | \Phi[t - 1], \Deltay[2, t] | \Phi[t - 1])^T  ~  
N2(\mu \Deltay[1, t],  \mu \Deltay[2, t] , V).
}

The covariance matrix is assumed 
to have elements 
\eqn{\sigma_1^2, \sigma_2^2,}{
\sigma^2[1], \sigma^2[2],} and 
\eqn{\textrm{Cov}_{12}.}{Cov[12].}

Hence, the parameter vector is
\deqn{\boldsymbol{\theta} = 
    (\phi_{0, 1}, \gamma_1, \phi_{1, i}, \phi_{2, j},
                      \psi_{0, 1}, \gamma_2, \psi_{1, i}, \psi_{2, j},
                      \sigma_1^2, \sigma_2^2, \textrm{Cov}_{12})^T,}{
      \theta = (\phi[0, 1], \gamma[1], \phi[1, i], \phi[2, j],
      \psi[0, 1], \gamma[2], \psi[1, i], \psi[2, j],
      \sigma^2[1], \sigma^[2], Cov[12])^T,}
for \eqn{i = 1, \ldots, K} and \eqn{j = 1, \ldots, L}. 


The linear predictor is
\deqn{\boldsymbol{\eta} = (\mu_{\Delta y_{1, t}  },
                          \mu_{\Delta y_{2, t} },
                          {\color{blue}\texttt{loge}}~\sigma_1^2,
                          {\color{blue}\texttt{loge}}~\sigma_2^2,
                          \textrm{Cov}_{12})^T.}{
                          \eta = (\mu \Deltay[1],
                                  \mu \Deltay[2],
                                  log \sigma^2[1], log \sigma^2[2],
                                  Cov[12] )^T.}


The estimated cointegrated vector,
\eqn{\boldsymbol{\widehat{\beta^{\star}}}}{
\betahat*} = \eqn{(1, -\boldsymbol{\widehat{\beta})}^T}{
(1, -\betahat)^T} is obtained by linear regression
depending upon \code{resids.pattern},
as follows:

1) \eqn{y_{2, t} = \beta_0 + \beta_1 y_{1, t} + z_t}{
y[2, t] = \beta[0] + \beta[1] * y[1, t] + z[t]}, 
if \code{resids.pattern = "intercept"},

2) \eqn{y_{2, t} = \beta_1 y_{1, t} + \beta_2 t + z_t}{
y[2, t] = \beta[1] * y[1, t] + \beta[2] * t + z[t]}, 
if \code{resids.pattern = "trend"},

3) \eqn{y_{2, t} = \beta_1 y_{1, t} + z_t}{
y[2, t] = \beta[1] * y[1, t] + z[t]}, 
if \code{resids.pattern = "neither"}, or else,

4) \eqn{y_{2, t} = \beta_0 + \beta_1 y_{1, t} + \beta_2 t + z_t}{
y[2, t] = \beta[0] + \beta[1] * y[1, t] + \beta[2] * t + z[t]}, 
if \code{resids.pattern = "both"},


where \eqn{\boldsymbol{\widehat{\beta}} = 
(\widehat{\beta_0}, \widehat{\beta_1}, \widehat{\beta_2})^T,}{
\betahat* = (\beta[0], \beta[1], \beta[2])^T,
} and \eqn{z_t}{z[t]} assigns the error term.

Note, the \emph{estimated residuals}, 
\eqn{\widehat{z_t}}{zhat[t]} are (internally) computed
from any of the linear models 1) -- 4) selected, and then lagged
up to order \code{alg.res},
and embedded as explanatories in models 
\eqn{\Delta y_{1, t} |\Phi_{t - 1}}{\Deltay[1, t] |\Phi[t - 1]} and
\eqn{\Delta y_{3, t} |\Phi_{t - 1}}{\Deltay[2, t] |\Phi[t - 1]} above.
By default, \eqn{\widehat{z}_{t - 1}}{zhat[t - 1]}
are considered (as \code{lag.res} = 1), although it may be any lag
\eqn{\widehat{z}_{t - k}}{zhat[t - k]}, for \eqn{k > 0}. 
Change this through argument \code{lag.res}.


}
\value{

 An object of class \code{"vglmff"} 
(see \code{\link[VGAM:vglmff-class]{vglmff-class}}) to be
used by VGLM/VGAM modelling functions, e.g.,
\code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.
  

}
\note{

Reduced--Rank VGLMs  (RR-VGLMs) can be utilized to aid the increasing
number of parameters as \eqn{K} and \eqn{L} grows.
See \code{\link[VGAM:rrvglm]{rrvglm}}.

By default, \eqn{\sigma_1^2, \sigma_2^2}{
\sigma^2[1], \sigma^[2]} and  \eqn{\textrm{Cov}_{12}}{Cov[12]}
are intercept--only. Set argument \code{zero} accordingly to change this.

Package \pkg{tsDyn} also has routines to fit ECMs.
However, the bivariate--ECM handled (similar to that one 
above) differs in their parametrization:
\pkg{tsDyn} considers the \emph{current} estimated residual,
\eqn{\widehat{z}_t}{zhat[t]} instead of 
\eqn{\widehat{z}_{t - 1}}{zhat[t - 1]} in models
\eqn{\Delta y_{1, t} |\Phi_{t - 1}}{\Deltay[1, t] |\Phi[t - 1]} and
\eqn{\Delta y_{2, t} |\Phi_{t - 1}}{\Deltay[2, t] |\Phi[t - 1]}.

See Example 3 below which compares ECMs
fitted with \pkg{VGAMextra} and \pkg{tsDyn}.

}

\references{

Engle, R.F. and Granger C.W.J. (1987)
Co-integration and error correction: Representation, estimation and
testing.
\emph{Econometrica}, \pkg{55}(2), 251--276.


Pfaff, B. (2011)
\emph{Analysis of Integrated and Cointegrated Time Series with \code{R}.}
Seattle, Washington, USA: \emph{Springer}.

}


\author{

 Victor Miranda
 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{

  \code{\link[VGAMextra:MVNcov]{MVNcov}},
  \code{\link[VGAM:rrvglm]{rrvglm}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[VGAM:Links]{Links}},
  \code{\link[VGAM:vglm]{vglm}}.
  
  
}

\examples{

## Example 1. Comparing the Engle -- Granger procedure carried oud by two procedures.
##            ECM.EngleGran() makes easier the fitting process.
## Here, we will use:
## A) The R code 4.2, in Chapter 4, Pfaff (2011).
##    This code 1) generates artificial data and 2) fits an ECM, following
##    the Engle --Granger procedure. 
## B) The ECM.EngleGran() family function to fit the same model assuming
##    bivariate normal innovations. 
## The downside in the R code 4.2 is the assumption of no--correlation among
## the errors. These are generated indenpendently.
## A)
## STEP 1. Set up the data (R code as in Pfaff (2011)).
nn <- 100
set.seed(123456)
e1 <- rnorm(nn)   # Independent of e2
e2 <- rnorm(nn)
y1 <-  cumsum(e1)
y2 <- 0.6 * y1 + e2
lr.reg <- lm(y2 ~ y1)
error <- residuals(lr.reg)
error.lagged <- error[-c(nn - 1, nn)]
dy1 <- diff(y1)
dy2 <- diff(y2)
diff.dat <- data.frame(embed(cbind(dy1, dy2), 2))
colnames(diff.dat) <- c('dy1', 'dy2', 'dy1.1', 'dy2.1')

##  STEP 2. Fit the ECM model, using lm(), R code as in Pfaff (2011).
ecm.reg <- lm(dy2 ~ error.lagged + dy1.1 + dy2.1, data = diff.dat)

summary(ecm.reg)


## B) Now, using ECM.EngleGran() and VGLMs, the steps at A) can be skipped. 
## Enter the I(1)--variables in the response vector only, putting down the
## the dependent variable from the I(1) set, i.e. y2, in the second column.

coint.data <- data.frame(y1 = y1, y2 = y2)
fit.ECM <- vglm(cbind(y1, y2) ~ 1, ECM.EngleGran, data = coint.data, trace = TRUE)

## Check coefficients ##
coef(fit.ECM, matrix = TRUE)  ## Compare 'Diff2' with summary(ecm.reg)
coef(summary(ecm.reg))

head(depvar(fit.ECM))   # The estimated differences (first order)
vcov(fit.ECM)
constraints(fit.ECM, matrix = TRUE)

\dontrun{
### Example 2.  Here, we compare ECM.EngleGran() from VGAMextra with VECM() from
##              package "tsDyn" when fitting an ECM(1, 1). We will make use of
##              the argument 'ordtsDyn' so that the outcomes can be compared.

library("tsDyn")  # Need to be installed first.
fit.tsDyn1 <- with(coint.data, VECM(cbind(y2, y1), lag = 1, estim = "2OLS")) #  MODEL 1
summary(fit.tsDyn1)

### Fit same model using ECM.EngleGran(). NOTE: Set ordtsDyn = 1 !!          #  MODEL 2
fit.ECM.2 <- vglm(cbind(y1, y2) ~ 1, ECM.EngleGran(ecm.order = c(1, 1),
                  resids.pattern = "neither", ordtsDyn = 1),
                  data = coint.data, trace = TRUE)

coef.ECM.2 <- coef(fit.ECM.2, matrix = TRUE)
fit.tsDyn1$coefficients                      ## From pakage 'tsDyn'.
t(coef.ECM.2[, 1:2][c(2, 1, 4, 3), ][, 2:1]) ## FROM VGAMextra 


### Example 3. An ECM(2, 2), with residuals estimated by OLS, with NO intercept
###            and NO trend term. The data set is 'zeroyld', from package tsDyn.
###            ECM.EngleGran() and with VECM() willbe compared again.
data(zeroyld, package = "tsDyn")

# Fit a VECM with Engle-Granger 2OLS estimator:
vecm.eg <- VECM(zeroyld, lag=2, estim = "2OLS") 
summary(vecm.eg)

# For the same data, fit a VECM with ECM.EngleGran(), from VGAMextra.
# Set ordtsDyn = 1 for compatibility! 
fit.ECM.3 <- vglm(cbind(long.run, short.run) ~ 1, ECM.EngleGran(ecm.order = c(2, 2),
                                  resids.pattern = "neither", ordtsDyn = 1),
                  data = zeroyld, trace = TRUE)
coef.ECM.3 <- coef(fit.ECM.3, matrix = TRUE)

#### Compare results
vecm.eg$coefficients                               # From tsDyn
t(coef.ECM.3[, 1:2][c(2, 1, 5, 3, 6, 4 ),][, 2:1]) # FROM VGAMextra
}
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
