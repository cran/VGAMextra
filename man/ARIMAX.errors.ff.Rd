\name{ARIMAX.errors.ff}
\alias{ARIMAX.errors.ff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ VGLTSMs Family Functions:
Generalized integrated regression
with order--\eqn{(p, q)} ARMA errors}

\description{

  A VLTSMff for dynamic regression.
  Estimates regression models with order--\eqn{(p, d, q)} ARIMA
  errors by maximum likelihood.
  %Differencing of the series (response) 
  %is handled up to order \eqn{d}.
}
\usage{
      ARIMAX.errors.ff(order = c(1, 1, 1),
                       zero = "var",  # optionally, "mean".
                       order.trend = 0,
                       include.int = TRUE,
                       diffCovs  = TRUE,
                       xLag = 0,
                       include.currentX = TRUE,
                       lvar = "loglink",
                       lmean = "identitylink") 
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{order}{
  The usual \eqn{(p, d, q)} integer vector as in, e.g., 
  \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}. 
  By default, an order--\eqn{(p, q)} ARMA model is 
  fitted on the errors, whlist \eqn{d} is the 
  degree of differencing on the response.

  
  }
  \item{zero}{
  What linear predictor is modelled as intercept--only?
  See \code{\link[VGAM:zero]{zero}} and
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}
  for further details.
  
  
  }
  \item{order.trend}{
   Non--negative integer. Allows to incorporate a polynomial
   trend of order \code{order.trend} in the forecast mean
   function. 
  
  
  }
  \item{include.int}{
  Logical. Should an intercept (int) be included in
  the model for \eqn{y_t}{y[t]}? Default is
  \code{TRUE}. See below for details.
  
  
  }
  \item{diffCovs}{
  Logical. If \code{TRUE} (default)
  the order--\code{d} difference
  of the covariates is internally computed and then
  incorporated in the regression model. Else, 
  only the current values are included.
  
  
  }
  \item{xLag}{
   Integer. If entered, the covariates, say 
   \eqn{\boldsymbol{x}_t}{x[t]} are laggeg 
   (up to order \code{xLag}) and then embedded
   in the regression model. 
   See below for further details.
   
  }
  \item{include.currentX}{
  Logical. If \code{TRUE}, the \emph{actual} values, 
  \eqn{\boldsymbol{x}_t}{x[t]}, are included in the regression
  model. Else, this is ignored and only the lagged
  \eqn{\boldsymbol{x}_{t - 1}, \ldots, \boldsymbol{x}_{t - xLag}}{
  x[t - 1], \ldots, x[t - xLag]
  } will be included.
  
  
  }
  \item{lvar, lmean}{
  Link functions applied to conditional mean and the variance.
  Same as \code{\link[VGAM:uninormal]{uninormal}}.
  
  }
}

\details{

The generalized linear regression model with ARIMA errors
is another subclass of VGLTSMs (Miranda and Yee, 2018).

For a univariate time series, say \eqn{y_t}{y[t]},
and a \eqn{p}--dimensional vector of covariates
\eqn{\boldsymbol{x}_t}{x[t]} covariates,
the model described by this VGLTSM family function is
\deqn{
 y_t = \boldsymbol{\beta}^T \boldsymbol{x}_t + u_t,
}{
y[t] = \beta^T x_t + u[t],
}
\deqn{
u_t = \theta_1 u_{t - 1} + \cdots + \theta_p u_{t - p} + z_t +
    \phi_1 z_{t - 1} + \cdots + \phi_1 z_{t - q}.
}{
u[t] = \theta[1] u[t - 1] + \ldots +\theta[p] u[t - p] + z[t] +
\phi[1] z[t - 1] + \ldots +\phi[q] u[t - q].
}
%

The first entry in \eqn{x_t} equals 1, allowing 
an intercept, for every $t$. Set 
\code{include.int = FALSE} to set this to zero,
dimissing the intercept.

Also, if \code{diffCovs = TRUE}, then the differences up to order
\code{d} of the set
\eqn{\boldsymbol{x}_t}{x[t]} are embedded in the model 
for \eqn{y_t}{y[t]}.
If \code{xLag}\eqn{> 0}, the lagged values up to order
\code{xLag} of the covariates are also included.
%regardless the
%differencing degree \eqn{d} entered in \code{order}.

The random disturbances \eqn{z_t}{z[t]} are by default
handled as \eqn{N(0, \sigma^2_z)}{N(0, \sigma[z]^2)}. Then,
denoting \eqn{\Phi_{t}}{\Phi[t]} as the history of the
process \eqn{(x_{t + 1}, u_t)}{(x[t + 1], u[t])} up to time \eqn{t}, yields
%
\deqn{ E(y_t | \Phi_{t - 1}) = \boldsymbol{\beta}^T \boldsymbol{x}_t +
\theta_1 u_{t - 1} + \cdots + \theta_p u_{t - p} +
    \phi_1 z_{t - 1} + \cdots + \phi_1 z_{t - q}. }{
E(y[t] | \Phi[t - 1]) =   \beta^T x[t] +
\theta[1] u[t - 1] + \ldots +\theta[p] u[t - p] + 
\phi[1] z[t - 1] + \ldots +\phi[q] u[t - q].
}


Denoting \eqn{\mu_t =  E(y_t | \Phi_{t - 1}),}{
\mu[t] = E(y[t] | \Phi[t - 1]), }
the default linear predictor for this VGLTSM family function is
\deqn{ \boldsymbol{\eta} = ( \mu_t, \log \sigma^2_{z})^T.}{ 
\eta = (\mu[t], log \sigma[t]^2)^T.}




}

\value{

An object of class \code{"vglmff"} 
(see \code{\link[VGAM:vglmff-class]{vglmff-class}})
 to be used by VGLM/VGAM modelling functions, e.g.,
\code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.
  

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%\references{
%Miranda, V. and Yee, T.W. (2018)
%Vector Generalized Linear Time Series Models.
%\emph{In preparation}.
%}

\author{

 Victor Miranda 

}


\note{

If \code{d = 0} in \code{order}, then \code{ARIMAX.errors.ff}
will perform as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.



}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[VGAM:uninormal]{uninormal}},
  \code{\link[VGAM:vglm]{vglm}}.


}

\examples{
### Estimate a regression model with ARMA(1, 1) errors.
## Covariates are included up to lag 1.
set.seed(20171123)
nn <- 250
x2 <- rnorm(nn)                                    # One covariate
sigma2 <- exp(1.15); theta1 <- 0.5; phi1 <- 0.27   # True coefficients
beta0  <- 1.25; beta1 <- 0.25; beta2 <- 0.5

y <- numeric(nn)
u <- numeric(nn)
z <- numeric(nn)

u[1] <- rnorm(1)
z[1] <- rnorm(1, 0, sqrt(sigma2))

for(ii in 2:nn) {
  z[ii] <- rnorm(1, 0, sqrt(sigma2))
  u[ii] <- theta1 * u[ii - 1] + phi1 * z[ii - 1] + z[ii]
  y[ii] <- beta0 + beta1 * x2[ii] + beta2 * x2[ii - 1] + u[ii]
}

# Remove warm-up values.
x2 <- x2[-c(1:100)]
y  <- y[-c(1:100)]

\donttest{
plot(ts(y), lty = 2, col = "blue", type = "b")
abline(h = 0, lty = 2)
}

## Fit the model.
ARIMAX.reg.fit <- vglm(y ~ x2, ARIMAX.errors.ff(order = c(1, 0, 1), xLag = 1),
             data = data.frame(y = y, x2 = x2), trace = TRUE)
coef(ARIMAX.reg.fit, matrix = TRUE)
summary(ARIMAX.reg.fit, HD = FALSE)

\donttest{
# Compare to arima()
# arima() can't handle lagged values of 'x2' by default, but these 
# may entered at argument 'xreg'.
arima(y, order = c(1, 0, 1), xreg = cbind(x2, c(0, x2[-150])))
}


}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
