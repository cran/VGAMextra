\name{ARIMAXff}
\alias{ARIMAXff}
\alias{ARMAXff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ VGLTSMs Family functions: 
The Order--\eqn{(p, d, q)} Autoregressive Integrated Moving
Average Model (ARIMA(p, d, q)) with covariates}
\description{

  Maximum likelihood estimation of the drift, 
  standard deviation or variance of the random noise,
  and coefficients of an autoregressive integrated
  moving average process of order-\eqn{(p, d, q)}
  with covariates by MLE using Fisher scoring.
  No seasonal terms handled yet. No seasonal components
  handled yet.

}
\usage{
      ARIMAXff(order   = c(1, 1, 0),
               zero     = c("ARcoeff", "MAcoeff"),
               diffCovs = TRUE,
               xLag     = 0,
               include.current = FALSE,
               type.EIM = c("exact", "approximate")[1], 
               var.arg  = TRUE,
               nodrift  = FALSE,
               noChecks = FALSE,
               ldrift   = "identitylink", 
               lsd      = "loglink",
               lvar     = "loglink",
               lARcoeff = "identitylink",
               lMAcoeff = "identitylink", 
               idrift   = NULL,
               isd      = NULL,
               ivar     = NULL,
               iARcoeff = NULL, 
               iMAcoeff = NULL) 
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{order}{
  Integer vector with three components, 
  (\code{p, d, q}): The AR order 
  (\code{p}), the degree of differencing (\code{d}),
  and the MA order (\code{q}).
  
  }
  \item{zero}{
  Integer or character--strings vector. 
  Name(s) or position(s) of the parameters/linear predictors
  to be modeled as intercept-only. Details at 
  \code{\link[VGAM:zero]{zero}}.
  
  
  }
  \item{diffCovs}{
  Logical. The default is \code{diffCovs = TRUE}, which means that
  the order--\code{d} differences
  of the covariates (if entered) are internally computed and then
  incorporated in the conditional--mean model. Otherwise,
  only the current \code{actual} values of the covariates are included.
  
  
  }
  \item{xLag}{
  Integer, non--negative. If \code{xLag} > 0, the covariates
  at \emph{current} time, \eqn{\boldsymbol{x}_t}{x[t]}, plus
  the lagged covariates up to order \code{'xLag'} are
  embedded into the design matrix (as covariates too).
  Leave \code{xLag = 0} and only the
  \emph{current} value, \eqn{\boldsymbol{x}_t}{x[t]},
  will be considered. See more details below.
  
  }
  \item{include.current}{
  Logical. Same as 
  \code{\link[VGAMextra:ARIMAX.errors.ff]{ARIMAX.errors.ff}}.
  
  
  }
  \item{type.EIM}{
  The type of expected information matrix (EIM) of the ARMA process
  to be utilized in Fisher scoring.
  \code{type.EIM = "exact"} (default)
  enables the \emph{exact} IM (Porat, et.al. 1986),
  otherwise the \code{approximate} version is utilized.
  
  
  }
  \item{var.arg}{
  Logical. If \code{FALSE} (default), then the standard deviation
  of the random noise is estimated. Else, the variance estimate is 
  returned.
  
  
  }
  \item{nodrift}{
  Logical. \code{nodrift = TRUE} supresses estimation of the 
  intercept (the \emph{drift} in the ARMA case),
  which is set to zero internally.
  
  
  }
  \item{noChecks}{
  Logical. If \code{FALSE} (default), this family function internally
  checks \emph{stationarity} (AR case) and
  \emph{invertibility} (MA case) of the the estimated model.
  A \code{warning} is correspondingly displayed.
  
  
  }
  \item{ldrift, lsd, lvar, lARcoeff, lMAcoeff}{
  Link functions applied to the intercept, 
  the random noise standard deviation (or optionally,
  the variance), and the coefficients in the 
  ARMA--type conditional--mean model.
  
  
  }
  \item{idrift, isd, ivar, iARcoeff, iMAcoeff}{
  Optional initial values for the intercept (drift), noise SD
  (or variance), and ARMA coeffcients (a vector of length \eqn{p + q}). 
  If failure to converge occurs then try different values and monitor 
  convergence by using \code{trace = TRUE} in the \code{vglm()} call.
  
  
  }
  
}

\details{

Let \eqn{\boldsymbol{x}_t}{x[t]} be a (probably time--varying) vector of 
suitable covariates. The ARIMAX model handled by \code{ARIMAXff} is

\deqn{ \nabla^d Y_t = \mu^{\star} + \boldsymbol{\beta}^T 
              \nabla^d \boldsymbol{x}_t +
             \theta_1 \nabla^d Y_{t - 1} + \cdots +
             \theta_p \nabla^d Y_{t - p} +
             \phi_1 \varepsilon_{t - 1} + \cdots + 
             \phi_q \varepsilon_{t - q} + \varepsilon, }{%
    D(Y[t], d) = \mu[*] +  \beta^T * D(x[t], d)  +
             \theta[1] D(Y[t], d)[t - 1] + \ldots +
             \theta[p] D(Y[t], d)[t - p] +
             \phi[1] e[t - 1] + \ldots + 
             \phi[q] e[t - q] + e[t],
  }
  
  with \eqn{\nabla^d (\cdot)}{D( * , d)} the operator differencing a
  series \code{d} times. If \code{diffCovs = TRUE}, this function
  differencing the covariates \code{d} times too. 
  
Similarly, \code{ARMAXff} manages

\deqn{ \nabla^d Y_t = \mu^{\star} + \boldsymbol{\beta}^T \boldsymbol{x}_t +
             \theta_1 Y_{t - 1} + \cdots + \theta_p Y_{t - p} +
             \phi_1 \varepsilon_{t - 1} + \cdots + 
             \phi_q \varepsilon_{t - q} + \varepsilon, }{%
      Y[t] = \mu[*] +  \beta^T * x_t +
             \theta[1] Y[t - 1] + \ldots + \theta[p] Y[t - p] +
             \phi[1] e[t - 1] + \ldots + 
             \phi[q] e[t - q] + e[t],
  }

where
\deqn{    \varepsilon_{t | \Phi_{t - 1}}
         \sim N(0, \sigma_{\varepsilon_t | \Phi_{t - 1}}^2).}{% 
                       e[t] |\Phi[t - 1] ~ N(0, \sigma[e]^2 | \Phi[t - 1]).
                       } 
Note, \eqn{\sigma_{\varepsilon | \Phi_{t - 1}}^2}{%
\sigma[e]^2 | \Phi[t - 1]} is \emph{conditional} on 
\eqn{\Phi_{t - 1}}{\Phi[t - 1]}, 
the information
of the joint process
\eqn{\left(Y_{t - 1}, \boldsymbol{x}_t \right)}{
(Y[t - 1], x[t])
},
at time \eqn{t}, and hence may be 
modelled in terms of \eqn{\boldsymbol{x}_t}{x[t]},
if required.


\code{ARIMAXff()} and \code{ARMAXff()}
handle \emph{multiple responses}, thus a
matrix can be used as the response. 
Note, no seasonal terms handled. This feature is to be
incorporated shortly.

The default linear predictor is
\deqn{\boldsymbol{\eta} = \left( 
\mu, \log \sigma^{2}_{\varepsilon_{t | \Phi_{t - 1}}},
\theta_1, \ldots, \theta_p,
\phi_1, \ldots, \phi_q
\right)^T.}{
\eta = (\mu, log \sigma[e]^2 | \Phi[t - 1],
\theta[1], \ldots, \theta[p],
\phi[1], \ldots, \phi[q] )^T.}
Other links are also handled. See \code{\link[VGAM:Links]{Links}}.

Further choices for the random noise, besides Gaussian, will
be implemented over time.

As with \code{\link[VGAMextra:ARXff]{ARXff}}
and \code{\link[VGAMextra:MAXff]{MAXff}}, 
choices for the EIMs are \code{"exact"} and \code{"approximate"}.
Covariates may be incorporated in the fit for any linear
predictor above. Hence, \code{ARIMAXff} supports non--stationary 
processes (\eqn{\sigma_{\varepsilon_t | \Phi_{t - 1}}^2}{
\sigma[e]^2 | \Phi[t - 1]}) may depend on 
\eqn{\boldsymbol{X}_t}{X[t]}.
Also, constraint matrices 
on the linear predictors may be
entered through \code{\link[VGAMextra:cm.ARMA]{cm.ARMA}} or
using the argument \code{constraints}, from
\code{\link[VGAM:vglm]{vglm}}.

Checks on stationarity and
invertibility on the esitmated process are performed by default.
Set \code{noChecks = TRUE} to dismiss this step.


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

\references{


Miranda, V. and Yee, T.W. (2018)
Vector Generalized Linear Time Series Models.
\emph{In preparation}.


Porat, B., and Friedlander, B. (1986)
Computation of the Exact Information Matrix of Gaussian
Time Series with Stationary Random Components.
\emph{IEEE Transactions on Acoustics, Speech and Signal
Processing}. \bold{ASSp-34(1)}, 118--130.

}

\author{

 Victor Miranda and T. W. Yee

}


\note{
  
  No seasonal components handled yet.
 
 If no covariates, \eqn{\boldsymbol{x}_t}{x[t]},
 are incorporated in the analysis,
 then {\code{ARIMAXff}} fits an ordinary ARIMA model.
 Ditto with \code{ARMAXff}.

If \code{nodrift = TRUE}, then the 'drift' is removed from the
vector of parameters and is not estimated.

By default, an ARMA model of order--\eqn{c(1, 0)} with
order--1 differences is fitted. When initial
values are entered (\code{isd}, \code{iARcoeff}, etc.),
they are recycled
according to the number of responses.

Also, the ARMA coefficients 
are intercept--only (note, \code{zero = c("ARcoeff", }
\code{ "MAcoeff")})
This may altered via \code{zero}, or by
constraint matrices (See \code{\link[VGAM:constraints]{constraints}})
using \code{\link[VGAMextra:cm.ARMA]{cm.ARMA}}.


Checks on stationarity and/or
invertibility can be manually via
\code{\link[VGAMextra:checkTS.VGAMextra]{checkTS.VGAMextra}}.

}

\section{Warning}{


\code{zero} can be a \bold{numeric} or a \bold{character--strings} 
vector specifying the position(s) or the name(s) of 
the parameter(s) modeled as intercept--only.
Numeric values can be set as 
usual (See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}).
If names are entered, the parameter names in this 
family function are:

\code{c("drift.mean", "noiseVar" || "noiseSD", "ARcoeff", "MAcoeff")}.

Manually modify this if required.
For simplicity, the second choice is recommended.


}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[VGAMextra:ARXff]{ARXff}},
  \code{\link[VGAMextra:MAXff]{MAXff}},
  \code{\link[VGAMextra:checkTS.VGAMextra]{checkTS.VGAMextra}},
  \code{\link[VGAMextra:cm.ARMA]{cm.ARMA}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[VGAM:constraints]{constraints}},
  \code{\link[VGAM:vglm]{vglm}}.


}

\examples{
set.seed(3)
nn      <- 90
theta   <- c(0.12, 0.17)  # AR coefficients
phi     <- c(-0.15, 0.20)  # MA coefficients.
sdWNN <- exp(1.0)          # SDs
mu    <- c(1.25, 0.85)     # Mean (not drift) of the process.
covX  <- runif(nn + 1)     # A single covariate.
mux3  <- mu[1] + covX
##
## Simulate ARMA processes. Here, the drift for 'tsd3' depends on covX.
##
tsdata <- data.frame(TS1 = mu[1] + arima.sim(model = list(ar = theta, ma = phi,
                             order = c(2, 1, 2)), n = nn, sd = sdWNN ),
                      TS2 = mu[2] + arima.sim(model = list(ar = theta, ma = phi,
                            order = c(2, 1, 2)), n = nn, sd = exp(2 + covX)),
                      TS3 =  mux3 + arima.sim(model = list(ar = theta, ma = phi,
                            order = c(2, 1, 2)), n = nn, sd = exp(2 + covX) ),
                      x2 = covX)

### EXAMPLE 1. Fitting a simple ARIMA(2, 1, 2) using vglm(). 
# Note that no covariates involved. 
fit.ARIMA1 <-  vglm(TS1 ~ 1, ARIMAXff(order = c(2, 1, 2), var.arg = FALSE,
                       # OPTIONAL INITIAL VALUES
                       # idrift = c(1.5)*(1 - sum(theta)), 
                       # ivar = exp(4), isd = exp(2),
                       # iARcoeff = c(0.20, -0.3, 0.1),
                       # iMAcoeff = c(0.25, 0.35, 0.1),
                        type.EIM = "exact"), 
                 data = tsdata, trace = TRUE, crit = "log")
Coef(fit.ARIMA1)
summary(fit.ARIMA1)
vcov(fit.ARIMA1, untransform = TRUE)
#------------------------------------------------------------------------#
# Fitting same model using arima().
#------------------------------------------------------------------------#
# COMPARE to EXAMPLE1 
( fitArima  <- arima(tsdata$TS1, order = c(2, 1, 2)) ) 


### EXAMPLE 2. Here only the ARMA coefficients and drift are intercept-only.
# The random noise variance is not constant.
fit.ARIMA2 <-  vglm(TS2 ~ x2,  ARIMAXff(order = c(2, 1, 2), var.arg = TRUE,
                     lARcoeff = "rhobitlink", lMAcoeff = "identitylink",
                     type.EIM = c("exact", "approximate")[1], 
                     # NOTE THE ZERO ARGUMENT. 
                     zero = c("drift.mean", "ARcoeff", "MAcoeff")),
              data = tsdata, trace = TRUE)

coef(fit.ARIMA2, matrix = TRUE)
summary(fit.ARIMA2)
constraints(fit.ARIMA2)



### EXAMPLE 3. Here only ARMA coefficients are intercept-only.
# The random noise variance is not constant.
# Note that the "drift" and the "variance" are "generated" in 
# terms of 'x2' above for TS3.

fit.ARIMA3 <- vglm(TS3 ~ x2,  ARIMAXff(order = c(1, 1, 2), var.arg = TRUE,
                     lARcoeff = "identitylink", lMAcoeff = "identitylink",
                     type.EIM = c("exact", "approximate")[1], nodrift = FALSE,
                     zero = c( "ARcoeff", "MAcoeff")), # NOTE THE ZERO ARGUMENT. 
              data = tsdata, trace = TRUE)
              
coef(fit.ARIMA3, matrix = TRUE)
summary(fit.ARIMA3)
constraints(fit.ARIMA3)

}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
