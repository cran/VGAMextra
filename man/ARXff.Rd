\name{ARXff}
\alias{ARXff}
\alias{AR1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ VGLTSMs family functions: Order--p Autoregressive Model 
with covariates}
\description{
  Maximum likelihood estimation of the order--p autoregressive model
  (AR(p)) with covariates.
  Estimates the drift, standard deviation (or variance) of the 
  random noise (not necessarily constant),
  and coefficients of the conditional--mean model.
  %of order-\eqn{p} by maximum 
  %likelihood estimation using Fisher scoring.

}
\usage{
      ARXff(order    = 1,
            zero     = c(if (nodrift) NULL else "ARdrift", "ARcoeff"),
            xLag     = 0,
            type.EIM = c("exact", "approximate")[1],
            var.arg  = TRUE, 
            nodrift  = FALSE,
            noChecks = FALSE,
            ldrift   = "identitylink", 
            lsd      = "loge",
            lvar     = "loge",
            lARcoeff = "identitylink",
            idrift   = NULL,
            isd      = NULL,
            ivar     = NULL,
            iARcoeff = NULL)     %-c(3:(order + 2)) 
           
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{order}{
  The order (i.e., 'p') of the AR model, which is recycled if needed.
  See below for further details.
  By default, an autoregressive model of order-\eqn{1} is fitted.
  %By default \code{order = NULL}. Thus, a valid value must be entered
  %initially.
  
  
  }
  \item{zero}{
  Integer or character--strings vector. 
  Name(s) or position(s) of the parameters/linear predictors
  to be modeled as intercept-only. Details at 
  \code{\link[VGAM:zero]{zero}}.
  
  
  }
  \item{xLag}{
  Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  
  }  
  \item{type.EIM, var.arg, nodrift, noChecks}{
  Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  
  }
  \item{ldrift, lsd, lvar, lARcoeff}{
  Link functions applied to the \emph{drift},
  the standar deviation (or variance) of the noise, and the
  AR coefficients.
  Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  
  Further details on 
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  
  }
  \item{idrift, isd, ivar, iARcoeff}{
  Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  
  
  }

  
}

\details{

This family function describes an autoregressive model of order-\eqn{p}
with covariates (ARX(p)). It is a special case of the subclass VGLM--ARIMA
(Miranda and Yee, 2018):

\deqn{ Y_t | \Phi_{t - 1} = \mu_t + \theta_{1} Y_{t - 1} + \ldots + 
\theta_p Y_{t - p} + \varepsilon_t,}{% 
Y[t] | \Phi[t - 1] = \mu[t] + 
\theta[1] Y[t - 1] + \ldots + 
\theta[p] Y[t - p] + e[t],
}

where  \eqn{\boldsymbol{x}_t}{x[t]}
a (possibly time--varying) covariate vector and
\eqn{\mu_t = \mu^{\star} + \boldsymbol{\beta}^T \boldsymbol{x}_t}{
\mu[t] = \mu[*] + \beta^T * x[t]
} is a (time--dependent) scaled--mean, known as \emph{drift}.

At this stage, conditional Gaussian white noise, 
\eqn{\varepsilon_t| \Phi_{t - 1}}{e[t] | \Phi [t - 1]} 
is handled, in the form
\deqn{\varepsilon_t | \Phi_{t - 1} \sim 
  N(0, \sigma^2_{\varepsilon_t | \Phi_{t - 1}}).}{
   e[t] | \Phi[t - 1] ~ N(0 , \sigma[e]^2 | \Phi [t - 1] ). }
   The distributional assumptions on the observations are then
   
   \deqn{Y_t | \Phi_{t - 1} \sim 
   N(\mu_{t | \Phi_{t - 1}}, 
   \sigma^2_{\varepsilon_t | \Phi_{t - 1}}), }{
   Y[t] | \Phi[t - 1] ~ N(\mu[t] | \Phi [t - 1] ,
   \sigma[e]^2 | \Phi [t - 1]  ),}
   involving the conditional mean 
   equation for the ARX(p) model:
\deqn{\mu_{t | \Phi_{t - 1}} = \mu_t + 
 \boldsymbol{\beta}^T * \boldsymbol{x}_t
   \theta_{1} Y_{t - 1} + \ldots + 
      \theta_p Y_{t - p}.}{
      \mu[t] | \Phi [t - 1] = \mu[t] + 
      \beta^T x_t + \theta[1] Y[t - 1] + 
      \ldots + \theta[p] Y[t - p].
      }

\eqn{\Phi_{t}}{\Phi[t]} denotes the information
of the joint process
\eqn{\left(Y_{t}, \boldsymbol{x}_{t + 1}^T \right)}{
(Y[t], x[t + 1])
},
at time \eqn{t}.

%With no covariates,
%\eqn{\mu_t = \mu^{\star}}{\mu[t] = \mu[*]} 
%is a scaled mean known as \emph{drift}.

The loglikelihood is computed by \code{\link[VGAMextra:dARp]{dARp}},
at each Fisher scoring iteration.

The linear predictor is
\deqn{\boldsymbol{\eta} = \left( 
\mu_t, \log \sigma^{2}_{\varepsilon_{t | \Phi_{t - 1}}},
\theta_1, \ldots, \theta_p
\right)^T.}{
\eta = (\mu[t], log \sigma[e]^2 | \Phi[t - 1],
\theta[1], \ldots, \theta[p] )^T.} 

Note, the covariates may also intervene in the conditional
variance model
\eqn{\log \sigma^{2}_{\varepsilon_{t | \Phi_{t - 1}}}.}{
log \sigma[e]^2 | \Phi[t - 1].
} Hence, this family function
does not restrict the noise to be \emph{strictly} white noise
(in the sense of \emph{constant variance}). 


The unconditional mean, 
\eqn{ E(Y_{t}) = \mu}{%
E(Y[t]) = \mu }, satisfies
\deqn{\mu \rightarrow \frac{\mu^{\star}}{1 - (\theta_1 +
\ldots + \theta_p)} }{%
\mu -> \mu[*] / ( 1 - \theta[1] - \ldots - \theta[p] )
}
when the process is stationary, and no covariates are involved.


This family function currently handles
\emph{multiple responses} so that a matrix can be used as the response.
Also, for further details on VGLM/VGAM--link functions refer to 
\code{\link[VGAM:Links]{Links}}.

Further choices for the random noise, besides Gaussian, will be
implemented over time.


}

\value{

An object of class \code{"vglmff"} 
(see \code{\link[VGAM:vglmff-class]{vglmff-class}}). The
object is used by VGLM/VGAM modelling functions, such as
\code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.


%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{

Madsen, H. (2008)
Time Series Analysis. 
Florida, USA: \emph{Chapman & Hall}(Sections 5.3 and 5.5).


Porat, B., and Friedlander, B. (1986)
Computation of the Exact Information Matrix of Gaussian
Time Series with Stationary Random Components.
\emph{IEEE Transactions on Acoustics, Speech and Signal
Processing}. \bold{ASSp-34(1)}, 118--130.

}

\author{

 Victor Miranda and T. W. Yee

}


\section{Note}{

\code{zero} can be either an \emph{integer} vector
or a vector of \bold{character strings}
specifying either the position(s) or name(s) (partially or not) of the 
parameter(s) modeled as intercept-only. Numeric values can be set as usual 
(See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}).
Character strings can be entered as per
parameter names in this family function, given by:

\code{c("drift", "noiseVar" or "noiseSD", "ARcoeff")}.

Users can modify the \code{zero} argument according to their needs.



By default, \eqn{\mu_t}{\mu[t]} and the coefficients
\eqn{\theta_1, \ldots, \theta_p}{\theta[1], \ldots, \theta[p]}
are intercept--only. That is, 
\eqn{\log \sigma^{2}_{\varepsilon_{t | \Phi_{t - 1}}}}{
log \sigma[e]^2 | \Phi[t - 1]} is modelled in terms of any explanatories
entered in the \code{formula}.

Users, however, can modify this
according to their needs via
{\code{\link[VGAM:zero]{zero}}}. For example, set the covariates in the
\code{drift} model, \eqn{\mu_t}{\mu[t]}.
In addition, specific constraints 
for parameters are handled through the function 
\code{\link[VGAMextra:cm.ARMA]{cm.ARMA}}.


If \code{var.arg = TRUE}, this family function estimates
\eqn{\sigma_{\varepsilon_t | \Phi_{t - 1}}^2}{ \sigma[e]^2 | 
\Phi[t - 1]}.
Else, the
\eqn{\sigma_{\varepsilon_t | \Phi_{t - 1}}}{ \sigma[e] |\Phi[ t- 1]}
estimate is returned.
%By default, an autoregressive model of order-\eqn{1} is fitted.

For this family function the \code{order} is recycled. That is,
\code{order} will be replicated up to the number of responses
given in the \code{vglm} call is matched.

}

\section{Warning}{


Values of the estimates may not correspond
to stationary ARs, leading to low accuracy in the MLE estimates,
e.g., values very close to 1.0.
\emph{Stationarity} is then examined, via
\code{\link[VGAMextra:checkTS.VGAMextra]{checkTS.VGAMextra}},
if
\code{noChecks = FALSE} (default) 
and \bold{no} constraint matrices are set
(See \code{\link[VGAM:constraints]{constraints}} for further
details on this).
If the estimated model very close to be non-stationary, then
a \code{warning} will be outlined. 
Set \code{noChecks = TRUE} to completely ignore this.


NOTE: Full details on these 'checks' are shown within the 
\code{summary()} output.


}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

  \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}},
  \code{\link[VGAMextra:ARMAXff]{ARMAXff}},
  \code{\link[VGAMextra:MAXff]{MAXff}},
  \code{\link[VGAMextra:checkTS.VGAMextra]{checkTS.VGAMextra}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[VGAM:Links]{Links}},
  \code{\link[VGAM:vglm]{vglm}},
  %\code{\link[VGAMextra:cm.ARMA]{cm.ARMA}},


}
\examples{
set.seed(1)
nn     <- 150
tsdata <- data.frame(x2 =  runif(nn))             # A single covariate.
theta1 <- 0.45; theta2 <- 0.31; theta3 <- 0.10     # Coefficients
drift  <- c(1.3, -1.1)                             # Two responses.
sdAR   <- c(sqrt(4.5), sqrt(6.0))                  # Two responses.

# Generate AR sequences of order 2 and 3, under Gaussian noise.
# Note, the drift for 'TS2' depends on x2 !
tsdata  <-  data.frame(tsdata, TS1 = arima.sim(nn, 
              model = list(ar = c(theta1, theta1^2)),  rand.gen = rnorm, 
              mean = drift[1], sd = sdAR[1]),
                                TS2 = arima.sim(nn,  
              model = list(ar = c(theta1, theta2, theta3)), rand.gen = rnorm, 
              mean = drift[2] + tsdata$x2 , sd = sdAR[2]))

# EXAMPLE 1. A simple AR(2), maximizing the exact log-likelihood
# Note that parameter constraints are involved for TS1, but not 
# considered in this fit. "rhobit" is used as link for AR coeffs.

fit.Ex1 <- vglm(TS1 ~ 1, ARXff(order = 2, type.EIM = "exact",
                      #iARcoeff = c(0.3, 0.3, 0.3), # OPTIONAL INITIAL VALUES
                      # idrift = 1, ivar = 1.5, isd = sqrt(1.5),
                      lARcoeff = "rhobit"), 
              data = tsdata,  trace = TRUE, crit = "loglikelihood")
Coef(fit.Ex1)
summary(fit.Ex1)
vcov(fit.Ex1, untransform = TRUE)       # Conformable with this fit.
AIC(fit.Ex1)
#------------------------------------------------------------------------#
# Fitting same model using arima(). 
#------------------------------------------------------------------------#
(fitArima <- arima(tsdata$TS1, order = c(2, 0, 0)))
# Compare with 'fit.AR'. True are theta1 = 0.45; theta1^2 = 0.2025
Coef(fit.Ex1)[c(3, 4, 2)]    # Coefficients estimated in 'fit.AR'


# EXAMPLE 2. An AR(3) over TS2, with one covariate affecting the drift only.
# This analysis makes sense as the TS2's drift is a function ox 'x2', 
# i.e., 'x2' affects the 'drift' parameter only. The noise variance 
# (var.arg = TRUE) is estimated, as intercept-only. See the 'zero' argument.

#------------------------------------------------------------------------#
# This model CANNOT be fitted using arima()
#------------------------------------------------------------------------#
fit.Ex2 <- vglm(TS2 ~ x2,  ARXff(order = 3, zero = c("noiseVar", "ARcoeff"), 
                             var.arg = TRUE), 
                  ## constraints = cm.ARMA(Model = ~ 1, lags.cm = 3, Resp = 1),
              data = tsdata,  trace = TRUE, crit = "log")

# True are theta1 <- 0.45; theta2 <- 0.31; theta3 <- 0.10
coef(fit.Ex2, matrix = TRUE)
summary(fit.Ex2)    
vcov(fit.Ex2)
BIC(fit.Ex2)
constraints(fit.Ex2)


# EXAMPLE 3. Fitting an ARX(3) on two responses TS1, TS2; intercept-only model with
#  constraints over the drifts. Here, 
# a) No checks on invertibility performed given the use of cm.ARMA().
# b) Only the drifts are modeled in terms of 'x2'. Then,  'zero' is
# set correspondingly.
#------------------------------------------------------------------------#
# arima() does not handle this model.
#------------------------------------------------------------------------#
fit.Ex3 <- vglm(cbind(TS1, TS2) ~ x2, ARXff(order = c(3, 3), 
                            zero = c("noiseVar", "ARcoeff"), var.arg = TRUE), 
             constraints = cm.ARMA(Model = ~ 1 + x2, lags.cm = c(3, 3), Resp = 2),
             trace = TRUE, data = tsdata, crit = "log")

coef(fit.Ex3, matrix = TRUE)
summary(fit.Ex3)
vcov(fit.Ex3)
constraints(fit.Ex3)


%}

}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
