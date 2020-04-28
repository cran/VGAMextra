\name{MAXff}
\alias{MAXff}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{VGLTSMs
Family Functions: Order--q Moving Average Model with covariates }
\description{


  Estimates the intercept, standard deviation (or variance) of the 
  random noise (not necessarily constant), and
  the conditional--mean model coefficients of an
  order--q moving average (MA) process with covariates (MAX(q))
  by maximum likelihood estimation using Fisher scoring.
  

}
\usage{
      MAXff(order    = 1,
            zero     = c(if (nomean) NULL else "Mean", "MAcoeff"),
            xLag     = 0,
            type.EIM = c("exact", "approximate")[1],
            var.arg  = TRUE, 
            nomean   = FALSE,
            noChecks = FALSE,
            lmean    = "identitylink", 
            lsd      = "loglink",
            lvar     = "loglink",
            lMAcoeff = "identitylink",
            imean    = NULL,
            isd      = NULL,
            ivar     = NULL,
            iMAcoeff = NULL)      % -c(3:(order + 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{order}{
  The order 'q' of the MA model, recycled if needed.
  By default \code{q = 1}.
  
  
  }
  \item{zero}{
  Integer or character--string vector.
  Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  Details at \code{\link[VGAM:zero]{zero}}.
  
  }
  \item{xLag}{
  Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  
  }
  \item{type.EIM, var.arg}{
  Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  
  }
  \item{nomean}{
  Logical. \code{nomean = TRUE} supresses estimation of the 
  \emph{mean} (intercept of the conditional--mean model).
  
  
  }
  \item{noChecks}{
  Logical. Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
  
  
  }
  \item{lmean, lsd, lvar, lMAcoeff}{
  Link functions applied to the mean (intercept),
  the standard deviation or variance of the random noise, and
  the MA coefficients (conditional--mean model).
   Note, \code{lmean} plays the role of \code{ldrift}.
  
  
  }
  \item{imean, isd, ivar, iMAcoeff}{
    Same as \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}}.
    Note, \code{imean} plays the role of \code{idrift}.
  
  }
}

\details{
Similar to \code{\link[VGAMextra:ARXff]{ARXff}},
this family function fits an order--\eqn{q} moving average model
with covariates (MAX(q)),
another special case of the class VGLM--ARIMA
(Miranda and Yee, 2018). Observations, \eqn{Y_t}{Y[t]}, are
modelled as

\deqn{ Y_t | \Phi_{t - 1} = \mu_t + \phi_{1} \varepsilon_{t - 1} +
\ldots + \phi_q \varepsilon_{t - q} + \varepsilon_t,}{% 
Y[t] | \Phi[t - 1] = \mu[t] + \phi[1] e[t - 1] + \ldots + 
\phi[q] e[t - p] + e[t],
}

where \eqn{\mu_t}{\mu[t]} is the (possibly time--dependent) intercept,
modelled as
\eqn{\mu_t = \mu + \boldsymbol{\beta}^T \boldsymbol{x}_t}{
mu[t] = \mu + \beta^T x[t]}, and
the errors are mean--zero Gaussian:
\eqn{\varepsilon_t | \Phi_{t - 1} \sim 
  N(0, \sigma^2_{\varepsilon_t | \Phi_{t - 1}})}{
   e[t] | \Phi[t - 1] ~ N(0 , \sigma[e]^2 | \Phi [t - 1] ) }.
   The symbol
    \eqn{\Phi_{t}}{\Phi[t]} denotes the history
of the joint process
\eqn{\left(Y_{t}, \boldsymbol{X}_{t + 1}^T \right)}{
(Y[t], X[t + 1])
},
at time \eqn{t} for a time--varying covariate 
vector \eqn{\boldsymbol{x}_t}{x[t]}.


At each step of Fisher scoring, the \emph{exact} log-likelihood
based on model above 
is computed through \code{\link[VGAMextra:dMAq]{dMAq}}.


The linear predictor by default is 
\deqn{\boldsymbol{\eta} = \left( 
\mu_t, \log \sigma^{2}_{\varepsilon_{t | \Phi_{t - 1}}},
\phi_1, \ldots, \phi_q
\right)^T.}{
\eta = (\mu[t], log \sigma[e]^2 | \Phi[t - 1],
\phi[1], \ldots, \phi[q] )^T.} 



The unconditional mean of the process is simply
\eqn{ E(Y_{t}) = \mu}{%
E(Y[t]) = \mu }, provided no covariates added.

This family function
is not restricted to the noise to be strictly white noise
(in the sense of \emph{constant variance}). 
That is, covariates may be incorporated in the linear predictor
\eqn{\log \sigma^{2}_{\varepsilon_{t | \Phi_{t - 1}}}.}{
log \sigma[e]^2 | \Phi[t - 1].
} Also, it handles
\emph{multiple responses} so that a matrix can be used as the response.
For further details on VGLM/VGAM--link functions,
such as \code{\link[VGAM:logit]{logit}}, refer to 
\code{\link[VGAM:Links]{Links}}.

}
\value{

 An object of class \code{"vglmff"} 
(see \code{\link[VGAM:vglmff-class]{vglmff-class}}) to be
used by VGLM/VGAM modelling functions, e.g.,
\code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.
  

}

\note{

Further choices for the random noise, besides Gaussian, will be
implemented over time.


\code{zero} can be either an \emph{integer} vector
or a vector of \bold{character strings}
specifying either the position(s) or name(s) (partially or not) of the 
parameter(s) modeled as intercept-only. For \code{MAXff},
the parameters are placed and named as follows (by convention):


\code{c("Mean", "noiseVar" or "noiseSD", "MAcoeff")}.

Users can modify the \code{zero} argument accordingly.
For simplicity, the second choice recommended.
See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}} 
for further details on \code{zero}.

If \bold{no} constraints are entered in the fitting process,
(e.g., via \code{\link[VGAMextra:cm.ARMA]{cm.ARMA}})
this family function internally verifies by default
whether the estimated series
is invertible (since \code{noChecks = FALSE}).
To ignore this step, set \code{noChecks = TRUE}.
If the estimated MA process is non-invertible, MLE
coefficients will conform with the corresponding
invertible MA model.

Further details about these checks are shown within 
the \code{summary()} output.





}

\references{

Miranda, V. and Yee, T.W. (2018)
Vector Generalized Linear Time Series Models.
\emph{In preparation}.

Madsen, H. (2008)
\emph{Time Series Analysis}
Florida, USA: \emph{Chapman & Hall}.
(Sections 5.3 to 5.5).


Tsay, R. (2013) 
An Introduction to Analysis of Financial data with \code{R}. 
New Jersey, USA: \emph{Wiley} Sections 2.2 to 2.4.


}


\author{

 Victor Miranda and Thomas W. Yee.
 
}


\section{Warning}{

By default, a moving-average model of order-\eqn{1} is fitted.

If different, the \code{order} is recycled up to the number of responses
entered in the \code{vglm \ vgam} call has been matched.

Successful convergence depends on reasonably setting initial values. If 
initial values computed by the algorithm are not adequate, make use of 
the the optional initial values (\code{imean}, \code{isd}, etc.)

For constraints on the paramaters see
\code{\link[VGAMextra:cm.ARMA]{cm.ARMA}}.



}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{

  \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}},
  \code{\link[VGAMextra:ARXff]{ARXff}},
  \code{\link[VGAMextra:checkTS.VGAMextra]{checkTS.VGAMextra}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[VGAM:Links]{Links}},
  \code{\link[VGAM:vglm]{vglm}},
  \code{\link[VGAM:vglm]{vgam}}.
  
  
}

\examples{
set.seed(2)
nn    <- 130
### Coefficients
phi1  <-  0.34; phi2 <- -1.19; phi3 <- 0.26
### Intercept
mu    <- c(-1.4, 2.3)
### Noise standar deviations (Two responses)
sdMA  <- c(sqrt(6.5), sqrt(4.0))
### A single covariate.
Xcov <- runif(nn)

# Generating two MA processes, TS1 and TS2, Gaussian noise.
# Note, the SD noise for TS2 is function of Xcov.

y1   <- mu[1] + arima.sim(nn, 
                          model = list( ma = c(phi1, phi1^2)), 
                          rand.gen = rnorm, sd = exp(sdMA[1]) ) 
y2   <- mu[2] + arima.sim(nn, 
                          model = list( ma = c(phi1, phi2, phi3) ), 
                          rand.gen = rnorm, sd = exp(Xcov + sdMA[2]) )
# OUR DATA
tsdata <- data.frame(x2 = Xcov , TS1 = y1, TS2 = y2)

#------------------------------------------------------------------------#
# 1. A simple MA(3) to compare with 'arima()'.

myfit0 <- vglm(TS1 ~ 1,
               MAXff(order = 3, type.EIM = "exact",
                    var.arg = FALSE),
               #constraints = cm.ARMA(Model = ~ 1, 
               #                      lags.cm = 2, 
               #                      Resp = 1),
               data = tsdata, trace = TRUE) 

Coef(myfit0)[c(3, 4, 1)]
fitArima <- arima(tsdata$TS1, order = c(0, 0, 2)) 
coef(fitArima)

AIC(myfit0); BIC(myfit0)

# ------------------------------------------------------------------------#
# 2. Estimate an MA(3), intercept-only, using initial values.

myfit <- vglm(TS2 ~ 1,
              MAXff(order = 3, type.EIM = c("exact", "approx")[1],
                   # Optional initial values.
                    imean = 3,
                    iMAcoeff = c(0.3, -0.2, 0.25),
                   var.arg = TRUE),
              data = tsdata, trace = TRUE)

Coef(myfit)
summary(myfit)
constraints(myfit)


#----------------------------------------#
# Same model fitted using arima()
#----------------------------------------#

fitArima <- arima(tsdata$TS2, order = c(0, 0, 3)) 
coef(fitArima)


#------------------------------------------------------------------------#
# 3. An MAX(3) with one covariate, testing its effect over the
#    standard deviation of the Gaussian noise. Note the 'zero' argument.

myfit1 <- vglm(TS2 ~ x2,
               # Or Multiple responses! 
               # cbind(TS1, TS2) ~ 1,
               MAXff(order = 3, type.EIM = "exact", xLag = 1,
                    # Optional initial values:
                    # idev.mean = 1.4, 
                    # iMAcoeff = c(2.3, -1.2, 0.25), isd = 1.6,
                    
                    # NOTE THE ZERO ARGUMENT:
                    zero = c("Mean", "MAcoeff"),
                    
                    var.arg = TRUE),
               data = tsdata, trace = TRUE)

coef(myfit1, matrix = TRUE) 
summary(myfit1)
vcov(myfit1)

constraints(myfit1)

#------------------------------------------------------------------------#
# Model above CANNOT be fitted using arima()
#------------------------------------------------------------------------#

}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
