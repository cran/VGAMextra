\name{summaryS4VGAMextra-methods}
\docType{methods}
\alias{summaryS4VGAMextra}
\alias{summaryS4VGAMextra, vgltsmff-method}
%\alias{summaryS4VGAMextra,ANY,vgltsm-method}


\title{Summary methods for Vector Generalized Time Series Models
       
       }
       
\description{ S4 summary methods for models fitted with
              time series family functions from 
              \pkg{VGAMextra}. 
              
              These function are all 
              \code{\link[utils:methods]{methods}} for objects of class
              \code{\link[VGAM:vglm]{vglm}} with signature 
              \code{\link[VGAMextra:vgltsmff-class]{vgltsmff-class}}.
              

    }
       
\details{

Implementation of vector generalized time series (TS) family functions 
(\emph{vgltsmff}) in \pkg{VGAMextra} is entirely based on the structure
of family functions of class 
\code{\link[VGAM:vglmff-class]{vglmff-class}} 
from \pkg{VGAM}. More precisely, \emph{vgltsmff} family functions
can be created by calls of the form \code{new("vgltsmff",...)},
following the structure \code{\link[VGAM:vglmff-class]{vglmff-class}}.
See \code{\link[VGAM:vglmff-class]{vglmff-class}} for 
additional details.
%In this line, specific S4 dispatching methods for TS family functions
%at \pkg{VGAMextra} are called by those methods a

In this line, specific S4 dispatching methods are currently 
implemented at \pkg{VGAM} to show (or plot) essential 
statistical information about the model fitted.
%Additionally, S4 dispatching methods for TS family 
%functions (signature \code{vgtsff}) at 
%\pkg{VGAMextra} have been implemented to show 
%additional information commonly required by time 
%series practicioners. 

For the generic \code{summary}, specifically, S4 methods for 
objects with signature \code{vgtsff} are incorporated
in \pkg{VGAMextra} to display supplementary analyses commonly
required by TS practicioners. That is, additional information
to the default output shown by 
\code{\link[VGAM:summaryvglm]{summaryvglm}} for family functions
at \pkg{VGAM}, as follows: 


%Particularly, summary methods at \pgk {VGAMextra} for TS 
%family functions the base-generic \code{summary}, follow the same 
%structure shown by S4 methods for 
%\code{\link[VGAM:vglmff-class]{vglmff-class}}, as described at
%\code{\link[VGAM:summaryvglm]{summaryvglm}}. There are, however,
%some extra outputs which may be of relevance for time series
%practitioners. Here, two important analyses are additionally 
%incorporated, as follows:

a) The standard errors, which are computed from the 
asymptotic distribution of the MLE estimates, unlike the asymptotic 
approach (z-value) from \pkg{VGAM}.

b) Checks on stationarity and/or invertibility for autoregressive 
(AR), moving average (MA), and autoregressive moving-average (ARMA)
models via the polynomial roots.

c) The AIC, AICC and BIC criteria for model identification.

Notice that, for intercept-only models in the 'vglm' context, 
the asypmtotic distribution of the estimates, either conditional
or unconditional, will coincide with the theoretical distributions
as long as \eqn{n} increases. In particular, for the AR(\eqn{p})
process, the MLEs and the Yule-Walker estimates will concur 
asymptotically.

%Notice that, for the AR and MA processes, the Yule - Walker and
%maximum likelihood estimates (either conditional or unconditional)
%have the same asymptotic distribution. This remark, however, is valid
%for intercept-only models in the 'vglm' context. 
Where covariates or parameter constraints are involved, 
the standard errors for the estimates from time series 
family functions at \pkg{VGAMextra} are calculated from the 
predicted values allocated in the slot \code{@predictors},
when \code{summary(...)} is called.
%to objects of class \code{vglm} at the final iteration of IRLS, specifically, 
%in the slot \code{@predictors}. 
In this case, the \emph{conditional} mean, 
\eqn{\textrm{E}[\eta_j | \textbf{x}]}{
E[\eta(j) |} \bold{x} \eqn{]}, is considered as the estimate, where:
\deqn{\eta_j  = \sum_{k = 1}^{p} \beta_{(j)k} \times x_{k},}{
      \eta[ j ] = \sum \beta[ (j) k] * x[ k ],
}
for \eqn{j = 1, \ldots, M}.

}
\note{
As for the intercept, notice that this is called \emph{drift-term}
at \code{\link[VGAMextra:ARXff]{ARXff}} and
\code{\link[VGAMextra:ARMAXff]{ARMAXff}}, 
whilst it is refered as \emph{intercept} 
in \code{\link[VGAMextra:MAXff]{MAXff}}. 
This parameter is also estimated by TS family 
functions in \pkg{VGAMextra}. In the MA model, particularly, 
it is the mean of the process.

The drift-term, 
denoted as \eqn{\mu^*}{\mu^{*} },
is linearly linked to the mean,
\eqn{\mu}{\mu},
of the AR and ARMA processes in turn, as follows:
\deqn{ \mu \to \frac{\mu^{*} }{1 - \sum \theta_i}.  }{
\mu -> \mu^* /(1 - \sum \theta[i]).
}
Here, \eqn{\theta_i}{\theta[i]} are the AR coefficients.
Hence, the standard error for the \emph{drift-term} is 
accordingly computed based on the asymptotic distribution
of the mean. More precisely, the relation

\deqn{V(\mu^{*}) = (1 - \sum \theta_i)^{-2} \times 
                  \frac{\sigma_{\varepsilon}^2 }{n},}{
                  (1 - \sum \theta[i])^(-2) * \sigma[e]^2 / n,
                  } 
is considered, 
where \eqn{\sigma_{\varepsilon}^2}{\sigma[e]^2} is the variance
of the random errors.

Finally, the AIC, AICC and BIC criteria are computed of the well-known
expressions

\deqn{AIC = (-2) \times Log-likelihood + 2 \times k}{
      AIC = (-2)  * Log-likelihood + 2 * k,
}
\deqn{AICC = AIC + \frac{2~k~(k + 1)}{n - k - 1}}{
      AICC = AIC + 2 * k * (k + 1) /( n - k - 1),
}
and
\deqn{BIC = (-2) \times Log-likelihood + k~\times~ln(n)}{
      BIC = (-2)  * Log-likelihood + k * ln(n),
}

with \eqn{k} denoting the number of parameters.

}
\value{

An object of class 
\code{summary.vglm} printed by specific methods defined
at \pkg{VGAMextra} for objects with signature
vgltsff-class.

%\code{\link[VGAM:summary.vglm-class]{summary.vglm-class}}.

}
\seealso{

 \code{\link[VGAMextra:vgtsff-class]{vgtsff-class}},
 \code{\link[VGAM:summaryvglm]{summaryvlgm}},
 \code{\link[VGAMextra:ARXff]{ARXff}},
 \code{\link[VGAMextra:MAXff]{MAXff}},
 \code{\link[VGAMextra:ARMAXff]{ARMAXff}}, 
 \code{\link[VGAM:vglm]{vglm}}.
 

}
\references{

Woodward, H., Gray, H. and Elliot A. (2012)
\emph{Applied Time Series Analysis}.
Taylor & Francis/CRC, Florida, USA.


}
\author{
V. Miranda and T.W. Yee.

}
\examples{

\donttest{
#------------------------------------------------------------------------#
# Fitting a simple Moving Average model to compare with arima().
#------------------------------------------------------------------------#
set.seed(0628)
nn    <- 300
theta <- c(0.2, -0.37)  # Autoregressive coefficients
phi   <- c(0.25)        # MA coefficients.
mu    <- c(1.5, 0.85)   # Mean (not drift) of the MA process.
x2 <- runif(nn)

tsd1 <- mu[1]/(1 - sum(theta)) + 
                  arima.sim(n = nn, 
                            model = list(order = c(2, 0, 0), 
                                          ar = theta),
                            sd = exp(1.5))
tsd2 <- mu[2]/(1 - sum(theta)) + 
                  arima.sim(n = nn, 
                            model = list(order = c(2, 0, 1),
                                         ar = theta, ma = phi), 
                            sd = exp(1 + 2 * x2))

tsdata <- data.frame(TS1 = tsd1, TS2 = tsd2, x2 = x2)
head(tsdata)

    ###    An ARIMA(2, 0, 0) model, that is an AR(2) model    ###
    
#fit1 <- vglm(TS1 ~ 1, 
#             ARIMAXff(order = c(2, 0, 0), var.arg = FALSE, type.EIM = "exact"), 
#             data = tsdata,  crit = "log", trace = TRUE)

fit1 <- vglm(TS1 ~ 1, 
             ARXff(order = 2, var.arg = FALSE, type.EIM = "exact"), 
             data = tsdata,  crit = "log", trace = TRUE)
m.coe <- Coef(fit1)

## Using arima to compare to summary(vgtsff)
summary(fit1)
arima(tsdata$TS1, order = c(2, 0, 0)) ## Similar SE's than VGAMextra.


m.coe[1] / (1 - sum(m.coe[-(1:2)]))  # THIS IS SIMILAR TO THE INTERCEPT 
                                     # ESTIMATED BY arima(): 1.1898

    ###    An ARIMA(2, 0, 1) models, that is an ARMA(2, 1)     ###
    ###   The errors standard deviation is a function of 'x2'  ###

### NOTICE: ARIMA and ARMA use the "identitylink" for coefficients ###
#fit2 <- vglm(TS2 ~ x2, 
#             ARMAXff(order = c(2, 1), var.arg = FALSE, type.EIM = "exact",
#                     zero = NULL), 
#            # constraints = list('x2' = rbind(0, 1, 0, 0, 0)),
#             data = tsdata,  crit = "loglikelihood", trace = TRUE)

#m.coe <- coef(fit2)
#coef(fit2, matrix = TRUE)

## Compare summary(vglm) to arima().
#summary(fit2)
#arima(tsdata$TS2, order = c(2, 0, 1))

}

}


