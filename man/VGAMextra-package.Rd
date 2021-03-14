\name{VGAMextra-package}
\alias{VGAMextra-package}
\alias{VGAMextra}
\docType{package}

\title{
Additions and extensions of the VGAM package.
%Extending the functionalities of the VGAM package. 
}
\description{
\pkg{VGAMextra} supplies additional functions and methods
to the package \pkg{VGAM}, under three main topics:


  ** \emph{Time series modelling}. A novel class 
  of VGLMs to model univariate time series, called
  \emph{vector generalized linear time series models}
  (VGLTSMs). It is
  characterized by incorporating
  \emph{past information} in the VGLM/VGAM loglikelihood.
 
  ** \emph{1--parameter distribution mean modelling}.
  We return full circle by developing new link functions for the
  mean of 1--parameter distributions. VGAMs, VGLMs and GAMLSSs are
  restricted to location, scale and shape, however, the
  VGLM/VGAM framework has infrastructure to accommodate
  new links as a function of the parameters.
 
  ** \emph{Quantile modelling of 1--parameter distributions}.
  Similarly, we have implemented link functions to model
  the quantiles of several 1--parameter distributions, for
  \emph{quantile regression}.
  %The framework is still under development.

}
\details{
The inference infrastructure of \pkg{VGAMextra} relies on the
VGLM/VGAM framework. Particularly, estimation is carried out via
IRLS using Fisher scoring, whilst additive models and reduced
rank regression are also
accommodated by all \pkg{VGAMextra} family functions.


At present, this package
allows the extent of VGLMs/VGAMs to operate
popular time series models as special cases, as well as
cointegrated time series (bivariate case), and 
modelling choices for volatility models incorporating
explanatories in the variance equation.
The central family functions in this respect are
\code{\link[VGAMextra:ARXff]{ARXff}},
\code{\link[VGAMextra:MAXff]{MAXff}},
\code{\link[VGAMextra:ARMAX.GARCHff]{ARMAX.GARCHff}}, and
\code{\link[VGAMextra:VGLM.INGARCHff]{VGLM.INGARCHff}}.

Regarding modelling the mean/quantile-functions, \pkg{VGAMextra}
affords links for several 1-parameter distributions,
e.g., \code{\link[VGAMextra:expMlink]{expMlink}},
\code{\link[VGAMextra:benini1Qlink]{benini1Qlink}}, or
\code{\link[VGAMextra:inv.chisqMlink]{inv.chisqMlink}}.
Collectively, the quantile-links represent an alternative
to quantile regression by directly modelling the quantile
function for distributions
beyond the exponential family (See Example 3 below).


The VGLM/VGAM framework is very large and encompasses a wide range of
multivariate response types and models, e.g., it includes
univariate and multivariate distributions,
categorical data analysis, and extreme values.
See \code{\link[VGAM:VGAM-package]{VGAM-package}} for a broad
description and further information on \pkg{VGAM}.


}
\section{Future work}{


* Implement VGLM time series family functions to handle
  error correction models (ECMs) for cointegrated time series.
  Upgrade this framework beyond the bivariate case, e.g., the
  the Vector ECM (VECMs).
  
* Upgrade VGLMs time series family functions to handle
  multivariate time series, e.g., the VAR model (Coming shortly).

* Incorporate VGLM/VGAM-links to model the
  mean and quantile functions of distributions with > 1 parameters.
 
* Develop the class of \emph{multiple} reduced--rank 
VGLMs towards time series, to handle, e.g.,
vector error correction models (VECMs),
for multiple cointegrated time series.

}
\references{

Pfaff, B. (2011)
Analysis of Integrated and Cointegrated Time Series with R.
Seattle, Washington, USA: \emph{Springer}.


Chan, N., Li, D., Peng, L. and Zhang, R. (2013)
Tail index of an AR(1) model with ARCH(1) errors.
\emph{Econometric Theory}, 29(5):920-940.


Yee, T. W. (2015)
Vector Generalized Linear and Additive Models:
With an Implementation in R.
New York, USA: \emph{Springer}.

Miranda, V. and Yee, T. W.
\emph{Vector generalized linear time series models.}
In preparation.


Miranda, V. and Yee, T. W.
\emph{On mean modelling of 1-parameter distributions using vector
generalized linear models}. In preparation



Miranda, V. and Yee, T. W.
\emph{Two--Parameter Link Functions with Applications to
Negative Binomail, Weibull, and Quantile 
Regression}. To be submitted to the 
Scandinavian Journal of Statistics.


Miranda, V. and Yee, T.W.
\emph{New Link Functions for Distribution-Specific Quantile
Regression Based on Vector Generalized Linear and Additive
Models}.  
Journal of Probability and Statistics,
Volume 2019, Article ID 3493628.




%Miranda, V. and Yee, T.W. (2018)
%Vector generalized linear and additive models towards 
%quantile modelling: An alternative to quantile regression.
%(\emph{Manuscript in preparation}).


Yee, T. W. (2008)
The \code{VGAM} Package.
\emph{R News}, \bold{8}, 28--39.


}
\section{Warning}{
\pkg{VGAMextra} is revised, altered, and/or upgraded on a regular basis.
Hence, be aware that any feature, e.g., function names, arguments, or
methods, may be modified without prior notice.
Check the \code{NEWS} for the latest changes and additions across
the different versions.


}
\author{
Victor Miranda, \email{victor.miranda@aut.ac.nz}.

Maintainer: Victor Miranda, \email{victor.miranda@aut.ac.nz}.

Contributor: Thomas Yee, \email{t.yee@auckland.ac.nz}.

}
\seealso{
  \code{\link[VGAM:vglm]{vglm}},
  \code{\link[VGAM:vgam]{vgam}},
  \code{\link[VGAM:rrvglm]{rrvglm}},
  \code{\link[VGAM:Links]{Links}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \url{https://www.stat.auckland.ac.nz/~vmir178/}.

}
\examples{
##### EXAMPLE 1. An AR(1) model with ARCH(1) errors.
# Chan et.al. (2013) proposed a long and technical methodology to
# estimate the tail index of an AR(1) with ARCH(1) errors involving 
# its estimation by QMLE. I fit this model straightforwardly by MLE
# using the family function ARXff() for time series, and constraining
# the effect of Y^2_{t - 1} to the conditional variance using
# constraint matrices.


# Generate some data
set.seed(1)
nn       <- ceiling(runif(1, 150, 160))
my.rho   <- rhobitlink(-1.0, inverse = TRUE)  # -0.46212
my.mu    <- 0.0
my.omega <- 1
my.b     <- 0.5
tsdata   <- data.frame(x2 =  sort(runif(n = nn)))
tsdata   <- transform(tsdata, index = 1:nn, TS1   = runif(nn))

for (ii in 2:nn)
  tsdata$TS1[ii] <- my.mu + my.rho * tsdata$TS1[ii-1] +
  sqrt(my.omega + my.b * (tsdata$TS1[ii-1])^2) * rnorm(1)
  
# Remove the burn-in data:
nnr <- ceiling(nn/5)
tsdata <- tsdata[-(1:nnr), ]
tsdata["index"] <- 1:(nn - nnr)

\donttest{
# The timeplot.
with(tsdata, plot(ts(TS1), lty = "solid", col = "blue", xlab ="Time", ylab = "Series"))
abline(h = mean(tsdata[, "TS1"]), lty = "dotted")
}
# The constraint matrices, to inhibit the effect of Y^2_{t - 1}
# over sigma^2 only.
const.mat <- list('(Intercept)' = diag(3), 'TS1l1sq' = cbind(c(0, 1, 0)))

# Set up the data using function WN.lags() from VGAMextra to generate
# our 'explanatory'
tsdata <- transform(tsdata, TS1l1sq = WN.lags(y = cbind(tsdata[, "TS1"])^2, lags = 1))

# Fitting the model
fit.Chan.etal <- vglm(TS1 ~ TS1l1sq, ARXff(order = 1,   # AR order 
                                          zero = NULL, noChecks = FALSE,
                                          var.arg = TRUE, lvar = "identitylink"),
             crit = "loglikelihood", trace = TRUE,
             constraints = const.mat, data = tsdata)  ## Constraints...
summary(fit.Chan.etal, lrt0 = TRUE, score0 = TRUE, wald0 = TRUE)
constraints(fit.Chan.etal)


###### EXAMPLE 2. VGLMs handling cointegrated (bivariate) time series.
# In this example, vglm() accommodates an error correction model
# of order (2, 2) to fit two (non-stationary) cointegrated time series.

# Simulating some data.
set.seed(2017081901)
nn <- 280
rho <- 0.75
s2u <- exp(log(1.5))  # Gaussian noise1
s2w <- exp(0)         # Gaussian noise2
my.errors <- rbinorm(nn, mean1 = 0, mean2 = 0, var1 = s2u, var2 = s2w, cov12 = rho)
ut <- my.errors[, 1]
wt <- my.errors[, 2]
yt <- xt <- numeric(0)

xt[1] <- ut[1]     # Initial value: error.u[0]
yt[1] <- wt[1]     # Initial value: error.w[0]
beta <- c(0.0, 2.5, -0.32)  # Coefficients true values.

for (ii in 2:nn) {
  xt[ii] <-  xt[ii - 1] + ut[ii]
  yt[ii] <- beta[1] + beta[2] * xt[ii]  + beta[3] * yt[ii - 1] + wt[ii]
}

# Regression of yt on xt, save residuals. Compute Order--1 differences.
errors.coint <- residuals(lm(yt ~ xt)) # Residuals from the static regression yt ~ xt
difx1 <- diff(ts(xt), lag = 1, differences = 1)  # First difference for xt
dify1 <- diff(ts(yt), lag = 1, differences = 1)  # First difference for yt

# Set up the dataset (coint.data), including Order-2 lagged differences.
coint.data <- data.frame(embed(difx1, 3), embed(dify1, 3)) 
colnames(coint.data) <- c("difx1", "difxLag1", "difxLag2",
                          "dify1", "difyLag1", "difyLag2")

# Remove unutilized lagged errors accordingly. Here, use from t = 3.
errors.cointLag1 <- errors.coint[-c(1:2, nn)]
coint.data <- transform(coint.data, errors.cointLag1 = errors.cointLag1)

\donttest{
# Plotting the data
plot(ts(yt[-c(1:3, NULL)]), lty = 4, type = "l", col = "red",
     main = "", xlab = "Time", las = 1, ylim = c(-32, 20),
     ylab = expression("Series"~x[t]~"and"~y[t]))
lines(ts(xt[-c(1:3, NULL)]), lty = 4, type = "l", col = "blue")
legend("bottomleft", c(expression("Series"~x[t]),
                       expression("Series"~y[t])),
       col = c("red", "blue"), lty = c(4, 4))
}

# Fitting an error correction model (2, 2), aka ECM(2, 2)
fit.coint <- vglm(cbind(dify1, difx1) ~ errors.cointLag1 + difxLag1 + difyLag1 + 
                                           difxLag2 + difyLag2,
          binormal(zero = c("sd", "rho")),  # 'sigma', 'rho' are intercept--only.
          trace = FALSE, data = coint.data)
summary(fit.coint)
coef(fit.coint, matrix = TRUE)


##### EXAMPLE 3. Quantile Modelling (QM).
# Here, the quantile function of the Maxwell distribution is modelled
# for percentiles 25\%, 50\% and 75\%. The resulting quantile-curves
# are plotted. The rate parameter is determined by an artificial covariate.

set.seed(123)
# An artificial covariate.
maxdata <- data.frame(x2 = sort(runif(n <- nn <- 120)))
# The 'rate' function.
mymu <- function(x) exp(2 -  6 * sin(2 * x - 0.2) / (x + 0.5)^2)
# Set up the data.
maxdata <- transform(maxdata, y = rmaxwell(n, rate = mymu(x2)))

# 25\%, 50\% and 75\% quantiles are to be modelled.
mytau <- c(0.25, 0.50, 0.75)
mydof <- 4

### Using B-splines with 'mydof'-degrees of freedom on the predictors
fit.QM <- vglm(Q.reg(y, pvector = mytau) ~ bs(x2, df = mydof),
                 family = maxwell(link = maxwellQlink(p = mytau), zero = NULL),
                 data = maxdata, trace = TRUE)
            
summary(fit.QM, lscore0 = TRUE)
head(predictors(fit.QM))         # The 'fitted values'

\donttest{
## The 25\%, 50\%, and 75\% quantile curves.
mylwd <- 1.5
with(maxdata, plot(x2, y, col = "orange",
                 main = "Example 1; Quantile Modelling",
                 ylab = "y", pch = "o", cex = 0.75))
with(maxdata, matlines(x2, predict(fit.QM)[, 1], col = "blue",
                       lty = "dotted", lwd = mylwd))
with(maxdata, matlines(x2, predict(fit.QM)[, 2], col = "chocolate",
                       lty = "dotted", lwd = mylwd))
with(maxdata, matlines(x2, predict(fit.QM)[, 3], col = "brown",
                       lty = "dotted", lwd = mylwd))
legend("topleft", c("percentile25", "percentile50", "percentile75"),
       lty = rep("dotted", 3), lwd = rep(mylwd, 3))
}

### Double check: The data (in percentage) below the 25\%, 50\%, and 75\% curves
round(length(predict(fit.QM)[, 1][(maxdata$y
                 <= predict(fit.QM)[, 1] )]) /nn, 5) * 100  ## Should be 25\% approx
round(length(predict(fit.QM)[, 2][(maxdata$y
                 <= predict(fit.QM)[, 2] )]) /nn, 5) * 100  ## Should be 50\% approx
round(length(predict(fit.QM)[, 3][(maxdata$y
                 <= predict(fit.QM)[, 3] )]) /nn, 5) * 100  ## Should be 75\% approx
}

\keyword{ package }% use one of  RShowDoc("KEYWORDS")
\keyword{ time series }% __ONLY ONE__ keyword per line
\keyword{ mean modelling }
\keyword{ quantile modelling }