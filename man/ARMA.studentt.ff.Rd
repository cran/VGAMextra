\name{ARMA.studentt.ff}
\alias{ARMA.studentt.ff}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{ VGLTSMs Family Functions:
Generalized autoregressive moving average model
with Student-t errors}

\description{
  
  For an ARMA model, estimates a 3--parameter
  Student-\eqn{t} distribution characterizing the errors
  plus the ARMA coefficients by MLE usign Fisher 
  scoring. Central Student--t handled currently.
}
\usage{
            ARMA.studentt.ff(order = c(1, 0),
                             zero = c("scale", "df"),
                             cov.Reg = FALSE,
                             llocation = "identitylink",
                             lscale    = "loge",
                             ldf       = "loglog",
                             ilocation = NULL,
                             iscale = NULL,
                             idf = NULL)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{order}{
   Two--entries vector, non--negative. 
   The order $u$ and $v$ of the ARMA model.
   
   
  }
  \item{zero}{
  Same as \code{\link[VGAM:studentt3]{studentt3}}.
  
  
  }
  \item{cov.Reg}{
  Logical. If covariates are entered, Should these be
  included in the ARMA model as a \code{Reg}ressand? Default is FALSE,
  then only embedded in the linear predictors.
  
  
  }
  \item{llocation, lscale, ldf, ilocation,
  iscale, idf}{
  Same as \code{\link[VGAM:studentt3]{studentt3}}.
  
  }
  
}

\details{

The normality assumption for time series analysis is relaxed to handle
heavy--tailed data, giving place to the ARMA model with shift-scaled
Student-\eqn{t} errors, another subclass of VGLTSMs. 


For a univariate time series, say \eqn{y_t}{y[t]},
the model described by this VGLTSM family function is
\deqn{
 \theta(B)y_t = \phi(B) \varepsilon_t,
}{
y[t] = \beta[0] + \beta[1] y[t - 1] + \ldots + \beta[p] y_[t - p] + e[t]+
  \phi[1] e[t - 1] + \ldots \phi[q] e[t - q],
}

where \eqn{\varepsilon_t}{e[t]} are distributed as 
a shift-scaled Student--\eqn{t}
with \eqn{\nu} degrees of freedom, i.e.,
\eqn{\varepsilon_t \sim t(\nu_t, \mu_t, \sigma_t)}{e[t] ~ t(\nu[t], \mu[t], \sigma[t])}.
This family functions estimates the location (\eqn{\mu_t}{mu[t]}),
scale (\eqn{\sigma_t}{\sigma[t]}) and degrees of freedom (\eqn{\nu_t}{\nu[t]})
parameters, plus the ARMA coefficients by MLE.

Currently only centered Student--t distributions are handled.
Hence, the non--centrality parameter is set to zero.

The linear/additive predictors are
\eqn{\boldsymbol{\eta} = (\mu, \log \sigma, \log \log \nu)^T,}{
\eta = (\mu, \sigma, log log \nu)^T,}
where \eqn{\log \sigma}{log \sigma} and \eqn{\nu} are intercept--only
by default.



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

If \code{order = 0}, then \code{AR.studentt.ff}
fits a usual 3--parameter Student--\eqn{t}, as with
\code{\link[VGAM:studentt3]{studentt3}}.

If covariates are incorporated in the analysis,
these are embedded in the location--parameter model.
Modify this through \code{zero}.
See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}} 
for details on \code{zero}.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[VGAMextra:ARIMAXff]{ARIMAXff}},
  \code{\link[VGAM:studentt]{studentt}},
  \code{\link[VGAM:vglm]{vglm}}.


}

\examples{
### Estimate the parameters of the errors distribution for an
## AR(1) model. Sample size = 50

set.seed(20180218)
nn <- 250
y  <- numeric(nn)
ncp   <- 0           # Non--centrality parameter
nu    <- 3.5         # Degrees of freedom.
theta <- 0.45        # AR coefficient
res <- numeric(250)  # Vector of residuals.

y[1] <- rt(1, df = nu, ncp = ncp)
for (ii in 2:nn) {
  res[ii] <- rt(1, df = nu, ncp = ncp)
  y[ii] <- theta * y[ii - 1] + res[ii]
}
# Remove warm up values.
y <- y[-c(1:200)]
res <- res[-c(1:200)]

### Fitting an ARMA(1, 0) with Student-t errors.
AR.stut.er.fit <- vglm(y ~ 1, ARMA.studentt.ff(order = c(1, 0)),
                       data = data.frame(y = y), trace = TRUE)

summary(AR.stut.er.fit)
Coef(AR.stut.er.fit)

\donttest{
plot(ts(y), col = "red", lty = 1, ylim = c(-6, 6), main = "Plot of series Y with Student-t errors")
lines(ts(fitted.values(AR.stut.er.fit)), col = "blue", lty = 2)
abline( h = 0, lty = 2)
}



}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
