\name{WN.InitARMA}
\alias{WN.InitARMA}
%\alias{InitialsMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimated White Noise (WN) from the autoregressive moving-average model
of order-(\eqn{p}, \eqn{q}) [ARMA(\eqn{p}, \eqn{q})].

}
\description{

 Estimates the unobserved white noise of the ARMA(\eqn{p}, \eqn{q}) 
 model via the corresponding inverted process.
 
 Also, provides the initial values of \code{\link[VGAMextra:ARXff]{ARXff}},
 \code{\link[VGAMextra:MAXff]{MAXff}}, and
 \code{ARMAXff} family functions.
 
 
}
\usage{
     WN.InitARMA(tsData    = NULL, 
                 order     = c(1, 0, 1),
                 whiteN    = FALSE, 
                 moreOrder = 0,
                 updateWN  = FALSE)
              
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tsData}{
  A univariate data frame containing the time series to be fitted according 
  to an ARMA(\eqn{p}, \eqn{q}) process. Data must be of class "\code{ts}".
  
  
  }
  \item{order}{
  A vector with three integer components. It is order of the ARMA model to be 
  inverted. These entries, \eqn{c(p, d, q)}, are the AR order, the 
  degree of differencing, and the MA order, respectively.
  %Specifically, if the initial values for the \code{MAXff} family function
  %are required, the \code{order} must match that order established
  %in the \code{MAXff()} call.
  
  
  }
  \item{whiteN}{
  Logical. If \code{TRUE}, then the estimated white noise computed from
  the inverted ARMA model is returned. This option is enabled only for 
  \code{\link[VGAMextra:MAXff]{MAXff}}, \code{ ARMAXff } family functions.
  
  
  }
  \item{moreOrder}{
  A non-negative integer (might be zero) used to increment the order of 
  the AR model initially fitted to estimate the residuals, i.e.,
  an AR(\eqn{p} + \code{moreOrder}) model.
  Empirically, values of \code{moreOrder} \eqn{ > 2} do NOT improve
  accuracy of estimates. This assert, however, may vary for different 
  time series family functions.
  %For more accurate results, \eqn{k = 2} is suggested.
  %Thus far, choices 1 and 2 are
  %available. If \code{moreOrder = k}, then an AR(k * order) is fitted, 
  %where \eqn{k = 1, 2}. 
  
  
  }
  \item{updateWN}{
  Logical. if \code{TRUE}, the white noise is \emph{updated} through a 
  second regression of \eqn{Y_{t}}{Y[t]} on 
  \eqn{Y_{t -1}, \ldots, Y_{t -p}, \widehat{\varepsilon_{t - 1}}, \ldots,
  \widehat{\varepsilon_{t - q} }}{Y[t - 1], \ldots, Y[t - p],
   ^e[t - 1], \ldots, ^e[t - q] }.
  
  
  }
}

\details{

Overall, the autoregressive moving average process of order \eqn{c(p, q)}, 
shortly denoted as ARMA(\eqn{p}, \eqn{q}), with \emph{intercept} 
\eqn{\mu} can be expressed as

\deqn{y_{t} = \mu + \theta_{1} y_{t - 1} + \ldots + \theta_{p} y_{t - p} +
              \phi_1 \varepsilon_{t - 1} + \ldots + 
              \phi_q \varepsilon_{t - q} + \varepsilon_{t}.}{%
              y[t] = \mu + y[t - 1] + \ldots + y[t - p] +
              \phi[1] e[t - 1] + \ldots + \phi[q] e[t - q] + e[t].
}

It is well known that it can be expressed in terms of an autoregressive 
process of infinite order, AR(\eqn{\infty}), by
recursive substitutions. For instance, given a mean-zero ARMA(1, 1),
\deqn{y_{t} = \theta_1 y_{t - 1} +  \phi_1 \varepsilon_{t - 1} + 
               \varepsilon_{t},
               \quad \quad (1) }{%
      Y[t] = \mu + \theta[1] y[t -1] + \phi[1] e[t - 1] + e[t] 
              \dots (1), }
one may express
\deqn{\varepsilon_{t - 1} = Y_{t - 1} - ( \theta_{1} y_{t - 2} + 
                          \phi_{1} \varepsilon_{t - 2} }{%
                          e[t - 1]  = Y[t - 1] - 
                          (\theta[1] y[t - 2] + \phi[1] e[t - 2] ). }

Substituting this equation in (1) yields the initial inverted process, 
as follows:
\deqn{ y_{t} = \psi_{1} y_{t - 1} + \psi_{2} y_{t - 2} +
                f(\varepsilon_{t - 2}, \varepsilon_{t} ). }{%
                y[t] = \psi[1] y_{t - 1} + \psi[2] y_{t - 2} +
                f(e[t - 2], e[t]).
                }
where \eqn{f} is a function of \eqn{\varepsilon_{t - 2}}{e[t - 2]} and
 \eqn{\varepsilon_{t}}{e[t]}.
 
Repeated substitutions as above produces the so-called \emph{inverted
process},
\deqn{ y_{t} = \sum_{k = 1}^{\infty} \psi_{k} y_{t - k} + 
               \varepsilon_{t}. \quad \quad (2) }{%
               y[t] = \sum \psi[k] y[t - k] + e[t],   (2)}

\eqn{ k = 1, \ldots, \infty}.
Hence, setting an acceptable order (via the \code{moreOrder} 
argument, \eqn{1} or \eqn{2} for instance), an 
AR(\eqn{p} + \code{moreOrd}) 
\emph{inverted} model is internally fitted 
within \code{WN.InitARMA}. Consequently, the unobserved white noise,
\eqn{ \{ \varepsilon_{t} \} }{ {e[t]} }, is estimated by computing 
the \emph{residuals} in (2), after regression. 
\code{whiteN = TRUE} enables this option.

Finally, initial values of the \code{\link[VGAMextra:MAXff]{MAXff}}, and 
\code{ARMAXff} family functions can be computed by least squares from 
the estimated white noise above, \eqn{ \{ \varepsilon_{t} \} }{ {e[t]} } 
and the given data, \eqn{ \{ t_{t} \} }{ {y[t]} }. 

Initial values of \code{\link[VGAMextra:ARXff]{ARXff}} 
are also internally computed using \eqn{ \{ t_{t} \} }{ {y[t]} } only.

}
\section{Warning}{

  For some time series family functions, 
  \code{\link[VGAMextra:MAXff]{MAXff}} for instance, values of 
  \code{moreOrder} \eqn{ > 3} do NOT improve
  the accuracy of estimates, and may lead the algorithm to failure to
  converge. 
  
  

}
\value{
A list with the following components:
\item{Coeff}{The initial values of the VGLM/VGAM family function in turn:
\code{\link[VGAMextra:ARXff]{ARXff}}, \code{\link[VGAMextra:MAXff]{MAXff}},
or \code{ARMAXff}.}
  
\item{whiteN}{(Optional) Estimated white noise enabled only for 
\code{\link[VGAMextra:MAXff]{MAXff}}, \code{ ARMAXff }.
That sequence is returned if \code{whiteN = TRUE}.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

}
\seealso{
  \code{\link[VGAMextra:MAXff]{MAXff}}, 
  \code{\link[VGAMextra:ARMAXff]{ARMAXff}}.
  
}
\references{

Brockwell, P. and Davis, R. (2002) 
\emph{Introduction to Time Series and Forecasting}. 
Springer, New York, USA.

Durbin, J. (1959) 
Efficient Estimation of Parameters in Moving-Average Models. 
\emph{Biometrika}, \bold{46}, pp 306--316.

}

\author{

 Victor Miranda and T. W. Yee.

}

\examples{

%\dontrun{ 
# Generating some data -> an MA(3) 
set.seed(1004)
mydata <- arima.sim( n = 200, list(ma = c(0.3, 0.56 , 0.11)) )

# Computing initial values to be passed to MAXff()
WN.InitARMA(tsData = data.frame(y = mydata), 
            order = c(0, 0, 3), 
            moreOrder = 1)

# Returning initial values and white noise.
initMA <- WN.InitARMA(tsData = data.frame(y = mydata), 
                      order = c(0, 0, 3), 
                      moreOrder = 1, 
                      whiteN = TRUE) 
                      
# Initial values passed to MAXff()
initMA$Coeff

# Estimated white noise
head(initMA$WhiteNoise)
                      
%}
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
