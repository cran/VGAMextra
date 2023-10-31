\name{checkTS.VGAMextra}
\alias{checkTS.VGAMextra}
\alias{checkTS.ffs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polynomial roots based on \emph{transfer operators} in Vector 
Generalized Time Series Family Functions


}
\description{

  \code{checkTS.VGAMextra} computes the polynomial roots as per 
  \emph{transfer operator} in Vector 
  Generalized Time Series Family Functions in \pkg{VGAMextra}


}
\usage{
        checkTS.VGAMextra(thetaEst = NULL, 
                          tsclass  = c("AR", "MA"), 
                          chOrder  = 1,
                          NofS     = 1,
                          retmod   = TRUE,
                          pRoots   = TRUE)
                          
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{thetaEst}{
  A vector of coefficients. Its lenght must be \code{NofS * chOrder}.
  If \code{\link[VGAM:vglm]{vglm}} is called,
  then \code{thetaEst} contains the estimated coefficients of the model
  specified in \code{formula} .
  
  
  }
  \item{tsclass}{
  Character indicating the model class to be checked. Presently, 
  options \code{"AR"} and \code{"MA"} are handled.
  
  
  }
  \item{chOrder}{
  Positive integer. The order of polynomial associated to the
  underlying procees involved:
  either $p$ or $q$, which apply for 
  \code{"AR"} or \code{"MA"} rspectively.
  %MA(\eqn{q}), or ARMA(\eqn{p, q}). If \code{\link[VGAM:vglm]{vglm}} 
  %is called, it must match that order specified in the time series family 
  %function.
  
  
  }
  \item{NofS}{
  A positive integer denoting the number of Time Series to verify.
  In the \code{\link[VGAM:vglm]{vglm}} environment, \code{NofS} is 
  the number of responses given in \code{formula}.
  
  
  }
  \item{retmod}{
  Logical. If \code{TRUE} (default), the \emph{Module} of all roots as per 
  transfer operator in the process established in \code{tsclass} 
  is returned. Else, essentially the roots are returned.
  
  
  }
  \item{pRoots}{
  Logical. If \code{TRUE} (default), the roots computed from estimated 
  models are displayed along with the time series family function 
  execution.
  
  
  }
  

}
\details{

Stationarity and/or Invertibility of time series (TS) are usually verified 
via the roots of the polynomial derived from the \emph{transfer
operators}. 

In particular, \code{checkTS.VGAMextra} computes such roots via
the coefficients estimated by vector generalized TS family functions 
available in \pkg{VGAMextra} ( \code{\link[VGAMextra:ARXff]{ARXff}}, and 
\code{\link[VGAMextra:MAXff]{MAXff}}).
%The form of these polynomials depends on the TS considered. 

Specifically, \code{checkTS.VGAMextra} verifies whether the TS
analyzed via \code{\link[VGAM:vglm]{vglm}} is \emph{stationary} or 
\emph{invertible}, accordingly.

Note that an autoregressive process of order-\eqn{p} 
[AR(\eqn{p})]  with coefficients
\eqn{\theta_{1}, \ldots, \theta_{p}}{\theta[1], \ldots, \theta[p]} 
can be written in the form
\deqn{  \theta(B) Y_{t} = \varepsilon_{t}, }{% 
        \theta(B) Y[t] = e[t],
} 
%If there exists an inverse operator \eqn{\pi(B)}, such that
%\deqn{ \pi(B) = \psi^{-1}(B), } then the linear process (1) can be also 
%written in the form
where
\deqn{ \theta(B) = 1 - \sum_{k = 1}^{p} \theta_{k} B^{k} }{%
       \theta(B) = 1 - \sum \theta[k] B^[k], k = 1, \ldots p.

}
Here, \eqn{ \theta(B) } is referred to as 
the \emph{transfer operator} of the process, and 
\eqn{B^{k} Y_{t} = Y_{t - k},}{ B^[k] Y[t] = Y[t - k]},
for \eqn{k = 0, 1, \ldots,p}, is the lagged single-function.

In general, an autoregressive process of order-\eqn{p} is
\emph{stationary} if the roots of 
\deqn{ \theta(z) = 1 - \theta_{1} z - \ldots - \theta_{q} z^q  }{%
       \theta(z) = 1 - \theta[1] z - \ldots - \theta[q] z^q
}
lie \emph{outside} the unit circle, i.e. \eqn{|z| > 1}.

Similarly, a moving-average process of order-\code{q} can be formulated 
(without loss of generality \eqn{\mu = 0}{ \mu = 0}) 

\deqn{ Y_{t} = \psi(B) \varepsilon_{t}, }{ Y[t] =  \psi(B) e[t], }

where \eqn{ \psi(B) } is the \emph{transfer operator}, given by

\deqn{ \psi(B) = 1 + \sum_{k = 1}^{q} \psi_{k} B^{k}, }{%
       \psi(B) = 1 + \sum \psi[k] B^[k], k = 1, \ldots q,
}
Note that \eqn{ \psi_{0} = 1 }{ \psi[0] = 1 }, and 
\eqn{B^{k} \varepsilon_{t} = \varepsilon{t - k} }{ B^[k] e[t] = e[t - k] }.

Hence, a moving-average process of order-\eqn{q} [MA(\eqn{q})], 
generally given by (note \eqn{\mu = 0})
\deqn{Y_{t} = \phi_1 \varepsilon_{t - 1} + \ldots + 
              \phi_q \varepsilon_{t - q} + \varepsilon_{t},}{%
              Y[t] =  \phi[1] e[t - 1] + \ldots +
              \phi[q] e[t - q] + e[t],
}
is \emph{invertible} if all the roots of

\deqn{ \phi(B) = 1 + \phi_{1} B + \ldots + \phi_{q} B^q  }{%
       \phi(B) = 1 + \phi[1] B + \ldots + \phi[q] B^q
}
lie \emph{outside} the unit circle., i.e.m \eqn{|z| > 1}.

Parallel arguments can be stated for autoregressive moving
aberage processes (ARMA). See Box and Jenkins (1970) for 
further details.


}

\value{

A vector whose elements are the roots of polynomials \emph{inherited}
from \emph{transfer operators} according to the process analyzed. 

Alternatively, the modules of roots can by returned instead of merely roots
via the \code{retmod} argument. 
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Box, G.E.P. and Jenkins, G.M. (1970) 
\emph{Time Series Analysis: Forecasting and Control}. 
Holden-Day, San Francisco, USA.

Madsen, H (2007)
\emph{Time Series Analysis}.
Chapman and Hall/CRC, Boca Raton, Florida, USA.


}
\author{

Victor Miranda and T. W. Yee.
%%  ~~who you are~~


}

\note{

For TS family functions in the VGLM/VGAM framework, 
\code{checkTS.VGAMextra} is called at the final iteration of Fisher scoring.
It means that the MLE estimates are actually evaluated to 
verify whether the process is \emph{stationary} or \emph{invertible}.

If any root has module less than \eqn{1 + 1e-5}, a \code{warning} is displayed 
for informative purposes.

Argument \code{thetaEst} manages the parameters of the TS model 
in turn.


}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{

The argument \code{thetaEst} manages the coefficients of the TS
model(s) in turn. Then, it must be \code{NofS * chOrder} length, 
where \code{NofS} is the number of responses established in the 
\code{\link[VGAM:vglm]{vglm}} call.
Here the coefficients for each response must be sequentially groped.

A moving average process is always stationary (See Madsen (2007) for
further details). Consequently, the \code{\link[VGAMextra:MAXff]{MAXff}} 
in \pkg{VGAMextra} verifies (by default) only for invertibility.
To enable this option set \code{nowarning = FALSE} in the 
\code{\link[VGAMextra:MAXff]{MAXff}} 
call.

Similarly, \code{\link[VGAMextra:ARXff]{ARXff}} verifies whether the TS data 
fitted is stationary, whereas \code{ARMAXff()} verifies both properties.

}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
%\dontrun{

# A moving average process order-3 with coeffs --> c(2.4, -5.6, 0.83)
#-------------------------#
# This is NOT invertible !
#-------------------------#

MAcoeffs <- c(2.4, -5.6, 0.83)
checkTS.VGAMextra(thetaEst = MAcoeffs, 
                  tsclass = "MA", 
                  chOrder = 3,
                  retmod = FALSE)


# AR process order-3 with coeffs --> c( 0.45, 0.45^2, 0.45^3 )
#-------------------------#
# This is stationary !
#-------------------------#

ARcoeffs <- c( 0.45 , 0.45^2 , 0.45^3 )
checkTS.VGAMextra(thetaEst = ARcoeffs, 
                  tsclass = "AR", 
                  chOrder = 3,
                  retmod = TRUE,
                  pRoots = TRUE)  # DEFAULT for 'pRoots'
%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
