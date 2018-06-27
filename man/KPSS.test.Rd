\name{KPSS.test}
\alias{KPSS.test}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{KPSS tests for stationarity}
\description{
  
  The Kwiatkowski-Phillips-Schmidt-Shin (KPSS) 
  test for the null hypothesis that the series \eqn{x} is
  \code{level} or \code{trend} stationary

}
\usage{
           KPSS.test(x, type.H0 = c("level", "trend")[1],
                        trunc.l = c("short", "large")[1],
                        show.output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   Numeric. A univariate series.
  
  }
  \item{type.H0}{
  The null hypothesis to be tested: either \code{level} or
  \code{trend} stationarity.
  
  }
  \item{trunc.l}{
  The lag truncation parameter. See below for more details.
  
  
  }
  \item{show.output}{
  Logical. Should the results be displyed? Default is 
  \code{TRUE}.
  }

}
\details{
To test the null hypothesis that a univariate time series
is level--stationary or stationary around a deterministic
trend. The alternative states the existence of a unit root.

Under this methodology, the series, say
\eqn{ \{ y_t;~t = 1, \ldots, T\} }{ {y[t]; t = 1, \ldots, T}}
is assumed to be decomposed as
%
\deqn{y_t = \rho t + \xi_t + \varepsilon_t,}{
      y[t] = \rho t + \xi_t + e_t,
}
that is, as the sum of a deterministic trend, 
a random walk (\eqn{\xi_t}), and a stationary error
(\eqn{\varepsilon_t \sim N(0, \sigma^2_z)}{
e[t] ~ N(0, \sigma^2[z]))}. 
Hence, this test reduces to simply test the
hypothesis that \eqn{\{ \xi_t \} }{{\xi_t}} 
is stationary,
that is, \eqn{H_0: \sigma^2_z = 0.}{Ho: \sigma^2[z] = 0.}


The test statistic combines the one--sided Lagrange
multiplier (LM) statistic 
and the \emph{locally best invariant} (LBI) test
statistic (Nabeya and Tanaka, (1988)). 
Its asymptotic distribution is discussed
in Kwiatkowski et al. (1992), and depends on 
the `long--run' variance \eqn{\sigma^2}.
The test statistic is given by
\deqn{\eta = T^{-2} \sum_i S^2_i / \widehat{\sigma}^2= 
 T^{-2} \sum_i S^2_i / s^2(l).}{
\eta = \sum_i S^2_i / (T^2 * \sigmahat^2(l)) = 
     \sum_i S^2_i / (T^2 * s^2(l)). }
     

where \eqn{s^2(l)} is a consistent estimate of
\eqn{\sigma^2}{\sigmahat^2}, given by
\deqn{s^2(l) = (1/T)\sum_{t = 1}^T\varepsilon^2_t + 
   (2 / T) \sum_{s = 1}^l w(s, l) \sum_{t = s + 1}^T \varepsilon_t
   \varepsilon_{t - s}.}{
   s^2(l) = (1/T) \sum_t e^2[t] + (2 / T) \sum_{s = 1}^l
   w(s, l) \sum_{t = s + 1}^T e[t] e[t - s].
   }


Here, \eqn{w(s, l) = 1 - s/(l + 1)}, where
\code{l} is taken from \code{trunc.l}, the 
lag--truncation parameter. The choice \code{"short"} 
gives the smallest integer not less than
\eqn{3 \sqrt{T} / 11}{3 *  T^{1/2} / 11}, or else,
\eqn{9 \sqrt{T} / 11}{9 * T^{1/2} / 11}, if
\code{trunc.l = "large"}.


Note, here the errors, \eqn{\varepsilon_t}{e[t]}, are
estimated from the regression 
\eqn{x ~ 1} (level) or \eqn{x ~ 1 + t} (trend), depending upon 
the argument \code{type.H0}.

Unlike other software using linear interpolates, here
the p--values for both, trend and level stationarity,
are interpolated by cubic spline interpolations 
from the tail critical values given 
in Table 1 in Kwiatkowski et al. (1992).
The interpolation takes place on \eqn{\eta}. 

}
\note{

There is no standard methodology to select an 
appropriate value for \code{trunc.l}, however,
satisfactory results have been found for 
\code{trunc.l} proportional to \eqn{T^{1/2}}.
See Andrews, D.W.K. (1991) for a discussion on this.
Empirically, this parameter may be \emph{suggested}
by the problem in turn, and should be large enough to 
approximate the true dynamic behaviour of the series.

}
\value{
A list with the following:

1) Test statistic and P-value,

2) Critical values,

3) Residuals, \eqn{\varepsilon_t}{e[t]}.

}
\references{

Andrews, D.W.K. (1991)
Heteroskedasticiy and autocorrelation consistent
covariance matrix estimation.
\emph{Econometrica}, \pkg{59}, 817--858.


Kwiatkowski, D., Phillips, P.C.B., 
Schmidt, P., and Shin, Y. (1992)
Testing the null hypothesis of stationarity against
the alternative of a unit root.
\emph{Journal of Econometrics}, \pkg{54}, 159--178.


Nabeya, S. and Tanaka, K. (1988)
Asymptotic theory of a test for the constancy
regression coefficients against the random walk
alternative.
\emph{Annals of Statistics}, \pkg{16}, 218--235.


Phillips, P.C.B. and Perron, P. (1988)
Testing for a unit root in time series regression.
\emph{Biometrika}, \pkg{75}, 335--346.


Phillips, P.C.B. (1987)
Time series with unit roots.
\emph{Econometrica}, \pkg{55}, 277--301

}
\author{

 Victor Miranda.
 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{

  \code{\link[VGAMextra:checkTS.VGAMextra]{checkTS.VGAMextra}}.
  
}

\examples{
set.seed(2802)
test <- KPSS.test(rnorm(20), type.H0 = "trend")
class(test)

test$crit.value
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
