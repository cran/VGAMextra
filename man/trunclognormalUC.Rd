\name{truncLognormal}
\alias{trunclognormalDist}
\alias{dtrunclnorm}
\alias{ptrunclnorm}
\alias{qtrunclnorm}
\alias{rtrunclnorm}
\title{The Truncated Log-Normal Distribution}


% ##################################
% Edition and updates.
% File created on 20230829.

% Modified on:
% 20230829
% #################################


\description{
  Density, distribution function, quantile function and random 
  generation for the truncated log-normal distribution
}


\usage{
dtrunclnorm(x, meanlog = 0, sdlog = 1, min.support = 0, max.support = Inf, log = FALSE)
ptrunclnorm(q, meanlog = 0, sdlog = 1, min.support = 0, max.support = Inf) 
qtrunclnorm(p, meanlog = 0, sdlog = 1, min.support = 0, max.support = Inf, log.p = FALSE)
rtrunclnorm(n, meanlog = 0, sdlog = 1, min.support = 0, max.support = Inf)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, meanlog, sdlog}{
  Same as \code{\link[stats:Lognormal]{Lognormal}}. 
  
  
  }
  \item{min.support, max.support}{ 
  Lower and upper truncation limits.

  

  }
  \item{log, log.p}{
  Same as \code{\link[stats:Lognormal]{Lognormal}}.
  
  
  }
}
\details{
   Consider \eqn{Y \sim} 
   \code{Lognormal}\eqn{(\mu_Y, \sigma_Y  )} restricted
   to \eqn{(A, B )}, 
   that is, \eqn{0 < A  = \code{min.support} < X < B 
   = \code{max.support}}.
  The (conditional)
  random variable
  \eqn{Y = X \cdot I_{(A , B)} }{Y = X \cdot I_{(A , B)}}
  has a log--truncated normal distribution. Its p.d.f. is given by
  
  \deqn{
  f(y; \mu, \sigma, A, B) = (y^{-1} / \sigma) \cdot
  \phi(y^*) / [  \Phi(B^*) - \Phi(A^*) ],
  }{
  f(y; \mu, \sigma, A, B) = (1 / \sigma) \cdot
  \phi(y^*) / [  \Phi(B^*) - \Phi(A^*) ],
  
  }

  where \eqn{y^* = [\log(y) - \mu_Y]/ \sigma_Y}{y^* = 
  [log(y) - \mu_Y]/ \sigma_Y}, 
  \eqn{A^* = [\log(A) - \mu_Y] / \sigma_Y}{
  A^* = [log(A) - \mu_Y] / \sigma_Y}, 
  and
  \eqn{B^* = [\log(B) - \mu_Y] / \sigma_Y}{B^* = [\log(B) - 
  \mu_Y] / \sigma_Y}.
  
  Its mean is:
  \deqn{\exp(\mu + \sigma^2/2) \cdot \{\Phi[(\log(B) - \mu) / \sigma
  - \sigma] - \Phi[(\log(A) - \mu) / \sigma
  - \sigma] \} / \{ 
  \Phi[(\log(B) - \mu) / \sigma] - \Phi[(\log(A) - \mu) / \sigma] \}.
    }
  
  Here, \eqn{\Phi}{\Phi} is the standard normal c.d.f and 
  \eqn{\phi}{\phi} is the standard normal p.d.f.

}
\value{

 \code{dtrunclnorm()} returns the density, 
 \code{ptrunclnorm()} gives the 
  distribution function, 
  \code{qtrunclnorm()} gives the quantiles, and 
  \code{rtrunclnorm()} generates random deviates.
  
  %\code{dtruncnormal} is computed from the definition, as 
  %in 'Details'. \code{[pqr]}truncnormal are computed based
  %on their relationship to the normal distribution.
  
}

\references{

  Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) 
  \emph{Continuous Univariate Distributions}, Second Edition,
  (Chapter 13) Wiley, New York.
  %Kleiber, C. and Kotz, S. (2003) 
  %\emph{Statistical Size Distributions in Economics and Actuarial Sciences}.
  %Wiley Series in Probability and Statistics. Hoboken, New Jersey, USA.  
  %Abramowitz, M. and Stegun, I. A. (1972) 
  %\emph{Handbook of Mathematical Functions with Formulas, Graphs, and 
  %Mathematical Tables}, ch.6, p.255.
  %Dover, New York, USA.
  
}

\author{Victor Miranda and Thomas W. Yee.}

%\section{Warning}{
%
%  The of the arguments \emph{scale} and \emph{shape} does not match 
%  other lower/upper truncation limits, e.g.
%  \code{\link[VGAM:truncweibull]{truncweibull}}.
%  
%}
\seealso{

    \code{\link[stats:Lognormal]{Lognormal}},
    \code{\link[VGAMextra:truncnormal]{truncnormal}}.
  
}
\examples{

###############
## Example 1 ##

mymeanlog <- exp(0.5)    # meanlog
mysdlog   <- exp(-1.5)   # sdlog
LL   <- 3.5              # Lower bound
UL   <- 8.0              # Upper bound

## Quantiles:
pp <- 1:10 / 10
(quants <- qtrunclnorm(p = pp , min.support = LL, max.support = UL, 
                        mymeanlog, mysdlog))
sum(pp - ptrunclnorm(quants, min.support = LL, max.support = UL,
                      mymeanlog, mysdlog))     # Should be zero

###############
## Example 2 ##

set.seed(230723)
nn <- 3000

## Truncated log-normal data
trunc_data <- rtrunclnorm(nn, mymeanlog, mysdlog, LL, UL)

## non-truncated data - reference
nontrunc_data <- rtrunclnorm(nn, mymeanlog, mysdlog, 0, Inf)

\dontrun{
## Densities
plot.new()
par(mfrow = c(1, 2))
plot(density(nontrunc_data), main = "Non-truncated Log--normal", 
     col = "green", xlim = c(0, 15), ylim = c(0, 0.40))
abline(v = c(LL, UL), col = "black", lwd = 2, lty = 2)
plot(density(trunc_data), main = "Truncated Log--normal", 
     col = "red", xlim = c(0, 15), ylim = c(0, 0.40))


## Histograms
plot.new()
par(mfrow = c(1, 2))
hist(nontrunc_data, main = "Non-truncated Log--normal", col = "green", 
       xlim = c(0, 15), ylim = c(0, 0.40), freq = FALSE, breaks = 22,
       xlab = "mu = exp(0.5), sd = exp(-1.5), LL = 3.5, UL = 8")
abline(v = c(LL, UL), col = "black", lwd = 4, lty = 2)

hist(trunc_data, main = "Truncated Log--normal", col = "red",
     xlim = c(0, 15), ylim = c(0, 0.40), freq = FALSE, 
     xlab = "mu = exp(0.5), sd = exp(-1.5), LL = 3.5, UL = 8")
}


## Area under the estimated densities
# (a) truncated data
integrate(approxfun(density(trunc_data)), 
          lower = min(trunc_data) - 0.1, 
          upper = max(trunc_data) + 0.1)

# (b) non-truncated data
integrate(approxfun(density(nontrunc_data)), 
          lower = min(nontrunc_data), 
          upper = max(nontrunc_data))

  
}
%\keyword{distribution}
