\name{truncNormal}
\alias{truncnormalDist}
\alias{dtruncnorm}
\alias{ptruncnorm}
\alias{qtruncnorm}
\alias{rtruncnorm}
\title{The Truncated Normal Distribution}


% ##################################
% Edition and updates.
% File created on 20230823.

% Modified on:
% 20230823
% #################################


\description{
  Density, distribution function, quantile function and random numbers 
  generator for the truncated normal distribution
}


\usage{
dtruncnorm(x, mean = 0, sd = 1, min.support = -Inf, max.support = Inf, log = FALSE)
ptruncnorm(q, mean = 0, sd = 1, min.support = -Inf, max.support = Inf) 
qtruncnorm(p, mean = 0, sd = 1, min.support = -Inf, max.support = Inf, log.p = FALSE)
rtruncnorm(n, mean = 0, sd = 1, min.support = -Inf, max.support = Inf)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n, mean, sd}{
  Same as \code{\link[stats:Normal]{Normal}}. 
  
  
  }
  \item{min.support, max.support}{ 
  Lower and upper truncation limits.

  

  }
  \item{log, log.p}{
  Same as \code{\link[stats:Normal]{Normal}}.
  
  
  }
}
\details{

  Consider \eqn{X \sim}{\sim} N(\eqn{\mu}{\mu},
  \eqn{\sigma^2}{\sigma^2}), with \eqn{A < X < B}{A < X < B},
  i.e., \eqn{X} restricted to \eqn{(A , B)}{(A, B)}.
  We denote \eqn{A} = \code{min.support} and
  \eqn{B} = \code{max.support}.
   
   
   Then the conditional
  random variable
  \eqn{Y = X \cdot I_{(A , B)} }{Y = X \cdot I_{(A , B)}}
  has a truncated normal distribution. Its p.d.f. is given by

    
  \deqn{
  f(y; \mu, \sigma, A, B) = (1 / \sigma) \cdot
  \phi(y^*) / [  \Phi(B^*) - \Phi(A^*) ],
  }{
  f(y; \mu, \sigma, A, B) = (1 / \sigma) \cdot
  \phi(y^*) / [  \Phi(B^*) - \Phi(A^*) ],
  
  }

  where \eqn{y^* = (y - \mu)/ \sigma}{y^* = (y - \mu)/ \sigma}, 
  \eqn{A^* = (A - \mu) / \sigma}{A^* = (A - \mu) / \sigma}, 
  and
  \eqn{B^* = (B - \mu) / \sigma}{B^* = (B - \mu) / \sigma}.
  
   Its mean is 
  \deqn{\mu + \sigma \cdot [ \phi(A) - \phi(B) ] / 
  [\Phi(B) - \Phi(A)]. }{
  \mu + \sigma \cdot [ \phi(A) - \phi(B) ] / [\Phi(B) - \Phi(A)]. }
  
  Here, \eqn{\Phi}{\Phi} is the standard normal c.d.f and 
  \eqn{\phi}{\phi} is the standard normal p.d.f.

}
\value{

 \code{dtruncnorm()} returns the density, 
 \code{ptruncnorm()} gives the 
  distribution function, 
  \code{qtruncnorm()} gives the quantiles, and 
  \code{rtruncnorm()} generates random deviates.
  
  \code{dtruncnorm} is computed from the definition, as 
  in 'Details'. \code{[pqr]}truncnormal are computed based
  on their relationship to the normal distribution.
  
}

\references{

  Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995) 
  \emph{Continuous Univariate Distributions},
  Second Edition (Chapter 13). Wiley, New York.
  %Kleiber, C. and Kotz, S. (2003) 
  %\emph{Statistical Size Distributions in Economics and Actuarial Sciences}.
  %Wiley Series in Probability and Statistics. Hoboken, New Jersey, USA.  
  %Abramowitz, M. and Stegun, I. A. (1972) 
  %\emph{Handbook of Mathematical Functions with Formulas, Graphs, and 
  %Mathematical Tables}, ch.6, p.255.
  %Dover, New York, USA.
  
}

\author{Victor Miranda and Thomas W. Yee.}

\seealso{

    \code{\link[stats:Normal]{Normal}},
    \code{\link[VGAM:truncweibull]{truncweibull}}.
  
}
\examples{

###############
## Example 1 ##

mymu <- 2.1   # mu
mysd <- 1.0   # sigma
LL   <- -1.0  # Lower bound
UL   <- 3.0   # Upper bound

## Quantiles:
pp <- 1:10 / 10
(quants <- qtruncnorm(p = pp , min.support = LL, max.support = UL, 
                        mean = mymu,  sd = mysd))
sum(pp - ptruncnorm(quants, min.support = LL, max.support = UL,
                      mean = mymu, sd = mysd))     # Should be zero

###############
## Example 2 ##

## Parameters
set.seed(230723)
nn <- 3000
mymu <- 12.7    # mu
mysigma <- 3.5  # sigma
LL <- 6     # Lower bound
UL <- 17    # Upper bound

## Truncated-normal data
trunc_data <- rtruncnorm(nn, mymu, mysigma, LL, UL)

## non-truncated data - reference
nontrunc_data <- rnorm(nn, mymu, mysigma)


\dontrun{
## Densities
par(mfrow = c(1, 2))
plot(density(nontrunc_data), main = "Non-truncated ND", 
     col = "green", xlim = c(0, 25), ylim = c(0, 0.15))
abline(v = c(LL, UL), col = "black", lwd = 2, lty = 2)
plot(density(trunc_data), main = "Truncated ND", 
     col = "red", xlim = c(0, 25), ylim = c(0, 0.15))


## Histograms
plot.new()
par(mfrow = c(1, 2))
hist(nontrunc_data, main = "Non-truncated ND", col = "green", 
     xlim = c(0, 25), ylim = c(0, 0.15), freq = FALSE, breaks = 22,
     xlab = "mu = 12.7, sd = 3.5, LL = 6, UL = 17")
abline(v = c(LL, UL), col = "black", lwd = 4, lty = 2)
hist(trunc_data, main = "Truncated ND", col = "red",
     xlim = c(0, 25), ylim = c(0, 0.15), freq = FALSE,
     xlab = "mu = 12.7, sd = 3.5, LL = 6, UL = 17")

}


## Area under the estimated densities
# (a) truncated data
integrate(approxfun(density(trunc_data)), 
          lower = min(trunc_data) - 1, 
          upper = max(trunc_data) + 1)

# (b) non-truncated data
integrate(approxfun(density(nontrunc_data)), 
          lower = min(nontrunc_data), 
          upper = max(nontrunc_data))

  
}
%\keyword{distribution}
