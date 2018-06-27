\name{genbetaIIDist}
\alias{genbetaIIDist}
\alias{dgen.betaII}
\alias{pgen.betaII}
\alias{qgen.betaII}
\alias{rgen.betaII}
\title{The Generalized Beta Distribution of the Second King}
% ##################################
% Edition and updates.
% File created on 201500202.
% Modified on:
% 
% #################################

\description{
  
  Density, distribution function, inverse distribution (quantile function)
  and random generation for the Generalized Beta of the Second Kind (GB2).

}
\usage{
  dgen.betaII(x, scale = 1.0, shape1.a = 1.0, shape2.p = 1.0, shape3.q = 1.0, 
             log = FALSE)
  pgen.betaII(q, scale = 1.0, shape1.a = 1.0, shape2.p = 1.0, shape3.q = 1.0, 
             lower.tail = TRUE, log.p = FALSE)
  qgen.betaII(p, scale = 1.0, shape1.a = 1.0, shape2.p = 1.0, shape3.q = 1.0, 
             lower.tail = TRUE, log.p = FALSE)
  rgen.betaII(n, scale = 1.0, shape1.a = 1.0, shape2.p = 1.0, shape3.q = 1.0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
  Vector of quantiles.
  %Same as \code{\link[stats:Beta]{Beta}}
  
  
  }
  \item{p}{
  Vector of probabilities.
  
  
  }
  \item{n}{
  Number of observations. If \code{length(n) > 1}, its length is taken to be
  the numbre required.
  
  
  }
  \item{scale, shape1.a, shape2.p, shape3.q}{ 
  Strictly positive scale and shape parameters.
  
  
  }
  \item{log, log.p, lower.tail}{
  Same meaning as in \code{\link[stats:Beta]{Beta}}
  
  
  }
}
\details{
  The GB2 Distribution is defined by the probability density (pdf)
  
  \deqn{f(y) = \frac{a x^{ap - 1}}{b^{ap} B(p, q) [1 + (y/b)^{a}]^{p + q},}}{%
        f(y) = a * x^(ap - 1) / b^(ap) * B(p, q) * [1 + (x/b)^a]^(p + q),
        }
  for \eqn{y > 0}{y > 0}, and \eqn{b, a, p, q > 0}{b, a, p, q > 0}.
  Here, \eqn{B(p, q)}{B(p, q)} is the beta function as in 
  \code{\link[base:Special]{beta}}.
  
  The GB2 Distribution and the Beta Distribution 
  (see \code{\link[stats:Beta]{Beta}}) are linked, as follows:
  Let \eqn{X}{X} be a random variable with the Beta density and parameters 
  \eqn{p = shape_{1}}{p = shape1} and \eqn{q = shape_{2}}{q = shape2}. 
  Then, introducing additional \eqn{b = scale}{b = scale} and
  \eqn{a = shape}{a = shape} parameters, the variable
  \deqn{Y = \frac{(x/b)^{a}}{1 + (x/b)^{a}}}{%
        Y = (x/b)^a / [1 + (x/b)^a]
        } 
  has the GB2 Distribution, with parameters \eqn{b, a, p, q}{b, a, p, q}.
  
  The GB2 \eqn{k^{th}}{k^th} moment exists for \eqn{-ap < k < aq}{-ap < k < aq} 
  and is given by
  
  \deqn{E(Y^{k}) = \frac{b^{k} B(p + k/a, q - k/a)}{B(p, q)}}{%
        E(Y^k) = b^k * B(p+ k/a, q - k/a) / B(p, q)
        } 

  or, equivalently, 
  \deqn{E(Y^{k}) = \frac{b^{k} \Gamma(p + k/a) \Gamma(q - k/a)}
                   {\Gamma(p) \Gamma(q)}).}{%
        E(Y^k) = \Gamma(p + k/a) * \Gamma(q - k/a) / \Gamma(p) * \Gamma(q).} 
  Here, \eqn{\Gamma(\cdot)}{\Gamma(.)} is the gamma function as in
  \code{\link[base:Special]{gamma}}.

}

\value{

  \code{dgen.betaII()} returns the density (p.d.f), \code{pgen.betaII()} gives 
  the distribution function (p.d.f), \code{qgen.betaII()} gives the quantile
  function (Inverse Distribution function), and \code{rgen.betaII()} generates 
  random numbers from the GB2 distribution.
 
}
\references{

  Abramowitz, M. and Stegun, I. A. (1972) 
  \emph{Handbook of Mathematical Functions with Formulas, Graphs, and 
  Mathematical Tables}, 
  \bold{ch.6}, p.255. Dover, New York, USA.
  
  Kleiber, C. and Kotz, S. (2003) 
  \emph{Statistical Size Distributions in Economics and Actuarial Sciences}.
  Wiley Series in Probability and Statistics. 
  Hoboken, New Jersey, USA.
  
  McDonald, J. B. and Xu, Y. J. (1995)
  A generalization of the beta distribution with applications.
  \emph{Journal of Econometrics}, \bold{66}, p.133--152.
  
  McDonald, J. B. (1984)
  Some generalized functions for the size distribution of income.
  \emph{Econometrica}, \bold{52}, p.647--663.
  
  
}

\author{V. Miranda and T. W. Yee}

\note{

  Values of the \code{shape2.p} parameter moderately close to zero may imply 
  obtaning numerical values too close to zero or values represented as zero 
  in computer arithmetic from the function \code{rgen.betaII()}.
  
  Additionally, for specific values of the arguments \code{x, q, p} such as 
  \code{Inf, -Inf, NaN} and \code{NA}, the functions \code{qgen.betaII()}, 
  \code{pgen.betaII()} and \code{qgen.betaII()} will return the limit when 
  the argument tends to such value. 
  
  In particular, the quantile \code{qgen.betaII()} retunrs zero for negative
  values and \eqn{Inf}{Inf} for missed probabilities greater than 1.0.
  
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\section{Warning}{
%  The order of the arguments of \emph{[dpqr]-Inverse Weibull} does not match 
%  those in \code{\link[stats:Weibull]{Weibull}}.
%}
\seealso{

     \code{\link[stats:Beta]{Beta}},
     \code{\link[base:Special]{beta}}.

}

\examples{

  # Setting parameters to both examples below.
  b  <- exp(0.4)      # Scale parameter.
  a  <- exp(0.5)      # Shape1.a
  p  <- exp(0.3)      # Shape2.p
  q  <- exp(1.4)      # Shape3.q
  
  # (1) ______________  
  probs.y <- seq(0.0, 1.0, by = 0.01)
  data.1  <- qgen.betaII(p = probs.y, scale = b, shape1.a = a, 
                        shape2.p = p, shape3.q = q)
  max(abs(pgen.betaII(q = data.1, scale = b, shape1.a = a, 
                     shape2.p = p, shape3.q = q)) - probs.y) # Should be 0.
 \donttest{
  # (2)_________________
  xx    <- seq(0, 10.0, length = 200)
  yy    <- dgen.betaII(xx, scale = b, shape1.a = a, shape2.p = p, shape3.q = q)
  qtl   <- seq(0.1, 0.9, by = 0.1)
  d.qtl <- qgen.betaII(qtl, scale = b, shape1.a = a, shape2.p = p, shape3.q = q)
  plot(xx, yy, type = "l", col = "red", 
       main = "Red is the GB2 density, blue is the GB2 Distribution Function",
       sub  = "Brown dashed lines represent the 10th, ..., 90th percentiles",
       las = 1, xlab = "x", ylab = "", xlim = c(0, 3), ylim = c(0,1))
  abline(h = 0, col = "navy", lty = 2)
  abline(h = 1, col = "navy", lty = 2)
  lines(xx, pgen.betaII(xx, scale = b, shape1.a = a, 
                       shape2.p = b, shape3.q = q), col= "blue")
  lines(d.qtl, dgen.betaII(d.qtl, scale = b, shape1.a = a, 
                        shape2.p = p, shape3.q = q), 
                        type ="h", col = "brown", lty = 3)
 }
}
%\keyword{}
