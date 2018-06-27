\name{invgammaDist}
\alias{invgammaDist}
\alias{dinvgamma}
\alias{pinvgamma}
\alias{qinvgamma}
\alias{rinvgamma}
\title{The Inverse Gamma Distribution}


% ##################################
% Edition and updates.
% File created on 20141030.

% Modified on:
% 20141117
% #################################


\description{
  Density, distribution function, quantile function and random numbers 
  generator for the Inverse Gamma Distribution.
}


\usage{
dinvgamma(x, scale = 1/rate, shape, rate = 1, log = FALSE)
pinvgamma(q, scale = 1/rate, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)
qinvgamma(p, scale = 1/rate, shape, rate = 1, lower.tail = TRUE, log.p = FALSE)
rinvgamma(n, scale = 1/rate, shape, rate = 1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same as \code{\link[stats:GammaDist]{GammaDist}}. 
  
  
  }
  \item{scale, shape}{ 
  Scale and shape parameters, same as 
  \code{\link[stats:GammaDist]{GammaDist}}.
  Both must be positive.
  
  
  }
  \item{rate}{ 
  Same as \code{\link[stats:GammaDist]{GammaDist}}.


  }
  \item{log, log.p, lower.tail}{
  Same as \code{\link[stats:GammaDist]{GammaDist}}.
  
  
  }
}
\details{

  The Inverse Gamma density with parameters scale = \eqn{b}{b} and 
  shape = \eqn{s}{s} is given by
  \deqn{f(y) = \frac{b^{s}}{\Gamma(s)} y^{-s-1} e^{-b/y},}{%
        f(y) = (b^s / \Gamma(s)) * y^(-s-1) * e^(-b/y),}
  for \eqn{y > 0}{y > 0}, \eqn{b > 0}{b > 0}, and \eqn{s > 0}{ s > 0}. 
  Here, \eqn{\Gamma(\cdot)}{gamma()} is the gamma function as in 
  \code{\link[base:Special]{gamma}} 
  
  The relation between the Gamma Distribution and the Inverse Gamma 
  Distribution is as follows:
  
  Let \eqn{X}{X} be a random variable distributed as Gamma 
  (\eqn{b, s}{b , s}), then \eqn{Y = 1 / X}{Y = 1/X} is distributed as 
  Inverse Gamma (\eqn{1/b, s}{1/b , s}). 
  It is worth noting that the math relation between the \emph{scale} paramaters
  of both, the Inverse Gamma and Gamma distributions, is inverse.
  
  Thus, algorithms of \emph{dinvgamma(), pinvgamma(), qinvgamma()} and 
  \emph{rinvgamma()} underlie on the algorithms 
  \code{\link[stats:GammaDist]{GammaDist}}.
  %implemented 
  %in \R to compute \eqn{\Gamma(x)}, for all real \code{x} except zero and 
  %negative integers, using the algorithm provided by Abramowitz and Stegun 
  %(p. 255; 1972).
  
  Let \eqn{Y}{Y} distributed as Inverse Gamma (\eqn{b, s}{b , s}). Then the 
  \eqn{k^{th}}{k-th} moment of \eqn{Y} exists for \eqn{-\infty < k < s} and is 
  given by
  
  \deqn{E[Y^k] = \frac{\Gamma(s - k)}{\Gamma(s)} b^k.}{%
        E[Y^k] = (b^k) * (\Gamma(s - k)/\Gamma(s)).}

  The mean (if \eqn{s > 1}{s > 1}) and variance (if \eqn{s > 2}{s >2}) are
    \deqn{E[Y] = \frac{b}{(s - 1)}; \ \ \ Var[Y] = \frac{b^2}{(s - 1)^2 
    \times (s - 2)}.}{%
    E[Y] = b/(s - 1); Var[Y] = b^2 / [((s - 1)^2) * (s - 2)].}    
  % \deqn{E[Y] = \frac{b}{(s - 1)},}{%
  %  E[Y] = b/(s - 1),}
  %  \deqn{Var[Y] = \frac{b^{2}}{((s - 1)^{2}) \times (s - 2)}}{
  %  Var[Y] = b^2 / [((s - 1)^2) * (s - 2)],}
  
}
\value{

 \code{dinvgamma()} returns the density, \code{pinvgamma()} gives the 
  distribution function, \code{qinvgamma()} gives the quantiles, and 
  \code{rinvgamma()} generates random deviates.
  
}

\references{

  Kleiber, C. and Kotz, S. (2003) 
  \emph{Statistical Size Distributions in Economics and Actuarial Sciences}.
  Wiley Series in Probability and Statistics. Hoboken, New Jersey, USA.  
  %Abramowitz, M. and Stegun, I. A. (1972) 
  %\emph{Handbook of Mathematical Functions with Formulas, Graphs, and 
  %Mathematical Tables}, ch.6, p.255.
  %Dover, New York, USA.
  
}

\author{V. Miranda and T. W. Yee.}

\note{
  %\code{[dpqr]-Inverse Gamma} functions display \eqn{invalid arguments} 
  %and \eqn{warnings} according to \code{[dpqr]-Gamma} functions in
  %\code{\link[stats:GammaDist]{GammaDist}}.
  
  Unlike the \code{\link[stats:GammaDist]{GammaDist}}, small values
  of \eqn{a}{shape} (plus modest \eqn{\mu}{mu}) or very large values of 
  \eqn{\mu}{mu} (plus moderate \eqn{a > 2}{shape > 2}), generate Inverse Gamma 
  values so near to zero. Thus, \code{rinvgamma} in 
  \code{\link[VGAMextra:invgammaDist]{invgammaDist}} may return either values 
  too close to zero or values represented as zero in computer arithmetic.
  
  In addition, function \code{dinvgamma} will return zero for 
  \eqn{x = 0}{x = 0}, which is the limit of the Inverse Gamma density when 
  \eqn{'x'}{'x'} tends to zero.

}
\section{Warning}{

  The order of the arguments \emph{scale} and \emph{shape} does not match 
  \code{\link[stats:GammaDist]{GammaDist}}.
  
}
\seealso{

    \code{\link[stats:GammaDist]{GammaDist}},
    \code{\link[base:Special]{gamma}}.
  
}
\examples{

  # Example 1.______________________
  n        <- 20
  scale    <- exp(2)
  shape    <- exp(1)
  data.1   <- runif(n, 0, 1)
  data.q   <- qinvgamma(-data.1, scale = scale, shape = shape, log.p = TRUE)  
  data.p   <- -log(pinvgamma(data.q, scale = scale, shape = shape)) 
  arg.max  <- max(abs(data.p - data.1))     # Should be zero
  
  \donttest{
  # Example 2.______________________
  scale <- exp(1.0)
  shape <- exp(1.2)
  xx    <- seq(0, 3.0, len = 201)
  yy    <- dinvgamma(xx, scale = scale, shape = shape)
  qtl   <- seq(0.1, 0.9, by = 0.1)
  d.qtl <- qinvgamma(qtl, scale = scale, shape = shape)
  plot(xx, yy, type = "l", col = "orange", 
       main = "Orange is density, blue is cumulative distribution function",
       sub  = "Brown dashed lines represent the 10th, ... 90th percentiles",
       las = 1, xlab = "x", ylab = "", ylim = c(0, 1))
  abline(h = 0, col= "navy", lty = 2)
  lines(xx, pinvgamma(xx, scale = scale, shape = shape), col = "blue")
  lines(d.qtl, dinvgamma(d.qtl, scale = scale, shape = shape), 
        type ="h", col = "brown", lty = 3)
        }
}
%\keyword{distribution}
