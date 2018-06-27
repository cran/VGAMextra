\name{invweibullDist}
\alias{invweibullDist}
\alias{dinvweibull}
\alias{pinvweibull}
\alias{qinvweibull}
\alias{rinvweibull}
\title{The Inverse Weibull Distribution}

% ##################################
% Edition and updates.
% File created on 20141030.

% Modified on:
% 20141117
% #################################


\description{

  Density, distribution function, quantile function and random numbers 
  generator for the Inverse Weibull Distribution.

}
\usage{
    dinvweibull(x, scale = 1, shape, log = FALSE)
    pinvweibull(q, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
    qinvweibull(p, scale = 1, shape, lower.tail = TRUE, log.p = FALSE)
    rinvweibull(n, scale = 1, shape)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same as \code{\link[stats:Weibull]{Weibull}}. 
  
  
  }
  \item{scale, shape}{ 
  Scale and shape parameters, same as \code{\link[stats:Weibull]{Weibull}}.
  Both must be positive. 
  
  
  }
  \item{log, log.p, lower.tail}{
  Same as \code{\link[stats:Weibull]{Weibull}}. 
  
  
  }
}
\details{

  The Inverse Weibull density with parameters scale = \emph{b} and 
  shape = \eqn{s}{s}, is
  
  \deqn{f(y) = s b^s y^{-s-1} \exp{[-(y/b)^{-s}}],}{%
        f(y) = s * (b^s) * y^(-s-1) * e^[(-y/b)^(-s)],}
  for \eqn{y > 0}{y > 0}, \eqn{b > 0}{b > 0}, and \eqn{s > 0}{s > 0}. 
  
  The Weibull distribution and the Inverse Weibull distributions are related
  as follows:
  
  Let \eqn{X}{X} be a Weibull random variable with paramaters scale =\eqn{b}{b}
  and shape =\eqn{s}{s}. Then, the random variable \eqn{Y = 1/X}{Y = 1/X} has 
  the Inverse Weibull density with parameters scale = \eqn{1/b}{1/b} and 
  shape = \eqn{s}{s}. Thus, algorithms of \emph{[dpqr]-Inverse Weibull} 
  underlie on \code{\link[stats:Weibull]{Weibull}}.

  Let \eqn{Y}{Y} be a r.v. distributed as Inverse Weibull (\eqn{b, s}{b , s}).
  The \eqn{k^{th}}{k-th} moment exists for 
  \eqn{-\infty < k < s}{ -Inf < k < s} and is given by
  \deqn{E[Y^k] = b^{k} \ \Gamma(1 - k/s).}{%
        E[Y^k] = (b^k) * \Gamma(1 - k/s).}
   
   
  The mean (if \eqn{s > 1}{s > 1}) and variance (if \eqn{s > 2}{s >2}) are
  \deqn{E[Y] = b \ \Gamma(1 - 1/s); \ \ \ Var[Y] = b^{2} \ 
  [\Gamma(1 - 2/s) - (\Gamma(1 - 1/s))^2].}{%
        E[Y] = b * \Gamma(1 - 1/s),  Var[Y] = (b^2) * 
        [\Gamma(1 - 2/s) - (\Gamma(1 - 1/s))^2].}
  
  % \deqn{E[Y] = (b) \ast \Gamma(1 - 1/s),}{%
  %      E[Y] = b * \Gamma(1 - 1/s),}    
  % \deqn{Var[Y] = (b^{2}) [\Gamma(1 - 2/s) - (\Gamma(1 - 1/s))^2].}{%
  %      Var[Y] = (b^2) * [\Gamma(1 - 2/s) - (\Gamma(1 - 1/s))^2].}
  
  Here, \eqn{\Gamma(\cdot)}{gamma()} is the gamma function as in 
  \code{\link[base:Special]{gamma}}. 

}

\value{

  \code{dinvweibull()} returns the density, \code{pinvweibull()} computes the 
  distribution function, \code{qinvweibull()} gives the quantiles, and 
  \code{rinvweibull()} generates random numbers from the Inverse Weibull 
  distribution.


}

\references{

  Kleiber, C. and Kotz, S. (2003) 
  \emph{Statistical Size Distributions in Economics and Actuarial Sciences}.
  Wiley Series in Probability and Statistics. Hoboken, New Jersey, USA.
  
  Abramowitz, M. and Stegun, I. A. (1972) 
  \emph{Handbook of Mathematical Functions with Formulas, Graphs, and 
  Mathematical Tables}.
  ch.6, p.255. Dover, New York, USA.
  
}

\author{V. Miranda and T. W. Yee.}

\note{
  %\code{[dpqr]-Inverse Weibull} functions display \eqn{invalid arguments} 
  %and \eqn{warnings} according to functions 
  %\code{\link[stats:Weibull]{Weibull}}.
  
  Small values of \emph{scale} or \emph{shape} will provide Inverse Weibull 
  values too close to zero. Then, function \code{rinvweibull()} with such 
  characteristics will return either values too close to zero or values
  represented as zero in computer arithmetic.
  
  
  The Inverse Weibull distribution, which is that of \eqn{X}{X} where 
  \eqn{1/X}{1/X} has the Weibull density, is known as the log-Gompertz
  distribution. Thus, in order to emphazise the continuity concept of the 
  Inverse Weibull density, if \eqn{x = 0}{x = 0}, then \code{dinvweibull} 
  returns zero, which is the limit of such a density when \eqn{'x'}{'x'} 
  tends to zero.
  
  
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{

  The order of the arguments of \emph{[dpqr]-Inverse Weibull} does not match 
  those in \code{\link[stats:Weibull]{Weibull}}.

}

\seealso{

     \code{\link[stats:Weibull]{Weibull}},
     \code{\link[base:Special]{gamma}}.

}

\examples{
  #(1) ______________
  n        <- 20
  scale    <- exp(2)
  shape    <- exp(1)
  data.1   <- runif(n, 0, 1)
  data.q   <- qinvweibull(-data.1, scale = scale, shape = shape, log.p = TRUE)  
  data.p   <- -log(pinvweibull(data.q, scale = scale, shape = shape)) 
  arg.max  <- max(abs(data.p - data.1))     # Should be zero
  
\donttest{
  #(2)_________________
   scale  <- exp(1.0)
    shape <- exp(1.2)
    xx    <- seq(0, 10.0, len = 201)
    yy    <- dinvweibull(xx, scale = scale, shape = shape)
    qtl   <- seq(0.1, 0.9, by =0.1)
    d.qtl <- qinvweibull(qtl, scale = scale, shape = shape)
    plot(xx, yy, type = "l", col = "red", 
         main = "Red is density, blue is cumulative distribution function",
         sub  = "Brown dashed lines represent the 10th, ... 90th percentiles",
         las = 1, xlab = "x", ylab = "", ylim = c(0,1))
    abline(h = 0, col= "navy", lty = 2)
    lines(xx, pinvweibull(xx, scale = scale, shape = shape), col= "blue")
    lines(d.qtl, dinvweibull(d.qtl, scale = scale, shape = shape), 
          type ="h", col = "brown", lty = 3)
        }
}
%\keyword{distribution}
