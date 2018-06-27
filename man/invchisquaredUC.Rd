% 20170104
\name{inv.chisqDist}
\alias{inv.chisqDist}
\alias{dinv.chisq}
\alias{pinv.chisq}
\alias{qinv.chisq}
\alias{rinv.chisq}
%
\title{
          The Inverse Chi--squared Distribution
%
}
\description{
    Density, CDF, quantile function and random number generator for the
    Inverse Chi--squared distribution.
  
}
\usage{
    dinv.chisq(x, df, log = FALSE)
    pinv.chisq(q, df, lower.tail = TRUE, log.p = FALSE)
    qinv.chisq(p, df, lower.tail = TRUE, log.p = FALSE)
    rinv.chisq(n, df) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q, p, n}{
  Same as \code{\link[stats:Chisquare]{Chisquare}}. 
  
  
  }
  \item{df, lower.tail, log, log.p}{
  Same as \code{\link[stats:Chisquare]{Chisquare}}.
  
  
  }
}
\details{
  
  The inverse chi--squared distribution with non--negative
  \code{df} = \eqn{\nu} degrees
  of freedom implemented here has density
  %
  \deqn{f(x; \nu) = \frac{ 2^{-\nu / 2} x^{-\nu/2 - 1} 
               e^{-1 / (2x)} }{ \Gamma(\nu / 2) }, }{
        f(x; \nu) = 2^(-\nu / 2) x^(-\nu/2 - 1) e^(-1 / (2x)) /
                 \Gamma(\nu / 2),
               }
  %
  where \eqn{x > 0}{x > 0}, and
  \eqn{\Gamma} is the \code{\link[base:gamma]{gamma}} function.
  
  
  The mean is \eqn{1 / (\nu - 2)}, for \eqn{\nu > 2}, and the variance is
  given by \eqn{2 / [(\nu - 2)^2 (\nu - 4)]}, for \eqn{\nu > 4}.
  
  Also, as with \code{\link[stats:Chisquare]{Chisquare}}, the degrees
  of freedom can be non--integer.

}
\value{

 \code{dinv.chisq} returns the density, \code{pinv.chisq} returns the 
  distribution function, \code{qinv.chisq} gives the quantiles, and 
  \code{rinv.chisq} generates random numbers from this distribution.
  
}

\references{

  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995) 
  \emph{Continuous Univariate Distributions}.
  Chapters 18 (volume 1) and 29 (volume 2). Wiley, New York.  
  %Abramowitz, M. and Stegun, I. A. (1972) 
  %\emph{Handbook of Mathematical Functions with Formulas, Graphs, and 
  %Mathematical Tables}, ch.6, p.255.
  %Dover, New York, USA.
  
}

\author{V. Miranda}

\note{
  
  Yet to do: A \emph{non--central} parameter as an argument, 
  if amenable. 
  
  Two similar versions of the Inverse chi--squared distribution with
  \eqn{\nu} degrees of freedom may be found in the literature, as 
  follows:
  
  Let \eqn{Y \sim}{Y ~} chi--squared\eqn{(\nu)}, then
  
  I. \eqn{1 / Y \sim}{1 / Y ~} Inverse chi--squared\eqn{(\nu)}, and
  II. \eqn{ \nu / Y \sim}{\nu / Y ~} Inverse chi--squared\eqn{(\nu)}.
  
  Here, the former, which is the popular version, has been implemented.
  


}
\section{Source}{
      Specifically, it is the probability distribution of a random variable
  whose reciprocal follows a chi--squared distribution, i.e.,
  
  If \eqn{Y \sim}{Y ~} chi--squared\eqn{(\nu)}, then \eqn{1/Y \sim}{1/Y ~}
  Inverse chi--squared\eqn{(\nu)}.
  
  As a result, \code{dinv.chisq, pinv.chisq, qinv.chisq} and 
  \code{rinv.chisq} use the functions 
  \code{\link[stats:Chisquare]{Chisquare}} as a basis. Hence,
  invalid arguments will lead these functions to return
  \code{NA} or \code{NaN} accordingly.


}
\seealso{

    \code{\link[stats:Chisquare]{Chisquare}},
    \code{\link[base:gamma]{gamma}}.
  
}
\examples{

  ##  Example 1  ## 
  nn <- 50; df <- 1.4
  data.1   <- ppoints(nn)
  data.q   <- qinv.chisq(-data.1, df = df, log.p = TRUE)  
  data.p   <- -log(pinv.chisq(data.q, df = df)) 
  max(abs(data.p - data.1))     # Should be zero
  
  \donttest{
  ##  Example 2  ##

  xx    <- seq(0, 3.0, len = 301)
  yy    <- dinv.chisq(xx, df = df)
  qtl   <- seq(0.1, 0.9, by = 0.1)
  d.qtl <- qinv.chisq(qtl, df = df)
  plot(xx, yy, type = "l", col = "orange", 
       main = "Orange is density, blue is cumulative distribution function",
       sub  = "Brown dashed lines represent the 10th, ... 90th percentiles",
       las = 1, xlab = "x", ylab = "", ylim = c(0, 1))
  abline(h = 0, col= "navy", lty = 2)
  lines(xx, pinv.chisq(xx, df = df), col = "blue")
  lines(d.qtl, dinv.chisq(d.qtl, df = df), type ="h", col = "brown", lty = 3)
  }
}
%\keyword{distribution}
