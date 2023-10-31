\name{benini1Qlink}
\alias{benini1Qlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Link functions for the quantiles of several 1--parameter continuous
  distributions

}
\description{
  Computes the \code{benini1Qlink} transformation, its inverse and the
  first two derivatives.
  %
}
\usage{
  benini1Qlink(theta, p = stop("Argument 'p' must be entered."),
               y0 = stop("Argument 'y0' must be specified."),
               bvalue = NULL, inverse = FALSE,
               deriv = 0, short = TRUE, tag = FALSE)
}
\arguments{
  \item{theta}{
  Numeric or character. See below for further details.
   
  }
  \item{p}{
  Numeric. A single value between 0.0 and 1.0.
  It is the \eqn{p}--quantile to be modeled by this link function. 
  
  }
  \item{y0}{
  Same as \code{\link[VGAM:benini1]{benini1}}.
  
  }
  \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
  
  }
}
\details{
  This is a link function to model any \eqn{p}--quantile of the 
  1--parameter Benini distribution. It is called the
  \code{benini1Qlink} transformation defined as
  %
  \deqn{ \log y_0 + \sqrt{\frac{ -\log (1 - p) }{s}}  }{
          log y0  + ( -log(1 - p)/s )^(1/2),
  }
  %
  where \eqn{y_0 > 0}{y0 > 0} is a scale parameter and \eqn{s} is 
  a positive shape parameter, as in \code{\link[VGAM:benini1]{benini1}}. 

  Numerical values of \eqn{s} or \eqn{p} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  
  In particular, arguments \code{inverse} and \code{deriv} are
  disregarded if \code{theta} is character.
  

}
\value{
  For \code{deriv = 0}, the \code{benini1Qlink} transformation of
  \code{theta}, when \code{inverse = FALSE}. If \code{inverse = TRUE}, then
  the inverse transformation given by 
  \code{-log(1 - p) / (theta - log y0)^2} is returned.
  
  For \code{deriv = 1}, this function returns the derivative 
  \eqn{d} \code{eta} / \eqn{d} \code{theta}, 
  if \code{inverse = FALSE}. Else, the reciprocal 
  \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of \code{theta}.
  
  If \code{deriv = 2}, then the second order derivatives in terms of
  \code{theta} are accordingly returned. 
 
  
}
\author{

  V. Miranda and Thomas W. Yee.

}
\section{Warning}{
  The horizontal straight line \eqn{\log y0}{log y0} is a lower asymptote
  for this link function as \eqn{\theta}{theta} increases to \eqn{\infty}.
  Thus, when \code{inverse = TRUE} and \code{deriv = 0} entries at 
  \code{theta} becoming \eqn{\eta} must be greater than 
  \eqn{\log y0}{log y0}. Else, \code{NaN} will be returned. 
  See examples 2 and 3 below.
  
   
}
\note{

  Numerical instability may occur for values \code{theta} too close
  to zero or lower than \eqn{\log y0}. 
  Use argument \code{bvalue} to replace them before computing the link.

}
\seealso{
     \code{\link[VGAM:benini1]{benini1}}, 
     \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
  ## E1. benini1Qlink() and its inverse ##
   p <- 0.50; y0 = 1.25         ## Modeling the median
   my.s <- seq(0, 5, by = 0.1)[-1]
    max(my.s - benini1Qlink(benini1Qlink(my.s, p = p, y0 = y0), 
                            p = p, y0 = y0, inverse  =TRUE))    ## Zero

  ## E2. Plot of the benini1Qlink() transformation and its inverse     ##
  ## Note, inverse = TRUE implies that argument 'theta' becomes 'eta'. ##
  ## which must be greater than log(y0). Else, value less than log(y0) ##
  ## are replaced by NaN.                                              ##
  \donttest{
   #--- THE LINK
   my.b <- seq(0, 5, by = 0.01)[-1]
   plot(benini1Qlink(theta = my.b, p = p, y0 = y0) ~ my.b,
        type = "l", col = "blue", lty = "dotted", lwd = 3,
        xlim = c(-0.1, 6), ylim = c(-0.1, 5), las = 1,
        main = c("Blue is benini1Qlink(), green is the inverse"),
        ylab = "eta = benini1Qlink", xlab = "theta")
   abline(h = 0, v = 0, lwd = 2)
    
   #--- THE INVERSE
   lines(my.b, benini1Qlink(theta = my.b, p = p, y0 = y0, inv = TRUE),
         col = "green", lwd = 2, lty = "dashed")
   #--- Tracing the identity function for double--check
   lines(my.b, my.b)     
   }
    
   ## E3. WARNING! The first two values are less than log(y0)  ##
    benini1Qlink(theta = c(0.10, 0.15, 0.25, 0.35) , p = p, y0 = y0, inverse  = TRUE)
    
}