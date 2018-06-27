\name{rayleighQlink}
\alias{rayleighQlink}
\title{

  Link functions for the quantiles of several 1--parameter continuous
  distributions

}
\description{
  Computes the \code{rayleighQlink} transformation, its inverse and the
  first two derivatives.
  %
}
\usage{
  rayleighQlink(theta, p = stop("Argument 'p' must be specified."),
                bvalue = NULL, inverse = FALSE,
                deriv = 0, short = TRUE, tag = FALSE)
  %
}
\arguments{
  \item{theta}{
  Numeric or character. It is \eqn{\theta}{theta} by default, although 
  it may be \eqn{\eta}{eta}. See \code{\link[VGAM:Links]{Links}} for
  additional details about this.
   
  }
  \item{p}{
  Numeric. A single value between 0.0 and 1.0.
  It is the \eqn{p}--quantile to be modeled by this link function. 
  
  }
  \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
  
  }
}
\details{
  
  This link function directly models any \eqn{p}--quantile of the 
  Rayleigh distribution specified by the argument \code{p}.
  It is called the \code{rayleighQlink} transformation defined as
  %
  \deqn{b \sqrt{-2 \log(1 - p)},}{b * sqrt(-2 * log(1 - p)),}
  %
  where \eqn{b > 0} is a scale parameter as in 
  \code{\link[VGAM:rayleigh]{rayleigh}}. 

  Numerical values of \eqn{b} or \eqn{p} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  
  If \code{theta} is character, then arguments \code{inverse} and
  \code{deriv} are discarded.
  
}
\value{
  For \code{deriv = 0}, the \code{rayleighQlink} transformation of
  \code{theta}, when {inverse = FALSE}. If \code{inverse = TRUE}, then
  this function returns \code{theta / sqrt(-2 log(1 - p))}.

  For \code{deriv = 1}, then the function returns
  \eqn{d} \code{eta} / \eqn{d} \code{theta}, 
  if \code{inverse = FALSE}. If \code{inverse = TRUE}, then
  \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of \code{theta}.
  
  If {deriv = 2}, then the second order derivatives in terms of
  \code{theta}. 
 
  
}
\author{
  V. Miranda and Thomas W. Yee.


}
\note{

  Numerical instability may occur for values \code{theta} too close
  to zero. Use argument \code{bvalue} to replace them before
  computing the link.

}
\seealso{
     \code{\link[VGAM:rayleigh]{rayleigh}}, 
     \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
 ## E1. rayleighQlink() and its inverse ##
  p <- 0.50                           ## Modeling the median
  my.b <- seq(0, 5, by = 0.1)[-1]
  max(my.b - rayleighQlink(rayleighQlink(my.b, p = p), p = p, inverse  =TRUE)) ## Zero


 ## E2. Special values ##
  rayleighQlink(theta = c(Inf, -Inf, NA, NaN), p = p)
    
    
 ## E3. Use of argument 'bvalue' ##
  rayleighQlink(theta = seq(-0.2, 1.0, by = 0.1), p = p)  # WARNING: NaNs if theta <= 0
  rayleighQlink(theta = seq(-0.2, 1.0, by = 0.1), p = p, bvalue = .Machine$double.xmin)
    
}

