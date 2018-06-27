\name{maxwellQlink}
\alias{maxwellQlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Link functions for the quantiles of several 1--parameter continuous
  distributions

}
\description{
  Computes the \code{maxwellQlink} transformation, its inverse and the
  first two derivatives.
  %
}
\usage{
  maxwellQlink(theta, p = stop("Argument 'p' must be specified."),
               bvalue = NULL, inverse = FALSE,
               deriv = 0, short = TRUE, tag = FALSE)
}
\arguments{
  \item{theta}{
  Numeric or character. See below for further details.
   
  }
  \item{p}{
  Numeric. A single value between 0 and 1.
  It is the \eqn{p}--quantile to be modeled by this link function. 
  
  }
  \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
  
  }
}
\details{
  
  In the VGLM/VGAM quantile regression context, this link function 
  can be used to model any \eqn{p}--quantile of the 
  Maxwell distribition. It is the \code{maxwellQlink} transformation
  given by 
  %
  \deqn{\sqrt{ 2~\textrm{\texttt{qgamma(p, 1.5) / a.}}}}{ 
         sqrt(2  * qgamma(p, 1.5) / a ).} 
  %
  Here, \eqn{a} is positive parameter as in 
  \code{\link[VGAM:maxwell]{maxwell}} 
  whereas  \code{\link[stats:qgamma]{qgamma}} is the quantile 
  function of the gamma distribution. 
  
  Numerical values of \eqn{a} or \eqn{p} out of range will 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}
  correspondingly.
  
  In particular, arguments \code{inverse} and \code{deriv} are
  disregarded if \code{theta} is character. 
  Also, if  \code{inverse = TRUE} and \code{deriv = 0}, then argument
  \code{theta} becomes \eqn{\eta}. See \code{\link[VGAM:Links]{Links}}
  for further details about this.
  

}
\value{
  For \code{deriv = 0}, the \code{maxwellQlink} transformation of
  \code{theta}, when {inverse = FALSE}. If \code{inverse = TRUE}, then
  the inverse given by \code{2*qgamma(p, 1.5) / theta^2 } is returned.
  
  For \code{deriv = 1}, this function returns the derivative 
  \eqn{d} \code{eta} / \eqn{d} \code{theta}, 
  if \code{inverse = FALSE}. Else, the reciprocal 
  \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of \code{theta}.
  
  If {deriv = 2}, then the second order derivatives in terms of
  \code{theta} are accordingly returned. 
 
  
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
     \code{\link[VGAM:maxwell]{maxwell}}, 
     \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
  ## E1. maxwellQlink() and its inverse ##
   p <- 0.25        ## Modeling the first quartile
   my.a <- seq(0, 5, by = 0.1)[-1]
   max(my.a - maxwellQlink(maxwellQlink(my.a, p = p), p = p, inverse  =TRUE)) ## Zero
    
  ## E2. The first two values are negative, NaN is returned ##
    maxwellQlink(theta = c(-0.15, -0.10, 0.25, 0.35) , p = p, inverse  = FALSE)
    maxwellQlink(theta = c(-0.15, -0.10, 0.25, 0.35) , p = p, inverse  = TRUE)

  ## E3. Plot of the maxwellQlink() and its inverse                    ##
  ## Note, inverse = TRUE implies that argument 'theta' becomes 'eta'. ##
   
  \donttest{
   #--- THE LINK
   
   plot(maxwellQlink(theta = my.a, p = p) ~ my.a,
        type = "l", col = "blue", lty = "dotted", lwd = 3,
        xlim = c(-0.1, 10), ylim = c(-0.1, 5), las = 1,
        main = c("Blue is maxwellQlink(), green is the inverse"),
        ylab = "eta = maxwellQlink", xlab = "theta")
   abline(h = 0, v = 0, lwd = 2)
    
   #--- THE INVERSE
   lines(my.a, maxwellQlink(theta = my.a, p = p,inv = TRUE), 
          col = "green", lwd = 2, lty = "dashed")
   lines(my.a, my.a)     # Tracing the identity function for double--check
   }
    
}