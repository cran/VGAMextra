\name{gamma1Qlink}
\alias{gamma1Qlink}
\title{

  Link functions for the quantiles of several 1--parameter continuous
  distributions

}
\description{
  Computes the \code{gamma1Qlink} transformation, its inverse and the
  first two derivatives.
  %
}
\usage{
  gamma1Qlink(theta, p = stop("Argument 'p' must be specified."),
               bvalue = NULL, inverse = FALSE,
               deriv = 0, short = TRUE, tag = FALSE)
}
\arguments{
  \item{theta}{
  Numeric or character. It is \eqn{\theta}{theta} by default although
  it could be \eqn{\eta}{eta} depending upon other arguments. 
  See \code{\link[VGAM:Links]{Links}} for further details about this.
   
  }
  \item{p}{
  A numeric vector of \eqn{p}--quantiles (numbers between 0 and 1) 
  to be modeled by this link function. 
  
  }
  \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
  
  }
}
\details{
  
  This link function has been specifically designed to model any 
  \eqn{p}--quantile of the 1--parameter gamma distribution, 
  \code{\link[VGAM:gamma1]{gamma1}}, in the VGLM/VGAM context. 
  It is defined as
  %
  \deqn{\eta = \log {\tt{qgamma}}({\tt{p}}, {\tt{shape =}} s),}{
        \eta = log qgamma(p = p, shape = s),}
  %
  where \eqn{s} is a positive shape parameter as in 
  \code{\link[VGAM:gamma1]{gamma1}}, whilst
  \eqn{{\tt{qgamma()}}}{qgamma()} is the quantile function 
  \code{\link[stats:qgamma]{qgamma}}.
  %Notice, this link function is \emph{shifted} to the right by
  %\eqn{k.min} units to ensure a monotonic behavior and also to be 
  %an onto function, conforming with requirements for link functions 
  %in \pkg{VGAM}.
  
  The inverse of the \code{gamma1Qlink} cannot be expressed in closed
  form. Instead, the inverse image, \eqn{s_{\eta}}{s[\eta]},
  of \eqn{\eta} is numerically approximated by 
  \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}}.
  
  Numerical values of \eqn{s} or \eqn{p} out of range will 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}
  correspondingly.
  
  Arguments \code{inverse} and \code{deriv} are dismissed
  if \code{theta} is character.

}
\value{
  For \code{deriv = 0}, the \code{gamma1Qlink} transformation of
  \code{theta}, when {inverse = FALSE}. If \code{inverse = TRUE}, then
  \code{theta} becomes \eqn{\eta}, and therefore,
  the approximate inverse image of \eqn{\eta} is returned.
  
  For \code{deriv = 1}, the partial derivative 
  \eqn{d} \code{eta} / \eqn{d} \code{theta} is returned, 
  if \code{inverse = FALSE}. If \code{inverse = TRUE}, then the reciprocal 
  \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of \code{theta}.
  
  If {deriv = 2}, then the second order derivatives as a function of
  \code{theta}. 
 
  
}
\author{

  V. Miranda and Thomas W. Yee.

}
\note{

  Numerical instability may occur for values \code{theta} too large, or
  too close to 0.0. Use argument \code{bvalue} to replace them before
  computing the link.

}
\seealso{
     \code{\link[VGAM:gamma1]{gamma1}}, 
     \code{\link[stats:qgamma]{qgamma}},
     \code{\link[VGAM:Links]{Links}}.
     
}
\examples{

  ## E1. gamma1QLink() and values causing NaNs or out of range  ##
  
  p <- 0.75                            # The third quartile is of interest.
  my.s <- seq(0, 5, by = 0.1)[-1]
  
  max(my.s - gamma1Qlink(gamma1Qlink(my.s, p = p), p = p, inverse  =TRUE)) ## Zero
  
  ## E2. Special values of theta ##
  gamma1Qlink(theta = c(-0.15, -0.10, 0, 1:10) , p = p, inverse  = FALSE)  ## NaNs
  gamma1Qlink(theta = c(-5, -3, 0, 1:10) , p = p, inverse  = TRUE)         ## Out of range
   
   
  ## E3. Plot of gamma1QLink() and its inverse. ##
  \donttest{
  
    # gamma1Qlink()
    plot(gamma1Qlink(theta = my.s, p = p) ~ my.s,
         type = "l", col = "blue", lty = "dotted", lwd = 3,
         xlim = c(-0.1, 5), ylim = c(-5, 15), las = 1,
         main = c("Blue is gamma1Qlink(), green is the inverse"),
         ylab = "gamma1Qlink transformation", xlab = "theta")
     abline(h = 0, v = 0, lwd = 2)

     # The inverse
     lines(my.s, gamma1Qlink(theta = my.s, p = p, inverse = TRUE),
           col = "green", lwd = 2, lty = "dashed")
           
      # The identity function, for double-checking.
     lines(my.s, my.s, lty = "dotted")
  }
}