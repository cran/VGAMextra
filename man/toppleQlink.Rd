\name{toppleQlink}
\alias{toppleQlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Link functions for the quantiles of several 1--parameter continuous
  distributions

}
\description{
  Computes the \code{toppleQlink} transformation, its inverse and the
  first two derivatives.
  %
}
\usage{
  toppleQlink(theta, p = stop("Argument 'p' must be specified."),
              bvalue = NULL, inverse = FALSE,
              deriv = 0, short = TRUE, tag = FALSE)
}
\arguments{
  \item{theta}{
  Numeric or character. It is \eqn{\theta}{theta} by default although
  it could be \eqn{\eta}{\eta} depending upon other arguments. 
  See \code{\link[VGAM:Links]{Links}} for further details about this.
   
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
  This link function comforms with requirements of \pkg{VGAM} in order
  to be compatible within the VGLM/VGAM framework. That is, monotonic,
  onto, among other features.
  In this line, the \code{toppleQlink} transformation arises as the
  proper link to model any quantile of the Topp--Leone distribution
  \code{\link[VGAM:topple]{topple}}. It is defined as
  %
  \deqn{\frac{1 - \sqrt{1 - p^{1/s}}}{m.max}.}{
        [1 - sqrt(1 - p^(1/s))] / m.max.}
  %
  Here, \eqn{s} is a shape parameter lying in \eqn{(0, 1)} 
  as in \code{\link[VGAM:topple]{topple}},
  whereas \eqn{m.max} stands for the maximum in \eqn{(0, 1)} of
  \deqn{1 - \sqrt{1 - p^{1/s}}}{1 - sqrt(1 - p^(1/s))}
  as a function of \eqn{s}. Note, \eqn{p} is prespecified (fixed) between 
  \eqn{0} and \eqn{1}. 
  
  Numerical values of \eqn{s} or \eqn{p} out of range will 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}
  correspondingly.
  
  Arguments \code{inverse} and \code{deriv} will be ignored
  if \code{theta} is character.

}
\value{
  For \code{deriv = 0}, the \code{toppleQlink} transformation of
  \code{theta}, when \code{inverse = FALSE}. If \code{inverse = TRUE}, then
  the inverse transformation 
  \code{ log(p)/log(1 - (1 - theta * }\eqn{m.max}\code{)^2)} 
  is returned.
  
  For \code{deriv = 1}, this function returns 
  \eqn{d} \code{eta} / \eqn{d} \code{theta}, 
  if \code{inverse = FALSE}. If \code{inverse = TRUE}, then the reciprocal 
  \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of \code{theta}.
  
  If \code{deriv = 2}, then the second order derivatives as a function of
  \code{theta}. 
 
  
}
\section{Warning}{
 The expression \eqn{p^{1/s}}{p^(1/s)} tends rapidly to zero specially
 for values of \eqn{s} less than \eqn{0.005}. Therefore, in such cases
 numerical values represented as zero may be returned when computing
 this link function, regardless the value of argument \code{inverse}.

}
\author{

  V. Miranda and Thomas W. Yee.

}
\note{

  Numerical instability may occur for values \code{theta} too close
  to 0.0 or 1.0. Use argument \code{bvalue} to replace them before
  computing the link.

}
\seealso{
     \code{\link[VGAM:topple]{topple}}, 
     \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
  ## E1. toppleQLink() and its inverse. ##
  my.s <- ppoints(50); p <- 0.75
  max(my.s - toppleQlink(toppleQlink(my.s, p = p), p = p, inverse  =TRUE))   ## Zero
  toppleQlink(theta = c(-0.15, -0.10, 0.25, 0.35) , p = p, inverse  = FALSE) ## NaNs
  toppleQlink(theta = c(-0.15, -0.10, 0.25, 0.35) , p = p, inverse  = TRUE)  ## NaNs
   
   
  ## E2. toppleQLink() for different avlues of 'p'. ##
  \donttest{
  plot(toppleQlink(theta = my.s, p = 0.05) ~ my.s,
       type = "l", col = "blue", lty = "dotted", lwd = 3,
       xlim = c(-0.1, 1.65), ylim = c(-0.1, 1.1), las = 1,
       main = c("The toppleQlink() transformation"),
       ylab = "eta = toppleQlink", xlab = "theta")
  abline(h = 0, v = 0, lwd = 2)
  abline(h = 1, v = 1, lty = "dotted", col = "green")
  lines(toppleQlink(theta = my.s, p = 0.25) ~ my.s, lwd = 2, lty = "dashed", col = "gray")
  lines(toppleQlink(theta = my.s, p = 0.50) ~ my.s, lwd = 2, lty = "dashed", col = "brown")
  lines(toppleQlink(theta = my.s, p = 0.75) ~ my.s, lwd = 2, lty = "dashed", col = "orange")
  lines(toppleQlink(theta = my.s, p = 0.95) ~ my.s, lwd = 2, lty = "dashed", col = "gray50")
  legend(1.1, 1.0, c("p = 0.05", "p = 0.25", "p = 0.50", "p = 0.75", "p = 0.95"),
         lwd = rep(2, 5), lty = rep("dashed", 5), 
         col = c("blue", "gray", "brown", "orange", "gray50"))
  }
}