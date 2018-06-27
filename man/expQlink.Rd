\name{expQlink}
\alias{expQlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

   Link functions for the quantiles of several 1--parameter 
   continuous distributions. 
   
}
\description{
   %
   Computes the \code{expQlink} transformation, its inverse and
   the first two derivatives.
   %A number of transformations, their inverse and the first
   %two derivatives are incorporated here to model
   %any quantile of the exponential, Benini, Rayleigh, 
   %Maxwell, Topp--Leone and the 1--parameter gamma
   %distributions.
   %
}
\usage{
  expQlink(theta, p = stop("Argument 'p' must be entered."),
           bvalue = NULL, inverse = FALSE,
           deriv = 0, short = TRUE, tag = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{theta}{
     Numeric or character. This is \eqn{\theta}{theta} although may 
     be \eqn{\eta}{eta} depending on the other parameters.
     See below for further details.
    
    
    }
    \item{p}{
    Numeric. A prespecified number between 0 and 1.
    The particular \eqn{p}--quantile to be modelled.
    For example, \eqn{p = 0.5} means that the median is considered 
    by this link function.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{

  This is a link function to model any fixed quantile, say 
  \eqn{\xi_p}{\xi[p]}, of the exponential distribution.
  It is called the \code{expQlink} transformation and is defined as
  %
  \deqn{ \log (1 - p)^{-1/\lambda},}{ 
         -log(1 - p) / rate,}
  %
  where \eqn{\lambda}{rate} is positive as in 
  \code{\link[VGAM:exponential]{exponential}}.
  
  Numerical values of \eqn{\lambda}{rate} or \eqn{p} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.
  
}
\value{

   With \code{deriv = 0}, the \code{expQlink} transformation of
   \code{theta} for prespecified \eqn{p} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then the inverse \code{-log(1 - p)/theta}.
   
   For \code{deriv = 1}, this link function returns
   \eqn{d} \code{eta} / \eqn{d} \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then 
   \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of 
   \code{theta}.
   
   Similarly, when \code{deriv = 2}, the second derivatives in
   terms of \code{theta} are returned.

}
\section{Note}{

  Numerical instability may occur for values \code{theta} too close
  to zero. Use argument \code{bvalue} to replace them before
  computing the link.
  
  If \code{theta} is character, then arguments \code{inverse} and
  \code{deriv} are ignored. See \code{\link[VGAM:Links]{Links}}
  for further details.
  
  
}
\author{
  V. Miranda and Thomas W. Yee.

}
\seealso{
    \code{\link[VGAM:exponential]{exponential}}, 
     \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
 ## E1. expQlink() and its inverse ##
    p <- 0.25            # Modelling the first quartile.
    my.theta <- seq(0, 5, by = 0.1)[-1]
    my.diff <- my.theta - expQlink(expQlink(my.theta, p = p), p = p, inverse  =TRUE)
    summary(my.diff)     # Zero


  ## E2. Special values ##
    expQlink(theta = c(Inf, -Inf, NA, NaN), p = p)
    
\donttest{
  ##  E3. Plot of expQlink() for different quantiles  ##
    plot(expQlink(my.theta, p = p) ~ my.theta, 
         type = "l", lty = "dotted", col = "blue", lwd = 2,
         main = "expQink(p) transformation", xlab = "theta", ylab = "expQLink",
         xlim = c(-0.5, 5), ylim = c(-0.5, 5))
    abline(h = 0, v = 0, lwd = 2)
    lines(my.theta, expQlink(my.theta, p = 0.50),  col = "green", lty = "dotted", lwd = 2)
    lines(my.theta, expQlink(my.theta, p = 0.75),  col = "red", lty = "dotted", lwd = 2)
    legend(2, 4,  c("p = 0.25", "p = 0.50", "p = 0.75"), col = c("blue", "green", "red"),
           lwd = c(2, 2, 2), lty = c("dotted", "dotted", "dotted"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
