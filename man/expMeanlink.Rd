\name{expMeanlink}
\alias{expMeanlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 
   Link functions for the mean of 1--parameter 
   continuous distributions: The exponential distribution.
    
}
\description{
   %
   Computes the \code{expMeanlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
         expMeanlink(theta, location = 0, bvalue = NULL, inverse = FALSE,
                     deriv = 0, short = TRUE, tag = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{theta}{
     Numeric or character. This is \eqn{\theta}{theta} although may 
     be \eqn{\eta}{\eta} depending on the other parameters.
     See below for further details.
    
    }
    \item{location}{
    This is a known location parameter.
    Same as \code{location} in \code{\link[VGAM:exponential]{exponential}}.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{
  
  This is a link function to model the mean of the exponential 
  distribution, \code{\link[VGAM:exponential]{exponential}}. 
  It is defined as
  %
  \deqn{\eta = \log (A + \lambda^{-1}), }{ \eta = log (A + rate^(-1)), }
  %
  where \eqn{\lambda > 0}{rate > 0} is a rate parameter and \eqn{A} is a 
  known location parameter, 
  same as \code{\link[VGAM:exponential]{exponential}}.
  
  Numerical values of \eqn{\lambda}{rate} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{expMeanlink} transformation of
   \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then the inverse \code{exp(theta - A)^(-1)}.
   
   For \code{deriv = 1},
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
 ##  E1. Modelling the mean of the exponential distribution  ##
    set.seed(17010402)
    nn <- 100
    edata <- data.frame(x2 = runif(nn) - 0.5, x3 = runif(nn) - 0.5)
    edata <- transform(edata, eta = 0.2 - 0.7 * x2 + 1.9 * x3)
    
            #----- The mean is a function of 'x2' and 'x3' ------#
    edata <- transform(edata, rate = expMeanlink(eta, inverse = TRUE))
    
    edata <- transform(edata, y = rexp(nn, rate = rate))
    with(edata, stem(y))
    with(edata, hist(y))

    exp.fit <- vglm(y ~ x2 + x3, exponential(link = "expMeanlink"), 
                    data = edata, zero = NULL, trace = TRUE, crit = "log")
   coef(exp.fit, matrix = TRUE)
   summary(exp.fit)
 
 ## E2. expMeanlink() and its inverse ##
    theta <- 0.1 + 1:5
    location <- 1.5
    my.diff <- theta - expMeanlink(expMeanlink(theta = theta, 
                  location = location), location = location, inverse  =TRUE)
    summary(my.diff)     # Zero


  ## E3. Special values in a matrix ##
    theta <- matrix(c(Inf, -Inf, NA, NaN, 1 , 2), ncol = 3, nrow = 2)
    expMeanlink(theta = theta, location = location)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
