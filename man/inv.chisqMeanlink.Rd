\name{inv.chisqMlink}
\alias{inv.chisqMlink}
\alias{inv.chisqMeanlink}
%\alias{inv.chisqMlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 
   Link functions for the mean of 1--parameter 
   continuous distributions: The inverse chi--squared distribution.
    
}
\description{
   %
   Computes the \code{inv.chisqMlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
        inv.chisqMlink(theta, bvalue = NULL, inverse = FALSE,
                        deriv = 0, short = TRUE, tag = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{theta}{
     Numeric or character. This is \eqn{\theta}{theta} by default but may 
     be \eqn{\eta}{\eta} depending on the other parameters.
     See below for further details.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{
  
  This link functions models the mean of the
  inverse chi--squared distribution, 
  \code{\link[VGAMextra:inv.chisqff]{inv.chisqff}}.
  
  It is defined as 
  \deqn{ \eta = -\log ( df - 2),}{ \eta = -log(df - 2)}
  where \eqn{df} denotes the (non--negative) degrees of freedom, as in 
   \code{\link[VGAMextra:inv.chisqff]{inv.chisqff}}.
   
  Notice, however, that \eqn{df > 2} is required for the mean of
  this distribution to be real. Consequently, the domain set for
  \code{df} for this link function is \eqn{(2, \infty)}.
  
  Numerical values of \eqn{df} out of range will 
  result in \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{inv.chisqMlink} transformation of
   \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then the inverse \code{exp(-theta) + 2}.
   
   For \code{deriv = 1},
   \eqn{d} \code{eta} / \eqn{d} \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then 
   \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of 
   \code{theta}.
   
   When \code{deriv = 2}, the second derivatives in
   terms of \code{theta} are returned.

}
\section{Note}{

  Numerical instability may occur for values \code{theta} too large or
  possibly, too close to 2. Use argument \code{bvalue} to replace them
  before computing the link.
  
  If \code{theta} is character, then arguments \code{inverse} and
  \code{deriv} are ignored. See \code{\link[VGAM:Links]{Links}}
  for further details.
  
  
}
\author{
  V. Miranda and Thomas W. Yee.

}
\seealso{
    \code{\link[VGAMextra:inv.chisqff]{inv.chisqff}}, 
     \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
 ##  E1. Modelling the mean of the exponential distribution  ##
    set.seed(17010502)
    dof <- 2.5 
    isq.data <- data.frame(x2 = runif(100, 0, 1))
    isq.data <- transform(isq.data, y = rinv.chisq(n = 100, df = dof + x2))
    
    \donttest{
    hist(isq.data$y)
    }
    
    fit.inv <- vglm(y ~ x2, family = inv.chisqff(link = "inv.chisqMlink"), 
                    data = isq.data, trace = TRUE )
    coef(fit.inv, matrix = TRUE)
    summary(fit.inv)

  ## E3. Special values in a matrix ##
    (theta <- matrix(c(Inf, -Inf, NA, NaN, 1 , 2, 3, 4), ncol = 4, nrow = 2))
    inv.chisqMlink(theta = theta)   ## NaNs for df = theta <= 2 
 
  ## E2. inv.chisqMlink() and its inverse ##
    theta <- 0.1 + 1:5  # dof = df
    my.diff <- theta - inv.chisqMlink(inv.chisqMlink(theta = theta), inverse  =TRUE)
    summary(my.diff)     # Zero

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
