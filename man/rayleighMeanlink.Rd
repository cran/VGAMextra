\name{rayleighMlink}
\alias{rayleighMlink}
%\alias{rayleighMeanlink}
\alias{maxwellMlink}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 
   Link functions for the mean of 1--parameter 
   continuous distributions: The Rayleigh and 
   the Maxwell distributions.
   
    
}
\description{
   %
   The \code{rayleighMlink} and the \code{maxwellMlink} 
   transformations, their inverse and the first two derivatives.
   %
}
\usage{
   rayleighMlink(theta, bvalue = NULL, inverse = FALSE,
                 deriv = 0, short = TRUE, tag = FALSE) 
                    
   maxwellMlink(theta, bvalue = NULL, inverse = FALSE,
                deriv = 0, short = TRUE, tag = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{theta}{
     Numeric or character. It is \eqn{\theta}{theta} by default, 
     but it may be \eqn{\eta}{eta} depending upon other 
     parameters.
     See \code{\link[VGAM:Links]{Links}} for further details.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{
  
  \code{rayleighMlink} and \code{maxwellMlink} are link functions
  to model  the mean of the Rayleigh distirbution,
  (\code{\link[VGAM:rayleigh]{rayleigh}}),
  and the mean of the Maxwell distribution,
  (\code{\link[VGAM:maxwell]{maxwell}}), 
  respectively.
  
  Both links are somehow defined as the 
  \eqn{ \log {\tt{theta}} }{ log(theta)} 
  plus an \emph{offset}. Specifcally,
  \deqn{ {\tt{rayleighMlink}}(b) = \log ( b * \gamma(0.5) / sqrt{2} ),}{
              rayleighMlink(b) = log ( b * \gamma(0.5) / sqrt(2) ),
  }
  where \eqn{b > 0} is a scale parameter as in
  \code{\link[VGAM:rayleigh]{rayleigh}}; and
  \deqn{ {\tt{maxwellhMlink}}(b) = \log ( a^{-1/2} * sqrt{8 / \pi}  ).}{
              maxwellMlink(b) = log ( a^(-1/2) * sqrt(pi/8) ).
  }
  
  Here, \eqn{a} is positive as in \code{\link[VGAM:maxwell]{maxwell}}.
  
  Non--positive values of \eqn{a} and/or \eqn{b} will result in 
  \code{NaN}, whereas values too close to zero will return
  \code{Inf} or \code{-Inf}.


  
}
\value{
  
   For \code{deriv = 0}, the corresponding transformation of
   \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then \code{theta} becomes \eqn{\eta}, and
   the inverse transformations 
   
   I)  \code{exp(theta) * sqrt}(2) / \code{gamma(0.5)} 
   for \code{rayleighMlink}, and
   
   
   II)  \eqn{8}\code{ * exp(-2 * theta) } / \code{gamma(0.5)}\eqn{^2} for 
   \code{maxwellMlink},
   
   are returned.
   
   
   For \code{deriv = 1},
   \eqn{d} \code{eta} / \eqn{d} \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then 
   \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of 
   \code{theta}.
   
   When \code{deriv = 2}, the second derivatives in
   terms of \code{theta} are returned.

}
\section{Note}{
 
  Values of \eqn{a} or \eqn{b} out of range, e.g. when covariates
  involved, may cause numerical instability. 
  Use argument \code{bvalue} to replace them before computing any link.
  
  If \code{theta} is character, then arguments \code{inverse} and
  \code{deriv} are ignored. See \code{\link[VGAM:Links]{Links}}
  for further details.
  
  
}
\author{
  V. Miranda and Thomas W. Yee.

}
\seealso{
    \code{\link[VGAM:maxwell]{maxwell}}, 
    \code{\link[VGAM:rayleigh]{rayleigh}}
    \code{\link[VGAM:Links]{Links}}.
     
     
}
\examples{
 ##  The link and its inverse ##
    theta <- 0.1 + 1:10
    eta   <- maxwellMlink(maxwellMlink(theta = theta), inverse  =TRUE)
    summary(eta - theta)     # Zero
    
    eta   <- rayleighMlink(rayleighMlink(theta = theta), inverse  =TRUE)
    summary(eta - theta)     # Zero


 ## Modelling the mean of the Maxwell distribution  ##
    set.seed(17010401)
    
    rate <- maxwellMlink(theta = 2, inverse = TRUE)   # ~ 0.046
    mdata <- data.frame(y = rmaxwell(1000, rate = rate ))
    
    fit <- vglm(y ~ 1, maxwell(link = "maxwellMlink"), 
                data = mdata, trace = TRUE, crit = "coef")
    
    coef(fit, matrix = TRUE)
    Coef(fit)
 

}
