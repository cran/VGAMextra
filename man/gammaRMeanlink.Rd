\name{gammaRMeanlink}
\alias{gammaRMeanlink}

\title{
 
   Link functions for the mean of 2--parameter 
   continuous distributions: The Gamma distribution.
    
}
\description{
   %
   Computes the \code{gammaRMeanlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
        gammaRMeanlink(theta, rate = NULL, wrt.param = NULL,
                       bvalue = NULL, inverse = FALSE,
                       deriv = 0, short = TRUE, tag = FALSE)
                       
}
\arguments{
    \item{theta}{
     Numeric or character. This is \eqn{\theta}{theta} although may 
     be \eqn{\eta}{\eta} depending on the other parameters.
     See below for further details.
    
    }
    \item{rate}{
    A rate parameter. Same as
    \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}}
    
    }
    \item{wrt.param}{
    Positive integer, either \eqn{1} or \eqn{2}. The partial derivatives
    are computed with respect to one of the two linear predictors 
    involved in this link. Further details listed below.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{
  
  This link function allows to model the mean of the 2--parameter
  gamma distribution.
  
  The \code{gammaRMeanlink} transformation, \emph{for given} \eqn{\beta},
  is defined as
  %
  \deqn{ \eta = \eta(\alpha; \beta) = \log \frac{\alpha}{\beta},}{
         \eta = \eta(\alpha; \beta) = log (\alpha/\beta),
         }
  %
  where \eqn{\alpha > 0} is a \emph{shape} parameter, and \eqn{\beta > 0}
  is a \emph{rate} parameter, as in
  \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}}.
  Here, the gamma distribution is \emph{shape--rate} parametrized, with
  mean \eqn{\mu = \alpha / \beta}.
  
  Note, this link is expressly a function of \eqn{\alpha}, i.e.
  \eqn{\theta}, therefore \eqn{\beta} must be passed down to the
  function at every call.
  
  Numerical values of \eqn{\alpha} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{gammaRMeanlink} transformation of
   \code{theta}, i.e. \eqn{\alpha}, when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then \eqn{\theta} becomes \eqn{\eta},
   and the inverse,
   \code{exp(theta)} * \eqn{\beta}, \emph{for given} \eqn{\beta}, is
   returned.
   
   For \code{deriv = 1}, \code{theta} becomes 
   \eqn{\theta = (\alpha, \beta)=}\eqn{(\theta1, \theta2)},
   \eqn{\eta = (\eta1, \eta2)},
   and then, the argument \code{wrt.param} must be
   considered:
   
   A) If \code{inverse = FALSE}, then 
   \eqn{d} \code{eta1} / \eqn{d} \code{theta1} when 
   \code{wrt.param = 1}, and
   \eqn{d} \code{eta1} / \eqn{d} \code{theta2} if
   \code{wrt.param = 2}.
   
   
   B) For \code{inverse = TRUE}, this function returns
   \eqn{d} \code{theta1} / \eqn{d} \code{eta1} and
   \eqn{d} \code{theta2} / \eqn{d} \code{eta1} conformably arranged
   in a matrix, if \code{wrt.param = 1},
   as a function of \eqn{\theta_i}{\thetai}, \eqn{i = 1, 2}.
   When \code{wrt.param = 2} only
   \eqn{d} \code{theta1} / \eqn{d} \code{eta2} is returned,
   since \eqn{d} \code{theta2} / \eqn{d} \code{eta2} is 
   conviniently computed via the second link function.
   See \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} for
   further details.
   should return \eqn{d} \code{theta} / \eqn{d} \code{eta2} which 
   is zero since 
   has no effect on the function since 
   %returns \eqn{d} \code{theta} / \eqn{d} \code{eta2},
   %as a function of \code{theta} for both cases.
   
   
   Similarly, when \code{deriv = 2}, the second derivatives in
   terms of \code{theta} are returned.
   
  By default, the linear predictors in
  \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} are
  \eqn{\eta_1 = \log \frac{\alpha}{\beta}}{\eta1 = log (\alpha / \beta)},
  and \eqn{\eta_2 = \log \beta}{\eta2 = log \beta}
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
    \code{\link[VGAMextra:gammaRMean]{gammaRMean}},
    \code{\link[VGAM:gammaR]{gammaR}},
    \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
 
    alpha <- 0.1 + 1:8  # theta1
    rate  <- exp(1)     # theta2
 
 ## E1. gammaRMeanlink() and its inverse ##
    eta1  <- gammaRMeanlink(theta = alpha, rate = rate, inverse = FALSE)
    my.diff <- alpha - gammaRMeanlink(theta = eta1, rate = rate, inverse = TRUE)
    
    summary(my.diff)     # Zero
    
  ## E2. Special values arranged in a matrix ##
    rate <- matrix(rate, ncol = 3, nrow = 3)  #Ensure equal dimensions. 
    alpha <- matrix(c(Inf, -Inf, NA, NaN, -1 , 1, 0, -2, 2), ncol = 3, nrow = 3)
    gammaRMeanlink(theta = alpha, rate = rate)   # NaNs produced.
    
     
  ## E3. Special values arranged in a array ##
    dim(rate) <- NULL
    rate <- array(rate, dim = c(2, 2, 2))  #Ensure equal dimensions. 
    alpha <- array(0.1 + 1:8, dim = c(2, 2, 2))
    gammaRMeanlink(theta = alpha, rate = rate)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
