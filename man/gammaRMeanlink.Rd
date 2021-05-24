\name{gammaRMlink}
\alias{gammaRMlink}
\alias{gammaRMeanlink}

\title{
 
   Link functions for the mean of 2--parameter 
   continuous distributions: The gamma distribution.
    
}
\description{
   %
   Computes the \code{gammaRMlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
        gammaRMlink(theta, shape = NULL, wrt.param = NULL,
                    bvalue = NULL, inverse = FALSE,
                    deriv = 0, short = TRUE, tag = FALSE)
                       
}
\arguments{
    \item{theta}{
     Numeric or character. This is \eqn{\theta}{theta}
     ('rate' parameter) but iy may 
     be \eqn{\eta}{\eta} depending on the other parameters.
     See below for further details.
    
    }
    \item{shape}{
    The shape parameter. Same as
    \code{\link[VGAMextra:gammaRff]{gammaRff}}.
    
    }
    \item{wrt.param}{
    Positive integer, either \eqn{1} or \eqn{2}. The partial derivatives
    are computed with respect to one of the two linear predictors 
    involved with this link. Further details listed below.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{
  
  The link to model the mean of the 2--parameter
  gamma distribution.
  
  The \code{gammaRMlink} transformation, \emph{for given}
  \eqn{\alpha} ('shape' parameter),
  is defined as
  %
  \deqn{ \eta = \eta(\alpha; \beta) = \log \frac{\alpha}{\beta},}{
         \eta = \eta(\alpha; \beta) = log (\alpha/\beta),
         }
  %
  where \eqn{\beta > 0} is a \emph{rate} parameter. 
  %as in
  %\code{\link[VGAMextra:gammaRMeanff]{gammaRMeanff}},
  %a wrapper of 
  %\code{\link[VGAM:gammaR]{gammaR}}, adapted to handle this link.
  %The density 
%  is \emph{shape--rate} parametrized, with
%  mean \eqn{\mu = \alpha / \beta}.
  This link is expressly a function of \eqn{\beta}, i.e.
  \eqn{\theta}, therefore \eqn{\alpha} (\emph{shape})
  must be entered at every call.
  
  Numerical values of \eqn{\alpha} or \eqn{\beta} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{gammaRMlink} transformation of
   \code{theta}, i.e. \eqn{\beta}, when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then \eqn{\theta} becomes \eqn{\eta},
   and the inverse,
   \eqn{\alpha} * \code{exp(-theta)}, \emph{for given} \eqn{\alpha}, is
   returned.
   
   For \code{deriv = 1}, \code{theta} becomes 
   \eqn{\theta = (\beta, \alpha)=}\eqn{(\theta1, \theta2)}, and
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
   Also, when \code{wrt.param = 2}, a matrix with columns
   \eqn{d}\code{theta1} / \eqn{d} \code{eta2} and 
   \eqn{d}\code{theta2} / \eqn{d} \code{eta2}
   is returned.
  % See \code{\link[VGAMextra:gammaRMlink]{gammaRMlink}} for
  % further details.
  % should return \eqn{d} \code{theta} / \eqn{d} \code{eta2} which 
  % is zero since 
  % has no effect on the function since 
   %returns \eqn{d} \code{theta} / \eqn{d} \code{eta2},
   %as a function of \code{theta} for both cases.
   
   
   Similarly, when \code{deriv = 2}, the second derivatives in
   terms of \code{theta} are returned.
   
 % By default, the linear predictors in
%  \code{\link[VGAMextra:gammaRMlink]{gammaRMlink}} are
%  \eqn{\eta_1 = \log \frac{\alpha}{\beta}}{\eta1 = log (\alpha / \beta)},
 % and \eqn{\eta_2 = \log \beta}{\eta2 = log \alpha}.
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
    \code{\link[VGAMextra:gammaRff]{gammaRff}},
    \code{\link[VGAM:gammaR]{gammaR}},
    \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
 
    eta <- seq(-3, 3, by = 0.1) # this is eta = log(mu(b, a)).
    shape  <- rep(exp(0.8), length(eta))    # 'shape' argument.
 
 ## E1. Get 'rate' values.
   theta <- gammaRMlink(theta = eta, shape = shape, inverse = TRUE)  # rate
   
 \dontrun{
 ## E2. Plot theta vs. eta, 'shape' fixed.
   plot(theta, eta, type = "l", las = 1, ylab = "", 
   main = "gammaRMlink(theta; shape)")
 }
 
 ## E3. gammaRMlink() and its inverse ##
    etabis  <- gammaRMlink(theta = theta, shape = shape, inverse = FALSE)
    my.diff <- eta - etabis
    summary(my.diff)     # Zero
    
  ## E4. Special values arranged in a matrix ##
    bbeta <- matrix(eta[1:9], ncol = 3, nrow = 3)  #Ensure equal dimensions. 
    alpha <- matrix(c(Inf, -Inf, NA, NaN, -1 , 1, 0, -2, 2), ncol = 3, nrow = 3)
    # The gammaRMlink transformation (log(a/b))
    gammaRMlink(theta = bbeta, shape = alpha, inv = FALSE)   # NaNs produced.
    # Same as
    log(alpha/bbeta)
    
}
% ## E3. Special values arranged in a array ##
%    rate <- beta; dim(rate) <- NULL
%    rate <- array(rate, dim = c(2, 2, 2))  #Ensure equal dimensions. 
%    alpha <- array(0.1 + 1:8, dim = c(2, 2, 2))
%    gammaRMeanlink(theta = rate, shape = shape)
  

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
