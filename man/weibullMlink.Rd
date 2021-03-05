\name{weibullMlink}
\alias{weibullMlink}

\title{
 
   Link functions for the mean of 2--parameter 
   continuous distributions: The Weibull distribution.
    
}
\description{
   %
   Computes the \code{weibullMlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
       weibullMlink(theta, shape = NULL, wrt.param = NULL,
                    bvalue = NULL, inverse = FALSE,
                    deriv = 0, short = TRUE, tag = FALSE)
                       
}
\arguments{
    \item{theta}{
     Numeric or character. This is \eqn{\theta}{theta}
     ('scale' parameter) but it may 
     be \eqn{\eta}{\eta} depending on the other parameters.
     See below for further details.
    
    }
    \item{shape}{
    The shape parameter. Required for this link to work.
     See
    \code{\link[VGAMextra:weibullRff]{weibullRff}}.
    
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
  
   This is the link for the mean of 
   the 2--parameter Weibull distribution, also known as the
   \code{weibullMlink} transformation.
    It can only be used within
  \code{\link[VGAMextra:weibullRff]{weibullRff}}
  and is defined as
  %This is the link for the mean of the 2--parameter
  %Weibull distribution, as in 
  %\code{\link[VGAM:weibullR]{weibullR}}.
  %
  %
  \deqn{ \tt{weibullMlink}(\beta; \alpha) = \eta(\beta; \alpha) = 
                \log [\beta \cdot \Gamma (1 + 1/\alpha)],}{
         weibullMlink(\beta; \alpha) = 
         \eta(\beta; \alpha) = log (\beta * \Gamma(1 + 1/\alpha)),
         }
  %
 for \emph{given} \eqn{\alpha} ('shape' parameter) where
 \eqn{\beta > 0} is the \emph{scale} parameter.
%
 % a wrapper of 
  %
  \code{weibullMlink} is expressly a function of \eqn{\beta}, i.e.
  \eqn{\theta}, therefore \eqn{\alpha} (\emph{shape})
  must be entered at every call.
  
  Numerical values of \eqn{\alpha} or \eqn{\beta} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{weibmeanlnik} transformation of
   \code{theta}, i.e., \eqn{\beta}, when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then \eqn{\theta} becomes \eqn{\eta},
   and the inverse, 
   \eqn{\exp\left( theta - \log \Gamma(1 + 1/ \alpha)\right),}{
  % \code{exp(theta} - 
  % \eqn{\log \Gamma(1 + 1/\alpha)}\code{)}{
   exp[theta - log \Gamma(1 + 1 / \alpha)],
   }
   \emph{for given} \eqn{\alpha}, is
   returned.
   
   When \code{deriv = 1} \code{theta} becomes 
   \eqn{\theta = (\beta, \alpha)= (\theta_1, \theta_2)}, and
   \eqn{\eta = (\eta_1, \eta_2)} with 
   \eqn{\eta_2 = \log~\alpha}{\eta_2 = log~\alpha},
   and the argument \code{wrt.param} must be
   considered:
   
   A) If \code{inverse = FALSE}, then 
   \eqn{d} \code{eta1} / \eqn{d} \eqn{\beta} is returned when 
   \code{wrt.param = 1}, and
   \eqn{d} \code{eta1} / \eqn{d} \eqn{\alpha} if
   \code{wrt.param = 2}.
   
   
   B) For \code{inverse = TRUE}, this function returns
   \eqn{d} \eqn{\beta} / \eqn{d} \code{eta1} and
   \eqn{d} \eqn{\alpha} / \eqn{d} \code{eta1} conformably arranged
   in a matrix, if \code{wrt.param = 1},
   as a function of \eqn{\theta_i}{\thetai}, \eqn{i = 1, 2}.
   When \code{wrt.param = 2}, a matrix with columns
   \eqn{d\beta} / \eqn{d} \code{eta2} and 
   \eqn{d\alpha} / \eqn{d} \code{eta2}
   is returned.  
  % When \code{wrt.param = 2} only
  % \eqn{d} \code{theta1} / \eqn{d} \code{eta2} is returned, because
  %\eqn{d} \code{theta2} / \eqn{d} \code{eta2} can be
  % computed via the second link function.
  % See \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} for
  % further details.
  % should return \eqn{d} \code{theta} / \eqn{d} \code{eta2} which 
  % is zero since 
  % has no effect on the function since 
   %returns \eqn{d} \code{theta} / \eqn{d} \code{eta2},
   %as a function of \code{theta} for both cases.
   
   
   For \code{deriv = 2}, the second derivatives in
   terms of \code{theta} are likewise returned.
   
  %By default, the linear predictors in
  %\code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} are
  %\eqn{\eta_1 = \log \frac{\alpha}{\beta}}{\eta1 = log (\alpha / \beta)},
  %and \eqn{\eta_2 = \log \beta}{\eta2 = log \alpha}.
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
    \code{\link[VGAMextra:weibullQlink]{weibullQlink}},
    \code{\link[VGAMextra:weibullRff]{weibullRff}},
    \code{\link[VGAM:weibullR]{weibullR}},
    \code{\link[base:lgamma]{lgamma}},
    \code{\link[VGAM:Links]{Links}}.
    %\code{\link[VGAM:Links]{Links}}.
}
\examples{
    eta <- seq(-3, 3, by = 0.1) # this is eta = log(mu(b, a)).
    shape  <- exp(1)    # 'shape' argument.
 
 ## E1. Get 'scale' values with A WARNING (not the same length)!
   theta <- weibullMlink(theta = eta, shape = shape, inverse = TRUE)  # Scale
   
 \dontrun{
 ## E2. Plot theta vs. eta, 'shape' fixed.
 plot(theta, eta, type = "l", ylab = "", col = "blue",
      main = paste0("weibullMlink(theta; shape = ",
                    round(shape, 3), ")"))
  abline(h = -3:3, v = 0, col = "gray", lty = "dashed")
 }
 
 ## E3. weibullMlink() and its inverse ##
    etabis  <- weibullMlink(theta = theta, shape = shape, inverse = FALSE)
    summary(eta - etabis)     # Should be 0
}

 % plot(theta, eta, type = "l", las = 1, ylab = "", lty = 2, lwd = 4,
%   main = "weibullMlink(theta; shape)")
%   abline(v = 0, h = 0, col = "red")
   
 

% ## E3. Special values arranged in a array ##
%    rate <- beta; dim(rate) <- NULL
%    rate <- array(rate, dim = c(2, 2, 2))  #Ensure equal dimensions. 
%    alpha <- array(0.1 + 1:8, dim = c(2, 2, 2))
%    gammaRMeanlink(theta = rate, shape = shape)
  

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


    %    weibullRMeanlink(theta, shape = NULL, wrt.param = NULL,
    %                     bvalue = NULL, inverse = FALSE,
    %                     deriv = 0, short = TRUE, tag = FALSE)
    %   