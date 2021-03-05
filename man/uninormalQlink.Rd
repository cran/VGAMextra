\name{uninormalQlink}
\alias{uninormalQlink}

\title{
 
   Quantile regression:
   Link function for the quantiles of the normal distribution.
    
}
\description{
   %
   Computes the \code{uninormalQlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
     uninormalQlink(theta, percentile = stop("Enter percentiles."),
                    sd = NULL, wrt.param = NULL,
                    bvalue = NULL, inverse = FALSE,
                    deriv = 0, short = TRUE, tag = FALSE)
                       
}
\arguments{
    \item{theta}{
     Numeric or character. This is \eqn{\theta}{theta}
     ('mean' parameter) but iy may 
     be \eqn{\eta}{\eta} depending on the other parameters.
     See below for further details.
    
    }
    \item{percentile}{
    Numeric. A vector of percentiles of interest, denoted as \code{perc}.
    
    }
    \item{sd}{
    Numeric, positive. The 'standard deviation' parameter (required),
    denoted as \eqn{\sigma}. 
    
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
   A 2-parameter link for the quantiles of the normal 
   distribution. It can only be used within 
   \code{\link[VGAMextra:uninormalff]{uninormalff}} as the first linear
   predictor. It is defined as
   %This is the link for the mean of the 2--parameter
  %Weibull distribution, as in 
  %\code{\link[VGAM:weibullR]{weibullR}}.
  %
  %
  \deqn{ \tt{uninormalQlink}(\mu; \sigma) = \eta_1(\mu; \sigma) = 
                \mu + \sigma \cdot \Phi^{-1}(perc), }{
               % \log \{\beta \cdot [(-\log(1 - perc))^{(1/\alpha)}]\},}{
         uninormalQlink(\mu; \sigma) = 
         \eta_1(\mu; \sigma) = 
         \mu + \sigma * \Phi^(-1)(perc) }
  %
where \eqn{\Phi} is the error function
(see, e.g., \code{\link[VGAM:erf]{erf}}), 
\eqn{\mu in (-\infty, \infty)}, and \eqn{\sigma > 0}.
This link is expressly a function of \eqn{\theta = \mu}, therefore
\eqn{sigma} must be entered at every call.
%
 % a wrapper of 
  Numerical values of \eqn{\sigma} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{uninormalQlink} transformation of
   \code{theta}, i.e. \eqn{\mu}, when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then \eqn{\theta} becomes \eqn{\eta},
   and the inverse, \eqn{\eta - \sigma \Phi^{-1}(perc)},
   \emph{for given} \eqn{\sigma}, is returned.
   
   When \code{deriv = 1} \code{theta} becomes 
   \eqn{\theta = (\mu, \sigma)= (\theta_1, \theta_2)}, and
   \eqn{\eta = (\eta_1, \eta_2)} with 
   \eqn{\eta_2 = \log~\sigma}{\eta_2 = log \sigma},
   and the argument \code{wrt.param} must be
   considered:
   
   A) If \code{inverse = FALSE}, then 
   \eqn{d} \code{eta1} / \eqn{d} \eqn{\mu} is returned when 
   \code{wrt.param = 1}, and
   \eqn{d} \code{eta1} / \eqn{d} \eqn{\sigma} if
   \code{wrt.param = 2}.
   
   
   B) For \code{inverse = TRUE}, this link returns
   \eqn{d} \eqn{\mu} / \eqn{d} \code{eta1} and
   \eqn{d} \eqn{\sigma} / \eqn{d} \code{eta1} conformably arranged
   in a matrix, if \code{wrt.param = 1},
   as a function of \eqn{\theta_i}{\thetai}, \eqn{i = 1, 2}.
   When \code{wrt.param = 2}, then
   \eqn{d\mu} / \eqn{d} \code{eta2} and 
   \eqn{d\sigma} / \eqn{d} \code{eta2}
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
   terms of \code{theta} are similarly returned.
   
  %By default, the linear predictors in
  %\code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} are
  %\eqn{\eta_1 = \log \frac{\alpha}{\beta}}{\eta1 = log (\alpha / \beta)},
  %and \eqn{\eta_2 = \log \beta}{\eta2 = log \alpha}.
}
\section{Note}{
  
  Numerical instability may occur for values of
  \code{sigma} too close
  to zero. Use argument \code{bvalue} to replace the former only before
  computing the link.
  
  If \code{theta} is character, then arguments \code{inverse} and
  \code{deriv} are ignored. See \code{\link[VGAM:Links]{Links}}
  for further details.
  
  
}
\author{
  V. Miranda and Thomas W. Yee.

}
\seealso{
    \code{\link[VGAMextra:uninormalff]{uninormalff}},
    \code{\link[VGAM:uninormal]{uninormal}},
  %  \code{\link[VGAM:weibullR]{weibmeanlink}},
    \code{\link[VGAM:Links]{Links}}.
}
\examples{
    eta <- seq(-3, 3, by = 0.1) # this is eta = log(Normal - Quantiles).
    sigma  <- exp(1)    # 'sigma' argument.
    percentile <- c(25, 50, 75, 95)  # some percentiles of interest.
 
 ## E1. Get 'mean' values.
   theta <- uninormalQlink(theta = eta, percentile = percentile,
                          sd = sigma, inverse = TRUE)  # Mu
   
 \dontrun{
 ## E2. Plot theta vs. eta, 'shape' fixed, for different percentiles.
plot(theta[, 1], eta, type = "l", las = 1, lty = 2, lwd = 3,
     ylim = c(-10, 10), xlim = c(-10, 10),
     main = "uninormalQlink(theta; shape), fixed 'shape'.",
     xlab = "Theta (scale)", ylab = "uninormalQlink")
abline(v = 0, h = 0, col = "red")
lines(theta[, 2], eta, lty = 2, lwd = 3, col = "blue")
lines(theta[, 3], eta, lty = 2, lwd = 3, col = "orange")
lines(theta[, 4], eta, lty = 2, lwd = 3, col = "red")
legend("bottomright", c("25th Perc", "50th Perc", "75th Perc", "95th Perc"),
      col = c("black", "blue", "orange", "red"), lty = c(2, 2, 2, 2),
      lwd = rep(3, 4))

 }
 
 ## E3. uninormalQlink() and its inverse ##
    etabis  <- uninormalQlink(theta = theta, percentile = percentile,
                              sd = sigma, inverse = FALSE)
    my.diff <- eta - etabis
    summary(my.diff)     # Zero
}
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