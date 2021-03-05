\name{weibullQlink}
\alias{weibullQlink}

\title{
 
   Weibull Quantile regression:
   Link function for the quantiles of the Weibull distribution.
    
}
\description{
   %
   Computes the \code{weibullQlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
       weibullQlink(theta, percentile = stop("Enter percentiles."),
                    shape = NULL, wrt.param = NULL,
                    bvalue = NULL, inverse = FALSE,
                    deriv = 0, short = TRUE, tag = FALSE)
                       
}
\arguments{
    \item{theta}{
     Numeric or character. Same as
     \code{\link[VGAMextra:uninormalQlink]{uninormalQlink}}
    
    }
    \item{percentile}{
    Same as \code{\link[VGAMextra:uninormalQlink]{uninormalQlink}}.
    Denoted below as \code{perc}.
    
    }
    \item{shape}{
    Numeric, positive. The shape parameter, required. 
    
    }
    \item{wrt.param}{
    Same as in \code{\link[VGAMextra:uninormalQlink]{uninormalQlink}} 
    %Positive integer, either \eqn{1} or \eqn{2}. The partial derivatives
    %are computed with respect to one of the two linear predictors 
    %involved with this link. Further details listed below.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{
   
   The ordinary 
   scale--shape Weibull quantiles are directly modelled by this link,
   aka \code{weibullQlink} transformation.
   It can only be used within
  \code{\link[VGAMextra:weibullRff]{weibullRff}}
  as the first linear predictor, \eqn{\eta_1},
  and is defined as
  %This is the link for the mean of the 2--parameter
  %Weibull distribution, as in 
  %\code{\link[VGAM:weibullR]{weibullR}}.
  %
  %
  \deqn{ \tt{weibullQlink}(\beta; \alpha) = \eta_1(\beta; \alpha) = 
                \log \{\beta \cdot [(-\log(1 - perc))^{(1/\alpha)}]\},}{
         weibullQlink(\beta; \alpha) = 
         \eta_1(\beta; \alpha) = 
         log {\beta * [(-log(1 - perc))^(1/\alpha)]},
         }
  %
 for \emph{given} \eqn{\alpha} ('shape' parameter) where
 \eqn{\beta > 0} is the \emph{scale} parameter.
%
 % a wrapper of 
  %
  \code{weibullQlink} is expressly a function of \eqn{\beta}, i.e.
  \eqn{\theta}, therefore \eqn{\alpha} (\emph{shape})
  must be entered at every call.
  
  Numerical values of \eqn{\alpha} or \eqn{\beta} out of range may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{weibullQlink} transformation of
   \code{theta}, i.e. \eqn{\beta}, when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then \eqn{\theta} becomes \eqn{\eta},
   and the inverse, \code{exp[theta} - 
   \eqn{(1 / \alpha) log(-log(1 - perc))}\code{]},
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
   
   
   B) For \code{inverse = TRUE}, this link returns
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
   terms of \code{theta} are similarly returned.
   
  %By default, the linear predictors in
  %\code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} are
  %\eqn{\eta_1 = \log \frac{\alpha}{\beta}}{\eta1 = log (\alpha / \beta)},
  %and \eqn{\eta_2 = \log \beta}{\eta2 = log \alpha}.
}
\section{Note}{
  
  See \code{\link[VGAMextra:weibullMlink]{weibullMlink}}.
  %Numerical instability may occur for values \code{theta} and
  %\code{shape} too close
  %to zero. Use argument \code{bvalue} to replace the former only before
  %computing the link.
  
  %If \code{theta} is character, then arguments \code{inverse} and
  %\code{deriv} are ignored. See \code{\link[VGAM:Links]{Links}}
  %for further details.
  
  
}
\author{
  V. Miranda and Thomas W. Yee.

}
\references{
Miranda & Yee (2021) 
\emph{Two--Parameter Link Functions,
With Application to Negative Binomial, Weibull and
Quantile Regression}. In preparation.

%Yee, T. W. (2015)
%\emph{Vector Generalized Linear and Additive Models: 
%With an Implementation in R}. 
%Springer, New York, USA.

}
\seealso{
    \code{\link[VGAMextra:weibullRff]{weibullRff}},
    \code{\link[VGAMextra:Q.reg]{Q.reg}},
    \code{\link[VGAM:weibullR]{weibullR}},
    \code{\link[VGAM:weibullR]{weibmeanlink}},
    \code{\link[VGAM:Links]{Links}}.
}
\examples{
    eta <- seq(-3, 3, by = 0.1) # this is eta = log(Weibull-quantiles).
    shape  <- exp(1)    # 'shape' argument.
    percentile <- c(25, 50, 75, 95)  # some percentiles of interest.
 
 ## E1. Get 'scale' values. Gives a warning (not of the same length) !
   theta <- weibullQlink(theta = eta, percentile = percentile,
                         shape = shape, inverse = TRUE)  # Scale
   
 \dontrun{
 ## E2. Plot theta vs. eta, 'shape' fixed, for different percentiles.
plot(theta[, 1], eta, type = "l", lwd = 3,
     ylim = c(-4, 4), 
     main = paste0("weibullQlink(theta; shape = ", round(shape, 3), ")"), 
     xlab = "Theta (scale)", ylab = "weibullQlink")
abline(h = -3:3, v = 0, col = "gray", lty = "dashed")
lines(theta[, 2], eta, lwd = 3, col = "blue")
lines(theta[, 3], eta, lwd = 3, col = "orange")
lines(theta[, 4], eta, lwd = 3, col = "red")
legend("bottomright", c("25th Perc", "50th Perc", "75th Perc", "95th Perc"),
      col = c("black", "blue", "orange", "red"),
      lwd = rep(3, 4))
 }
 
 ## E3. weibullQlink() and its inverse ##
    etabis  <- weibullQlink(theta = theta, percentile = percentile,
                            shape = shape, inverse = FALSE)
    summary(eta - etabis)     # Should be 0 for each colum (percentile)
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