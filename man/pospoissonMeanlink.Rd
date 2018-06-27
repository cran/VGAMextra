\name{posPoiMeanlink}
\alias{posPoiMeanlink}
\alias{pospoilink.inv.deriv0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Link functions for the mean of 1--parameter 
  discrete distributions: The Positive Poisson Distribution.

}
\description{
 
  Computes the \code{posPoiMeanlink} transformation, 
  its inverse and the first two derivatives.

}
\usage{
%
     posPoiMeanlink(theta, bvalue = NULL, 
                    alg.roots = c("Newton-Raphson", "bisection")[1],
                    inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
%
}
\arguments{
  \item{theta}{
  Numeric or character. This is \eqn{\theta}{theta} by default, although
  it becomes \eqn{\eta}{eta} sometimes, depending on the other parameters. 
  See below for further details.
  
  
  }
  \item{bvalue}{
  Details at \code{\link[VGAM:Links]{Links}}.
  
  }
  \item{alg.roots}{
  Character. The iterative method to find the inverse of this link
  function. Same as \code{\link[VGAMextra:zetaffMeanlink]{zetaffMeanlink}}.
  
  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link[VGAM:Links]{Links}}
  
  }
}

\details{
  
  This is a link function for the mean of the 
  positive Poisson distribution. It is defined as 
  %
  \deqn{ \eta = \tt{posPoiMeanlink}(\lambda) = 
           -\log (\lambda^{-1} - \lambda^{-1} e^{-\lambda}),}{
          \eta = posPoiMeanlink(\lambda) =  
             - log (\lambda^(-1)- \lambda^(-1) * exp(-\lambda)),
        }
  %
  where \eqn{\lambda > 0} stands for the single parameter of
  \code{\link[VGAM:pospoisson]{pospoisson}}, i.e. 
  \code{theta} in the VGLM/VGAM context.
  
  Notice, the mean of the positive Poisson is given by
  %
  \deqn{ \frac{\lambda}{1 - e^{-\lambda}}. }{
              \lambda / (1 - exp(-\lambda)).
  }
  %
  This link function comes up by taking the logarithm on both
  sides of this equation.
  
  The domain set for \eqn{\lambda} is \eqn{(0, \infty)}{(0, \infty)}. 
  Hence, non--positive values of \eqn{\lambda} will result in \code{NaN} 
  or \code{NA}. Use argument \code{bvalue} to properly replace 
  them before computing the link function.
  
  \code{posPoiMeanlink} tends to infinity as \eqn{\lambda} 
  increases. Specially, its inverse grows at a higher rate. Therefore, 
  large values of \eqn{\lambda} will result in \code{Inf} accordingly. 
  See example 2 below.
  
  If \code{theta} is a character, arguments \code{inverse} and
  \code{deriv} are disregarded.
  
}
\value{

  For \code{deriv = 0}, the \code{posPoiMeanlink} transformation of 
  \code{theta}, if \code{inverse = FALSE}. When \code{inverse = TRUE},
  \code{theta} becomes \eqn{\eta} and the inverse of 
  \code{posPoiMeanlink} is needed but cannot be written in 
  closed--form. Instead this link function returns the approximated 
  inverse image of \eqn{\eta}, say \eqn{\theta_\eta}{\theta[\eta]}, 
  such that 
  %
  \deqn{  \tt{posPoiMeanlink}(\theta_\eta) = \eta.}{
              posPoiMeanlink(\theta[\eta]) =\eta.
  }
  %
  Here, \eqn{\theta_\eta}{\theta[\eta]} is iteratively computed
  as the unique root of the auxiliary function
  %
  \deqn{ f(\theta; \eta) = \tt{posPoiMeanlink}(\theta) - \eta,}{
         f(\theta; \eta) =  posPoiMeanlink(\theta) - \eta,
  }
  %
  as a function of \eqn{\theta}.
  This work is performed via Newton--Raphson or bisection,
  as per argument \code{alg.roots}. 
  
  For \code{deriv = 1}, \emph{d} \code{eta} / \emph{d} \code{theta} 
  as a function of \code{theta} if \code{inverse = FALSE}, else 
  the reciprocal \emph{d} \code{theta} / \emph{d} \code{eta}.
  
  Similarly, when \code{deriv = 2} the second order derivatives 
  are returned in terms of \code{theta}.
  
}
\section{Warning}{
  
  This link function is monotonic increasing in \eqn{(0, \infty)} so that
  the horizontal axis is an asymptote. Then, in order to assure the root 
  of the auxiliary
  %
  \deqn{ f(\theta; \eta) = \tt{posPoiMeanlink}(\theta) - \eta}{
         f(\theta; \eta) =  posPoiMeanlink(\theta) - \eta
  }
  %
  to be real, \eqn{\eta} must be positive. As a result, 
  \code{posPoiMeanlink} is \emph{shited--down} and hence intersecting 
  the horizontal axis uniquely.
 
}
\author{
 V. Miranda and T. W. Yee.

}
\note{
 
  This link function is useful to model any parameter
  in \eqn{(0, \infty)}. Some numerical issues may arise if there are 
  covariates causing negative values the parameter. 
  Try \code{\link[VGAM:identitylink]{identitylink}} alternatively.
  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[VGAM:pospoisson]{pospoisson}},
 \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}},
 \code{\link[VGAM:bisection.basic]{bisection.basic}},
 \code{\link[VGAM:Links]{Links}},
 \code{\link[VGAM:identitylink]{identitylink}}.
 
}
\examples{

## Example 1. Special values for theta (or eta, accordingly)  ##
m.lambda <- c(0, 0.5, 1, 10, 20, 25, 1e2, 1e3, Inf, -Inf, NaN, NA) 

# The 'posPoiMeanlink' transformation and the first two derivatives.
print(rbind(m.lambda, 
  deriv1 = posPoiMeanlink(theta = m.lambda, inverse = FALSE, deriv = 1),
  deriv2 = posPoiMeanlink(theta = m.lambda, inverse = FALSE, deriv = 2)),
  digits = 2)
  
# The inverse of 'posPoiMeanlink' and the first two derivatives.
print(rbind(m.lambda, 
  Invderiv1 = posPoiMeanlink(theta = m.lambda, inverse = TRUE, deriv = 1),
  Invderiv2 = posPoiMeanlink(theta = m.lambda, inverse = TRUE, deriv = 2)),
  digits = 2)
  

## Example 2. The inverse of 'posPoiMeanlink' ##
m.lambda <- c(0, 1, 5, 10, 1e2, 1e3) 
posPoiMeanlink(theta = posPoiMeanlink(m.lambda, inverse = TRUE))
pospoi.inv <- posPoiMeanlink(posPoiMeanlink(m.lambda, inverse = TRUE)) - m.lambda
                           
summary(pospoi.inv)                      ## Should be zero.


## Example 3. Plot of 'posPoiMeanlink' and its first two derivatives ##
## inverse = FALSE, deriv = 0, 1, 2. ##
%
\donttest{
m.lambda <- seq(0, 35, by = 0.01)[-1]
y.lambda <- posPoiMeanlink(theta = m.lambda, deriv = 0)
der.1 <- posPoiMeanlink(theta = m.lambda, deriv = 1)
der.2 <- posPoiMeanlink(theta = m.lambda, deriv = 2)

plot(y.lambda ~ m.lambda, col = "black", 
     main = "log(mu), mu = E[Y], Y ~ pospoisson(lambda).",
     ylim = c(-1, 10), xlim = c(-1, 26),
     lty = 1, type = "l", lwd = 3)
abline(v = 0, h = 0, col = "gray50", lty = "dashed")

lines(m.lambda, der.1, col = "blue", lty = 5, lwd = 3)
lines(m.lambda, der.2, col = "chocolate", lty = 4, lwd = 3)
legend(5, 9, legend = c("posPoiMeanlink", "deriv = 1", "deriv = 2"),
       col = c("black", "blue", "chocolate"), lty = c(1, 5, 4), lwd = c(3, 3, 3))
 }
}