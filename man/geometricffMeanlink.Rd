\name{geometricffMeanlink}
\alias{geometricffMeanlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

   Link functions for the mean of 1--parameter 
   discrete distributions: The Geometric Distribution.

}
\description{
 
  Computes the \code{geometricffMeanlink} transformation, including its inverse
  and the first two derivatives.

}
\usage{
    geometricffMeanlink(theta, bvalue = NULL, inverse = FALSE, 
                        deriv = 0, short = TRUE, tag = FALSE)
}
\arguments{
  \item{theta}{
   Numeric or character.
   See below for further details.
  
  }
  \item{bvalue, inverse, deriv, short, tag}{
  Details at \code{\link[VGAM:Links]{Links}}
  
  }
}

\details{
  This is a natural link function to model the mean of the (discret)
  geometric distribution, \code{\link[VGAM:geometric]{geometric}},
  defined as the logarithmm of its mean, i.e.,
  \deqn{ \eta = -\log \frac{p}{1 - p} = -{\tt{logit}}(p).}{
         \eta = - log (prob / (1 - prob)) = -logit(prob).
  }
  
  Here, \eqn{p}{prob} is the probability of succes, as in 
  \code{\link[VGAM:geometric]{geometric}}. 
  
  While this link function can be used to model any parameter lying in 
  \eqn{(0, 1)}, it is particularly useful for event-rate geometric data 
  where the mean can be written in terms of some rate of events, 
  say \eqn{\lambda = \lambda(\mathbf{x})}{\lambda(x)}, as
  %
  \deqn{\mu = \lambda(\mathbf{x}) t,}{\mu = \lambda(x) t ,}
  %
  and the time \eqn{t} (as \eqn{\log t}{log(t)}) can be easily 
  incorporated in the analysis as an offset.
  
  Under this link function the domain set for \eqn{p}{prob} 
  is \eqn{(0, 1)}. Hence, values of \eqn{\rho}{rho} too
  close to the extremes, or out of range will result 
  in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}. 
  Use argument \code{bvalue} to adequately replace them before 
  computing the link function.
  
  If \code{theta} is a character, arguments \code{inverse} and 
  \code{deriv} are disregarded.
  
}
\value{

  For \code{deriv = 0}, the \code{geometricffMeanlink} transformation of
  \code{theta} when \code{inverse = FALSE}. When 
  \code{inverse = TRUE} then \code{theta} becomes \eqn{\eta}, and
  \code{exp(-theta) / (exp(-theta) - 1)} is returned.
  
  
  For \code{deriv = 1}, \emph{d} \code{eta} / \emph{d} \code{theta},
  if \code{inverse = FALSE}, else 
  the reciprocal \emph{d} \code{theta} / \emph{d} \code{eta}
  as a function of \code{theta}.


  For \code{deriv = 2} the second order derivatives 
  are correspondingly returned.

}
\note{

  This function may return \code{Inf} of \code{-Inf} for values of 
  \eqn{p}{prob} too close to \eqn{0} and \eqn{1} respectively.
 
}
\section{Warning}{
   Numerical instability may occur if covariates are used leading to 
   values of \eqn{p}{prob} out of range. Try to overcome this by using
   argument \code{bvalue}.
  
}
\author{

 V. Miranda and T. W. Yee

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[VGAM:geometric]{geometric}},
 \code{\link[VGAM:Links]{Links}}, 
 \code{\link[VGAM:logit]{logit}},
 \code{\link[VGAMextra:logffMeanlink]{logffMeanlink}}.
 
}
\examples{

### Example 1  ###
my.probs <- ppoints(100)
geol.inv <-
    geometricffMeanlink(theta = geometricffMeanlink(theta = my.probs), # the inverse
                                inverse = TRUE) - my.probs
summary(geol.inv)     ## zero

###  Example 2. Special values of 'prob'  ###
my.probs <- c(-Inf, -2, -1, 0, 0.25, 0.75, 1.0, 5, Inf, NaN, NA) 
rbind(probs = my.probs, 
      geoffMlink = geometricffMeanlink(theta = my.probs),
      inv.geoffl = geometricffMeanlink(theta = my.probs, inverse = TRUE))


###  Example 3 Some probability link functions  ###
%
\donttest{
my.probs <- ppoints(100)

par(lwd = 2)
plot(my.probs, logit(my.probs), xlim = c(-0.1, 1.1), ylim = c(-5, 8),
     type = "l", col = "limegreen", 
     ylab = "transformation", las = 1, main = "Some probability link functions")
lines(my.probs, geometricffMeanlink(my.probs), col = "gray50")
lines(my.probs, logffMeanlink(my.probs), col = "blue")
lines(my.probs, probit(my.probs), col = "purple")
lines(my.probs, cloglog(my.probs), col = "chocolate")
lines(my.probs, cauchit(my.probs), col = "tan")
abline(v = c(0.5, 1), lty = "dashed")
abline(v = 0, h = 0, lty = "dashed")
legend(0.1, 8, 
      c("geometricffMeanlink", "logffMeanlink","logit", "probit", "cloglog", "cauchit"), 
      col = c("gray50", "blue", "limegreen", "purple", "chocolate", "tan"), lwd = 1, cex = 0.5)
par(lwd = 1) 
 }
}
