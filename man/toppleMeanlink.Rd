\name{toppleMlink}
\alias{toppleMlink}
\alias{toppleMeanlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 
   Link functions for the mean of 1--parameter 
   continuous distribution: The Topp--Leone distribution.
    
}
\description{
   %
   Computes the \code{toppleMlink} transformation, its inverse and
   the first two derivatives.
   %
}
\usage{
      toppleMlink(theta, bvalue = NULL, inverse = FALSE,
                  deriv = 0, short = TRUE, tag = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{theta}{
     Numeric or character.
     See \code{\link[VGAM:Links]{Links}} and below for further details.
    
    }
    \item{bvalue, inverse, deriv, short, tag}{
    See \code{\link[VGAM:Links]{Links}}.
    
    }
}
\details{
  
  The \code{toppleMlink} transformation arises as a link function to
  model the mean of the Topp--Leone distribution,
  \code{\link[VGAM:topple]{topple}}. It is defined as
  
  \deqn{ \eta = {\tt{logit}} \left( \left( 1 - \frac{4^{s} \Gamma(1 + s)^2}{ 
                \Gamma(2 + 2s)} \right) / sup.tp \right).}{
        
          \eta = logit ( (1 - 4^s * \Gamma(1 + s)^2 / 
                 \Gamma(2 + 2s) ) / sup.tp ).}
                 
  Here, \eqn{0 < s < 1} is a shape parameter as in
  \code{\link[VGAM:topple]{topple}}, whereas 
  \eqn{sup.tp} is the \emph{supremum} of 
  \deqn{  1 - \frac{4^{s} \Gamma(1 + s)^2}{ \Gamma(2 + 2s)}, }{
          1 - 4^s * \Gamma(1 + s)^2 / \Gamma(2 + 2s), 
  }
  in \eqn{(0, 1)}, as a function of \eqn{s}.
  
  For numerical values of \eqn{s} out of \eqn{(0, 1)}, this link may 
  result in \code{Inf}, \code{-Inf}, \code{NA} or \code{NaN}.

  
}
\value{
  
   For \code{deriv = 0}, the \code{toppleMlink} transformation of
   \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then \code{theta} becomes \eqn{\eta}, and
   the inverse transformation is required. However, it 
   can't be expressed in close form. Therefore, the approximate 
   \emph{inverse image} of entered \code{theta} computed by 
   \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}} 
   is returned.
   
   For \code{deriv = 1},
   \eqn{d} \code{eta} / \eqn{d} \code{theta} when \code{inverse = FALSE}.
   If \code{inverse = TRUE}, then 
   \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of 
   \code{theta}.
   
   %Similarly, when \code{deriv = 2}, the second derivatives in
   %terms of \code{theta} are returned.

}
\section{Note}{

  Values of \eqn{s} too close to zero or 1.0 may cause numerical 
  instability. Use argument \code{bvalue} to replace them before
  computing the link.
  
  If \code{theta} is character, then arguments \code{inverse} and
  \code{deriv} are ignored. See \code{\link[VGAM:Links]{Links}}
  for further details.
  
  
}
\author{
  V. Miranda and Thomas W. Yee.

}
\seealso{
    \code{\link[VGAM:topple]{topple}}, 
    \code{\link[VGAM:Links]{Links}},
    \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}}. 
     
     
}
\examples{
 ## E1. The toppleMlink() and its inverse ##
    theta <- ppoints(10)
    eta   <- toppleMlink(toppleMlink(theta = theta), inverse  =TRUE)
    summary(eta - theta)     # Zero


  ## E2. Some probability link functions ##
  \donttest{
my.probs <- ppoints(100)

par(lwd = 2)
plot(my.probs, logitlink(my.probs), xlim = c(-0.1, 1.1), ylim = c(-5, 8),
     type = "l", col = "limegreen", 
     ylab = "transformation", las = 1, main = "Some probability link functions")
lines(my.probs, geometricffMlink(my.probs), col = "gray50")
lines(my.probs, logffMlink(my.probs), col = "blue")
lines(my.probs, probitlink(my.probs), col = "purple")
lines(my.probs, clogloglink(my.probs), col = "chocolate")
lines(my.probs, cauchitlink(my.probs), col = "tan")
lines(my.probs, toppleMlink(my.probs), col = "black")
abline(v = c(0.5, 1), lty = "dashed")
abline(v = 0, h = 0, lty = "dashed")
legend(0.1, 8, 
      c( "toppleMlink", "geometricffMlink", "logffMlink",
      "logitlink", "probitlink", 
      "clogloglink", "cauchitlink"), 
      col = c("black", "gray50", "blue", "limegreen", "purple", "chocolate", "tan"),
      lwd = 1, cex = 0.5)
par(lwd = 1) 
 }

}
