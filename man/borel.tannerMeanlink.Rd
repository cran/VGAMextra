\name{borel.tannerMlink}
\alias{borel.tannerMlink}
\alias{borel.tannerMeanlink}
\title{

   Link functions for the mean of 1--parameter 
   discrete distributions: The Borel--Tanner distribution.

}
\description{
 
  Computes the \code{borel.tannerMlink} transformation, 
  its inverse and the first two derivatives.

}
\usage{
%
    borel.tannerMlink(theta, Qsize = 1, 
                      bvalue = NULL, inverse = FALSE, 
                      deriv = 0, short = TRUE, tag = FALSE)
%
}
\arguments{
  \item{theta}{
  Numeric or character. See below for further details.
  
  }
  \item{Qsize}{
  A positive integer. It is called \eqn{Q}.
  Same as \code{\link[VGAM:borel.tanner]{borel.tanner}}.
  Default it 1.
  
  }
  \item{bvalue, inverse, deriv, short, tag}{
  Details at \code{\link[VGAM:Links]{Links}}
  
  }
}

\details{
  
  As with 
  \code{\link[VGAMextra:zetaffMlink]{zetaffMlink}} or
  \code{\link[VGAMextra:yulesimonMlink]{yulesimonMlink}},
  this link function is part of a set of link functions in
  \pkg{VGAM} developed under a common methodology: by taking the 
  logarithm of the mean of the corresponding distribution. 
  
  In particular, this link function emerges by computing the 
  logarithm of the mean of the Borel--Tanner distribution.
  It is defined as
  %
  \deqn{  \tt{borel.tannerMlink}(a) = - \log (Q^{-1} - aQ^{-1}),}{
          borel.tannerMlink(a) = - log ( Q^(-1) - a Q^(-1) ),
  }
  %
  where \eqn{a}, \eqn{0 < a < 1}, is a scale parameter 
  as in \code{\link[VGAM:borel.tanner]{borel.tanner}}. 
  
  The domain set of \code{borel.tannerMlink} is the open interval
  \eqn{(0, 1)}, except when \code{inverse = TRUE} and \code{deriv = 0}.
  See below for further details about this.
  Moreover, unlike \code{\link[VGAMextra:zetaffMlink]{zetaffMlink}}
  or \code{\link[VGAMextra:posPoiMlink]{posPoiMlink}}, the
  inverse of \code{borel.tannerMlink} can be written in
  closed--form.
  
  Values of \eqn{a} (i.e. \code{theta}) out of range will result in 
  \code{NaN} of \code{NA}.

  If \code{theta} is a character, arguments \code{inverse} and
  \code{deriv} are discarded.
  
}
\value{

  For \code{deriv = 0}, the \code{borel.tannerMlink} transformation of
  \code{theta}, if \code{inverse = FALSE}. When \code{inverse = TRUE},
  \code{theta} becomes \eqn{\eta} and the inverse of 
  \code{borel.tannerMlink}, given by
  %
  \deqn{ 1 - \frac{Q}{ e^{\eta} },}{1 - Q / exp(\eta),}
  %
  is returned. Here, the domain set changes to \eqn{(0, \infty)}.
  
  For \code{deriv = 1}, \emph{d} \code{eta} / \emph{d} \code{theta} 
  as a function of \code{theta} if \code{inverse = FALSE}, else 
  the reciprocal \emph{d} \code{theta} / \emph{d} \code{eta}.
  
  Similarly, when \code{deriv = 2} the second order derivatives
  in terms of \code{theta} are returned.
  %Specially the zetafflink function as well as its derivatives
  %are graphically delimited for specific asymptotes.
  %Consequently, the mathematical limit of this link function is 
  %returned for special values of \code{theta}, e.g. 
  %for \code{theta} = \eqn{\infty}. See example 2 below.
  
}
\author{
 V. Miranda and T. W. Yee

}
\note{
 
  The vertical line \eqn{a = 1} is an asymptote for this link 
  function, which sharply grows for values of \eqn{a} too close to 1.0
  from the left. For such cases, \code{Inf} might result when computing 
  \code{borel.tannerMlink}.
  
  This link function is useful to model any parameter
  in \eqn{(0, 1)}. Then, some problems may occur if there
  are covariates causing out of range values.
  
}

\section{References}{
Haight, F. and Brueuer, M. A. (1960) The Borel--Tanner distribution.
\emph{Biometrika}, \bold{47}, 143--150.

}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[VGAM:borel.tanner]{borel.tanner}},
  \code{\link[VGAMextra:yulesimonMlink]{yulesimonMlink}},
 \code{\link[VGAMextra:zetaffMlink]{zetaffMlink}},
 \code{\link[VGAMextra:posPoiMlink]{posPoiMlink}},
 \code{\link[VGAM:Links]{Links}}.
 
}
\examples{

## Example 1. Special values for theta (or eta, accordingly)  ##
a.par <- c(0, 1:10/10, 20, 1e1, Inf, -Inf, NaN, NA) 

# The borel.tannerMlink transformation and the first two derivatives.
print(rbind(a.par, 
 deriv1 = borel.tannerMlink(theta = a.par, inverse = FALSE, deriv = 1),
 deriv2 = borel.tannerMlink(theta = a.par, inverse = FALSE, deriv = 2)),
 digits = 2)
  
# The inverse of 'borel.tannerMlink()' and the first two derivatives.
# 'theta' turns into 'eta'.
print(rbind(a.par, 
    Invderiv1 = borel.tannerMlink(theta = a.par, inverse = TRUE, deriv = 1),
    Invderiv2 = borel.tannerMlink(theta = a.par, inverse = TRUE, deriv = 2)),
    digits = 2)

  
## Example 2 ##
a.param <- c(0, 1, 5, 10, 1e2, 1e3) 
rbind(a.values = a.param, 
      inv.BT = borel.tannerMlink(theta = a.param, inverse = TRUE))
    
data.inv <- borel.tannerMlink(borel.tannerMlink(a.param, inv = TRUE)) - a.param
summary(data.inv)                      ## Should be zero



## Example 3. Some link functions in VGAM with domain set (0, 1) ##
a.param <- ppoints(100)
%
\donttest{
par(lwd = 2)
plot(a.param, borel.tannerMlink(a.param), ylim = c(-5, 7), xlim = c(-0.01, 1.01),
     type = "l", col = "gray10", ylab = "transformation",
     las = 1, main = "Some probability link functions")
lines(a.param, logffMlink(a.param), col = "blue")
lines(a.param, logitlink(a.param), col = "limegreen")
lines(a.param, probitlink(a.param), col = "purple")
lines(a.param, clogloglink(a.param), col = "chocolate")
lines(a.param, cauchitlink(a.param), col = "tan")
abline(v = c(0.5, 1), lty = "dashed")
abline(v = 0, h = 0, lty = "dashed")
legend(0.05, 7, c("borel.tanneMlink", "logffMlink", "logitlink", "probitlink",
                  "clogloglink", "cauchitlink"),
       col = c("gray10", "blue", "limegreen", "purple", "chocolate", "tan"),
      lwd = 1)
par(lwd = 1) 
}
%
}