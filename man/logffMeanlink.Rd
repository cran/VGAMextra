\name{logffMeanlink}
\alias{logffMeanlink}
\alias{logfflink.inv.deriv0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

   Link functions for the mean of 1--parameter 
   discrete distributions: The Logarithmic Distribuion.

}
\description{
 
  Computes the \code{logffMeanlink} transformation, including its inverse
  and the first two derivatives.

}
\usage{
%
    logffMeanlink(theta, bvalue = NULL, 
                  alg.roots = c("Newton-Raphson", "bisection")[1],
                  inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
%
}
\arguments{
  \item{theta}{
  Numeric or character. This is \eqn{\theta} by default although could
  be \eqn{\eta} depending on other parameters. See below for details.
  
  }
  \item{bvalue}{
  This is a boundary value. See below. 
  Also refer to \code{\link[VGAM:Links]{Links}} for additional
  details.
  
  }
  \item{alg.roots}{
  Character. The iterative algorithm to find the inverse of this link
  function. Default is the first (Newton--Raphson). Optionally, the 
  bisection method is also available. See below for more details.
  
  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link[VGAM:Links]{Links}}
  
  }
}

\details{
  This link function arises as a natural link function for the mean,
  \eqn{\mu}, of the logarithmic (or log-series) distribution, 
  \code{\link[VGAM:logff]{logff}}. It is defined for any value 
  of the shape parameter \eqn{s} (i.e. \code{theta} in the 
  VGLM/VGAM context), \eqn{0 < s < 1}, as the
  logarithm of \eqn{\mu = \mu(s)}. It can be easily shown that 
  \code{logffMeanlink} is the difference of two common link functions: 
  \code{\link[VGAM:logitlink]{logitlink}} and 
  \code{\link[VGAM:clogloglink]{clogloglink}}.
  
  It is particularly usefull for event--rate data where the expected 
  number of events can be modelled as
  %
  \deqn{\mu = \mu(s) = \lambda t.}
  %
  Here \eqn{\lambda} is the standardized mean (or event-rate) 
  per unit time, \eqn{t} is the timeframe observed, 
  whereas \eqn{\mu} and \eqn{s} are the mean and the 
  shape parameter of the logarithmic distribution respectively.
  The logarithm is then applied to both sides so 
  that \eqn{t} can be incorporated in the analysis as an offset. 
  
  While \code{logffMeanlink} is not the canonical link function of 
  the logarithmic distribution, it is certainly part of the
  canonical link, given by the composite
  %
  \deqn{\log \circ~(g^{-1}) \circ \log,}{
               log o (g^(-1)) o log,}
  %
  where \eqn{g^{-1}}{g^(-1)} denotes the inverse of
  \code{logffMeanlink}.
  
  The domain set of this link function is \eqn{(0, 1)}. 
  Therefore, values of \code{theta} (that is \eqn{s}) too close to 
  \eqn{0} or to \eqn{1} or out of range will result in \code{Inf}, 
  \code{-Inf}, \code{NA} or \code{NaN}. Use argument \code{bvalue} 
  to adequately replace them before computing the link function.
  
  Particularly, if \code{inverse = TRUE} and \code{deriv = 0},
  then \eqn{s} becomes \eqn{\eta}{eta}, and therefore the domain set 
  turns to \eqn{(0, \infty)}{(0, \infty)}.
  
  If \code{theta} is a character, then arguments \code{inverse} and
  \code{deriv} are disregarded.

}
\value{

  For \code{deriv = 0}, the \code{logffMeanlink} transformation of 
  \code{theta}, i.e., \code{logitlink(theta) - clogloglink(theta)}, if 
  \code{inverse = FALSE}. 
  
  When \code{inverse = TRUE} the vector entered at \code{theta} 
  becomes \eqn{\eta} and, then, this link function 
  returns a unique vector \eqn{\theta_{\eta}}{\theta[\eta]} such that 
  %
  \deqn{ {\tt{logffMeanlink}} (\theta_{\eta}) = \eta,}{
          logffMeanlink(\theta[\eta]) = \eta, }
  %
  i.e., the inverse image of \eqn{\eta}.
  Specifically, the inverse of \code{logffMeanlink} cannot be written in 
  closed--form, then the latter is equivalent to search for the roots 
  of the function 
  %
  \deqn{ {\tt{logff.func}}(\theta) = {\tt{logffMeanlink}}(\theta) - \eta }{
              logff.func(\theta) = logffMeanlink(\theta) - \eta} 
  %
  as a function of \eqn{\theta}. To do this, the auxiliary 
  function \code{logff.func} is internally generated. 
  Then, with the method established at \code{alg.roots}, 
  either Newton--Raphson or bisection, this link function
  approximates and returns the inverse image 
  \eqn{\theta_{\eta}}{\theta[\eta]} (of given \eqn{\eta}), which 
  plays the role of the inverse of \code{logffMeanlink}. In particular, 
  for \eqn{\eta = 0} and \eqn{\eta =}\code{Inf}, it returns \eqn{0} 
  and \eqn{1} respectively.

  For \code{deriv = 1}, \emph{d} \code{eta} / \emph{d} \code{theta} 
  as a function of \code{theta} if \code{inverse = FALSE}, else 
  the reciprocal \emph{d} \code{theta} / \emph{d} \code{eta}.


  Similarly, when \code{deriv = 2} the second order derivatives 
  are correspondingly returned.

  Both, first and second derivatives, can be written in closed--form.

}
\section{Warning}{

  \code{\link[VGAMextra:logffMeanlink]{logffMeanlink}} is a monotonically
  increasing, convex, and strictly positive function in \eqn{(0, 1)} such
  that the horizontal axis is an asymptote. Therefore, when the inverse 
  image of \eqn{\eta} is required, each entry of \eqn{\eta} 
  (via argument \code{theta}) must be non-negative so that 
  \eqn{{\tt{logff.func}(\theta; \eta) = \tt{logffMeanlink}}(\theta) - \eta}{
        logff.func(\theta; \eta) = logffMeanlink(\theta) - \eta} 
  is \emph{shifted down}. This fact allows this function to uniquely 
  intersect the horizontal axis which guarantees to iteratively find 
  the corresponding root \eqn{\theta_{\eta}}{\theta[\eta]}, i.e., the
  inverse image of \eqn{\eta}. Else, \code{NaN} will be returned.
  
  See example 3. It is the plot of \code{logffMeanlink} in \eqn{(0, 1)}
  for \eqn{\eta = 1.5}.
  
  Besides, the vertical straight line \code{theta}\eqn{ = 1} is also
  an asymptote. Hence, this link function may grow sharply for
  values of \code{theta} too close to \eqn{1}.
  See Example 4 for further details.
 
}
\author{
 V. Miranda and T. W. Yee

}
\note{

  To find the inverse image \eqn{\theta_{\eta}}{\theta[\eta]} of 
  \eqn{\eta}, either
  \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}} 
  or \code{\link[VGAM:bisection.basic]{bisection.basic}} is
  called.
 
  This link function can be used for modelling any parameter lying 
  between 0.0 and 1.0. Consequently, when there are covariates, 
  some problems may occur. For example, the method entered 
  at \code{alg.roots} to approximate the inverse image may converge
  at a slow rate. Similarly if the sample size is small, less than 20 say.
  Try another link function, as \code{\link[VGAM:logit]{logit}}, in such
  cases.
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[VGAM:logff]{logff}},
 \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}},
 \code{\link[VGAM:bisection.basic]{bisection.basic}},
 \code{\link[VGAM:Links]{Links}}, 
 \code{\link[VGAM:clogloglink]{clogloglink}},
 \code{\link[VGAM:logitlink]{logitlink}}.
 
}
\examples{

## Example 1  ##
set.seed(0906)
Shapes <- sort(runif(10)) 
logffMeanlink(theta = Shapes, deriv = 1)    ## d eta/d theta, as function of theta

logldata.inv <-
       logffMeanlink(theta = logffMeanlink(theta = Shapes), inverse = TRUE) - Shapes
                      
summary(logldata.inv)                   ## Should be zero


## Example 2 Some probability link funtions ##
s.shapes <- ppoints(100)
%
\donttest{
par(lwd = 2)
plot(s.shapes, logit(s.shapes), xlim = c(-0.1, 1.1), type = "l", col = "limegreen", 
     ylab = "transformation", las = 1, main = "Some probability link functions")
lines(s.shapes, logffMeanlink(s.shapes), col = "blue")
lines(s.shapes, probitlink(s.shapes), col = "purple")
lines(s.shapes, clogloglink(s.shapes), col = "chocolate")
lines(s.shapes, cauchitlink(s.shapes), col = "tan")
abline(v = c(0.5, 1), lty = "dashed")
abline(v = 0, h = 0, lty = "dashed")
legend(0.1, 4.5, c("logffMeanlink","logitlink", "probitlink", "clogloglink",
"cauchitlink"),
       col = c("blue", "limegreen", "purple", "chocolate", "tan"), lwd = 1)
par(lwd = 1) 
}
%

## Example 3. Plot of 'logffMeanlink()' with eta = 1.5.  ##
m.eta1.5 <- logffMeanlink(theta = s.shapes, deriv = 0) - 1.5
%
\donttest{
 plot(m.eta1.5 ~ s.shapes, type = "l", col = "limegreen", 
      las = 1, lty = 2, lwd = 3, xlim = c(-0.1, 1.0), ylim = c(-2, 3),
      xlab = "shape parameter, s, in (0, 1).", 
      ylab = "logffMeanlink(s) - 1.5",
      main = "logff.func(s; 1.5) = logffMeanlink(s) - 1.5, in (0, 1)")
 abline(h = 0, v = 0)
 abline(v = 1.0, lty = 2)
 }
%

## Example 4. Special values of theta, inverse = FALSE ##
s.shapes <- c(-Inf, -2, -1, 0.0, 0.25, 0.5, 1, 10, 100, Inf, NaN, NA) 
rbind(s.shapes, logffMeanlink(theta = s.shapes))
}