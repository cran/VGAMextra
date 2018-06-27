\name{zetaffMeanlink}
\alias{zetaffMeanlink}
\alias{zetafflink.inv.deriv0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

 Link functions for the mean of 1--parameter 
 discrete distributions: The Zeta Distribution.

}
\description{
 
  Computes the \code{zetaffMeanlink} transformation, including its inverse
  and the first two derivatives.

}
\usage{
%
  zetaffMeanlink(theta, bvalue = NULL, 
                 alg.roots = c("Newton-Raphson", "bisection")[1],
                 inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE)
%             
}
\arguments{
  \item{theta}{
  Numeric or character. This is \eqn{\theta}{theta} by default, although
  it can be \eqn{\eta}{eta} sometimes, depending on the other parameters. 
  See below for further details.
  
  }
  \item{bvalue}{
  Details at \code{\link[VGAM:Links]{Links}}.
  
  }
  \item{alg.roots}{
  Character. The iterative method to find the inverse of this link
  function. Default is Newton--Raphson. Optionally, the bisection 
  method is also available.
  
  }
  \item{inverse, deriv, short, tag}{
  Details at \code{\link[VGAM:Links]{Links}}
  
  }
}

\details{

  This is a link function for the zeta distribution,
  \code{\link[VGAM:zetaff]{zetaff}}, which emerges by applying 
  the logarithm transformation to its mean. Specifically, 
  assume \eqn{Y} follows a zeta distribution with shape parameter 
  \eqn{s} (this is \code{theta} in the VGLM/VGAM framework). 
  Then, the mean of \eqn{Y} is 
  %
  \deqn{ \mu = \frac{\zeta (s)}{\zeta (s + 1)}, }{ 
         \mu = zeta(s) / zeta(s + 1), }
  %
  provided \eqn{s > 1}, where \eqn{\zeta}{zeta} is the Riemann`s zeta
  function computed by \code{\link[VGAM:zeta]{zeta}}. The notation 
  adopted here conforms with \code{\link[VGAM:zetaff]{zetaff}}
  in terms of the density of the zeta distribution.
  
  The \code{zetaffMeanlink} transformation is then given by
  %
  \deqn{\eta = \tt{zetaffMeanlink}(s) = \log \frac{\zeta (s)}{\zeta (s + 1)}.}{
        \eta = zetaffMeanlink(s) = log (zeta(s) / zeta(s + 1)). }
  %
  It is particularly useful when modelling event--rate data where 
  the expected number of events, \eqn{\mu}, can be modelled as
  %
  \deqn{\mu = \lambda t.}{\mu = \lambda t.}
  %
  Specifically, \eqn{\lambda} is a standardized mean per unit--time, and 
  \eqn{t} is the observed timeframe.
  
  The domain set for \eqn{s}, i.e. \eqn{\theta}{theta}, is 
  \eqn{(1, \infty)}{(1, \infty)}. Hence, either large values of 
  \eqn{s}, or those too close to \eqn{1} from the right, or out of 
  range will result in \code{Inf}, \code{-Inf}, \code{NA} or 
  \code{NaN}. Use argument \code{bvalue} to adequately replace 
  them before computing the link function.
  %See \code{\link[VGAM:Links]{Links}} for further details about this.
  
  WARNING: While in \code{\link[VGAM:zetaff]{zetaff}} the parameter 
  \eqn{s} lies in \eqn{(1, \infty)}{(0, \infty)}, 
  \code{zetaffMeanlink} will be real when \eqn{s > 1}. Consequently, 
  for any VGLM fitted via \code{\link[VGAM:zetaff]{zetaff}} using this 
  link function, numerical problems will take place if any 
  \eqn{s} value lies between 0.0 and 1.0 at any iteration. Use 
  optional link functions like \code{\link[VGAM:loge]{loge}}.
  
  When \code{inverse = TRUE} and \code{deriv = 0},
  \eqn{s} changes into \eqn{\eta}{eta}, and therefore the domain set 
  (only in this case) turns into \eqn{(0, \infty)}{(0, \infty)}.
  See below for further details.
  
  If \code{theta} is a character, arguments \code{inverse} and
  \code{deriv} are disregarded.
}
\value{

  For \code{deriv = 0}, the \code{zetaffMeanlink} transformation of 
  \code{theta}, if \code{inverse = FALSE}. When \code{inverse = TRUE},
  \code{theta} becomes \eqn{\eta}, and then the inverse of \code{zetaffMeanlink}
  is required. However, it cannot be written in closed--form. Instead,
  the inverse image of \eqn{\eta}, say \eqn{\theta_\eta}{\theta[\eta]}, 
  is returned. That is, a unique vector \eqn{\theta_\eta}{\theta[\eta]} 
  such that 
  %
  \deqn{ {\tt{zetaffMeanlink}} (\theta_\eta) = \eta.}{
          zetaffMeanlink(\theta[\eta]) = \eta. }
  %
  This process is equivalent to find the root, 
  \eqn{\theta_\eta}{\theta[\eta]}, of the function 
  \eqn{ {\tt{zetaffMeanlink}} (\theta) - \eta, }{
          zetaffMeanlink(\theta) - \eta, }
  which is internally carried out via the method entered
  at \code{alg.roots}. Options available are 
  ``\code{Newton-Raphson}'' and ``\code{bisection}''.
  
  For \code{deriv = 1}, \emph{d} \code{eta} / \emph{d} \code{theta} 
  as a function of \code{theta} if \code{inverse = FALSE}, else 
  the reciprocal \emph{d} \code{theta} / \emph{d} \code{eta}.
  
  Similarly, when \code{deriv = 2} the second order derivatives 
  are returned accordingly.
  
  The first two derivatives of the Riemman's zeta function 
  are computed by \code{\link[VGAM:zeta]{zeta}}.
  
  Besides, the \code{zetaffMeanlink} function as well as its derivatives
  are graphically delimited for specific asymptotes.
  Consequently, the mathematical limit of this link function is 
  returned for special values of \code{theta}, e.g. 
  for \code{theta} = \eqn{\infty}. See example 2 below.
  
}
\section{Warning}{
  
  Where the inverse image of \eqn{\eta}, \eqn{\theta_\eta}{\theta[\eta]},
  is required, values entered at \code{theta} (becoming \eqn{\eta}) must
  be non-negative. The reason is that the \code{zetaffMeanlink}
  transformation is decreasing but strictly positive in \eqn{(1, \infty)} 
  asymptotically approaching to the horizontal axis. In this way, the 
  \emph{shifted--down} \code{zetaffMeanlink} function
  %
  \deqn{{\tt{zetaff.func}}(\theta | \eta) = 
           {\tt{zetaffMeanlink}}(\theta) - \eta}{
           zetaff.func(\theta | \eta) = zetaffMeanlink(\theta) - \eta }
  %
  uniquely intersects the horizontal axis and hence the inverse image 
  computed by ``\code{Newton-Raphson}'' or ``\code{bisection}'' 
  will be a real
  number.

 
}
\author{
 V. Miranda and T. W. Yee

}
\note{
 
  Overall, this link function is useful to model any parameter
  lying in \eqn{(1, \infty)}, specially if the theoretical mean
  can be written as \eqn{\mu = \lambda t}, as stated above.
  As a result, some problems may arise if there are covariates.
  Try another link function if any issue, such 
  as \code{\link[VGAM:loglog]{loglog}}.
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[VGAM:zetaff]{zetaff}},
 \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}},
 \code{\link[VGAM:bisection.basic]{bisection.basic}},
 \code{\link[VGAM:zeta]{zeta}},
 \code{\link[VGAM:loge]{loge}},
 \code{\link[VGAM:Links]{Links}}.
 
}
\examples{

## Example 1  ##
Shapes <- 1:10 + 0.1 
zetaffMeanlink(theta = Shapes, deriv = 1)  ## d eta/d theta, as function of theta

zetafflk.inv <-
 zetaffMeanlink(theta = zetaffMeanlink(theta = Shapes), inverse = TRUE) - Shapes
 
summary(zetafflk.inv)                      ## Should be zero


## Example 2. Special values of theta, inverse = FALSE ##
Shapes <- c(-Inf, -1, 0.5, 1, 1.5, 10, 100, Inf, NaN, NA) 
print(rbind(Shapes, zetaffMeanlink  = zetaffMeanlink(theta = Shapes),
            inv.zfflink = zetaffMeanlink(theta = Shapes, inverse = TRUE)), digits = 3)


## Example 3. Plot of 'zetaffMeanlink()' and its first two derivatives ##
## inverse = FALSE, deriv = 0, 1, 2 ##
\donttest{
Shapes <- seq(1, 20, by = 0.01)[-1]
y.shapes <- zetaffMeanlink(theta = Shapes, deriv = 0)
der.1 <- zetaffMeanlink(theta = Shapes, deriv = 1)
der.2 <- zetaffMeanlink(theta = Shapes, deriv = 2)

plot(y.shapes ~ Shapes, 
     col = "black", main = "log(mu), mu = E[Y], Y ~ Zeta(s).",
     ylim = c(-5, 10), xlim = c(-0.1, 5), lty = 1, type = "l", lwd = 3)
abline(v = 1.0, col = "orange", lty = 2, lwd = 3)
abline(v = 0, h = 0, col = "gray50", lty = "dashed")

lines(Shapes, der.1, col = "blue", lty = 5)
lines(Shapes, der.2, col = "chocolate", lty = 4)
legend(2, 7, legend = c("zetaffMeanlink", "deriv = 1", "deriv = 2"),
       col = c("black", "blue", "chocolate"), lty = c(1, 5, 4), lwd = c(3, 1, 1))
 }

}