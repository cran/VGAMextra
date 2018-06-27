\name{normal1sdQlink}
\alias{normal1sdQlink}
\title{
  Link functions for the quantiles of several 1--parameter continuous
  distributions.

}
\description{
  Computes the \code{normal1sdQlink} transformation for the Univariate
  Normal Distribution, its inverse and the first two derivatives.

}
\usage{
   normal1sdQlink(theta,  mean = stop("Please, enter the fixed 'mean'."),
                  p = stop(" Please, enter argument 'p'."),
                  bvalue = NULL, inverse = FALSE, deriv = 0,
                  short = TRUE, tag = FALSE)
}

\arguments{
  \item{theta}{
    Numeric or character. 
    This is \eqn{\theta} by default although it could be \eqn{\eta}
    depending upon other arguments. See below for further details.
    
  }
  \item{mean}{
  A numeric vector or a matrix. It is the (known) fixed mean of the 
  Normal distribution of interest. See below for further details.
  
  }
  \item{p}{
  Numeric vector of \eqn{p}--quantiles to be modelled by this link
  function.
  
  }
  \item{bvalue, inverse, deriv, short, tag}{
  Details at \code{\link[VGAM:Links]{Links}}.
  
  }
}
\details{

  This link function is necessarily required by 
  \code{\link[VGAMextra:normal1sd]{normal1sd}} if quantile regression 
  is to be performed. It computes the \code{normal1sdQlink} 
  transformation, defined as
  %
  \deqn{ {\tt{mean}} + \sqrt{2}\sigma \cdot 
            {\tt{erf^{-1}}}(2 {\tt{p}} - 1).}{
            mean + sqrt(2 \sigma^2)  erf^{-1}(2p - 1).
            }
  %
  Here, \eqn{{\tt{erf^{-1}}}}{erf^(-1)} denotes the inverse 
  of the error function \code{\link[VGAM:erf]{erf}}, and 
  \eqn{\sigma} is the standard deviation (\code{theta}) as in
  \code{\link[VGAMextra:normal1sd]{normal1sd}}.
  Technically, \code{normal1sdQlink} can be used for quantile regression
  over any vector of \eqn{p}--quantiles of Normally distributed data
  with \emph{known} mean.
  
  See \code{\link[VGAMextra:normal1sd]{normal1sd}} for further details 
  about the latter.
  %As a link for the quantiles of the Normal distribution with 
  %\emph{known} mean, the argument \code{theta} refers to its standard 
  %deviation, \eqn{\sigma}, which must be positive.

  
  Values of \eqn{p} out of the open interval \eqn{(0, 1)} or non--positive
  values of \code{theta} will result in \code{Inf}, \code{-Inf},
  \code{NA} or \code{NaN}. 
  

}
\value{

 When \code{deriv = 0}, the \code{normal1sdQlink} transformation
 of \code{theta}, if \code{inverse = FALSE}. Conversely, if
 \code{inverse = TRUE}, \code{theta} becomes \eqn{\eta} and
 the inverse transformation given by 
 \code{(theta -  mean)/}\eqn{\sqrt{2}}{sqrt(2)}
 \code{erf}\eqn{^{-1}}{^(-1)}\eqn{(2}\code{p}\eqn{ - 1)}
 is returned.
 
 For \code{deriv = 1}, \eqn{d} \code{eta} / \eqn{d} \code{theta}
 if \code{inverse = FALSE}. Else, this function returns  
 \eqn{d} \code{theta} / \eqn{d} \code{eta} as a function of \code{theta}.
 
 For \code{deriv = 2}, the second order derivatives are accordingly
 returned.
 
}
\author{
 V. Miranda
 
}
\section{Warning}{
  If \code{p} is a vector, then the recycling rule applies \emph{only} 
  if \code{theta} is entered as a matrix. Else, only the first entry in 
  \code{p} is considered.

}
\note{
  When \code{inverse = TRUE}, the reciprocal of the error function,
  \code{\link[VGAM:erf]{erf}}, evaluated at \eqn{2}\code{p}\eqn{- 1} 
  is required. However, the result is \code{Inf} for \code{p}\eqn{= 0.5}. 
  Here, in consequence, the limit of \code{\link[VGAM:erf]{erf}} when
  \code{p} tends to \eqn{0.5} is returned to avoid numerical issues.

}
\seealso{
  \code{\link[VGAMextra:normal1sd]{normal1sd}},
  \code{\link[VGAM:erf]{erf}},
  \code{\link[VGAM:Links]{Links}}.

}
\examples{
    ###   Example 1  ###
    theta <- seq(0, 3, by = 0.1)[-1]  # This is sigma, then must be positive.
    mean  <- -2.5                     # Intentionally, a negative value for mu.
    p     <- 0.25                     # Modelling the first quartile.
    
    eta <- normal1sdQlink(theta = theta, p = p, mean = mean) 
    inv.eta <- normal1sdQlink(theta = eta, p = p, mean = mean, inverse = TRUE)
    summary(inv.eta - theta)          ## Should be 0

    
    ###  Example 2. Special values of theta, using argument 'bvalue'.  ###
    
    theta <- c(-Inf, -5, -2.5, 0, 2.5, 5, Inf, NA, NaN)
    my.matrix <- rbind(theta, normal1sdQlink(theta = theta, p = p, mean = mean),
                     normal1sdQlink(theta = theta, p = p, mean = mean, bvalue = 1e-5))
    rownames(my.matrix) <- c("theta", "No 'bvalue'", "With 'bvalue'")
    colnames(my.matrix) <- rep("", 9)
    
    my.matrix   # Second row has NAs, whilst third row has NO NAs except for theta = NA
}
