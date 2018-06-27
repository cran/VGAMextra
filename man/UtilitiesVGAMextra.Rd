\name{UtilitiesVGAMextra}
\alias{UtilitiesVGAMextra}
\alias{Is.Numeric}
\alias{is.FormulaAR}
\alias{cross.gammas}
\alias{extract.Residuals}
\alias{fittedVGAMextra}
\alias{weightsVGAMextra}
\alias{WN.lags}
\alias{isNA}
\alias{inspectVGAMextra}
\alias{ffff.help}
\alias{XLMmat}
%\alias{almfact}
%\alias{findFirstMethodTS}
%\alias{poratFIM}
%\alias{inter}
%\alias{reduceCol.VGAMextra}
\title{ Utility Functions for the \pkg{VGAMextra} Package }


\description{

  A set of common utility functions required by time series family 
  functions at 'VGAMextra'.


}
\usage{

  Is.Numeric(x, isInteger  = FALSE, length.arg = NULL, Nnegative  = NULL)
  is.FormulaAR(Model = ~ 1, Resp  = 1) 
  cross.gammas(x, y = NULL, lags = 1)
  WN.lags(y, lags, to.complete = NULL)
  extract.Residuals(object, TSprocess,...)
  fittedVGAMextra(object,...)
  weightsVGAMextra(object, type.w = "prior",...)
  XLMmat(object,...)
  %findFirstMethodTS(methodsts, charvec)
  %almfact(z)
  %poratFIM(x = NULL, m.mean, WNsd, ARcoeff1)

}
\arguments{
\item{x}{
A vector of quantiles. Particularly, for \code{Is.Numeric} it is a
single number (or vector) to be tested: Whether is numeric or not.

}
\item{y}{
Vector of quantiles to be lagged. Then, the \emph{cross - covariances}
are computed from \eqn{x} and \eqn{y_t}{y[t]}, \eqn{x} and 
\eqn{y_{t -1}}{y[t - 1]}, etcetera.

}
\item{isInteger}{
Logical. If \code{TRUE}, it verifies that quantiles \code{x} are integers.
Default is \code{FALSE}.

}
\item{lags}{
Integer indicating the number of lags or \emph{delays} to be applied to 
vector \code{y}. Then, calculate the cross-covariance between the pair
of signals \code{x} and delayed samples computed from \code{y}.

}
\item{length.arg}{
Integer. If \code{length.arg > 0}, it verifies that the length of
\code{x} matches \code{length.arg}.

}
\item{Model}{
Formula. A symbolic form of the models fitted by the \code{vglm} call.
See \code{\link[stats:formula]{formula}} for further details.

}
\item{Nnegative}{
Logical. If \code{TRUE}, it verifies that \code{x} (all entries) are
positive.

}
\item{Resp}{
Integer. The number of \emph{responses} in the \code{Model}. It must 
macth the number of respones entered in the \code{vglm} call.

}
\item{object}{
An object of class \code{'vglm'}.
See \code{\link[VGAM:vglm-class]{vglm-class}} for details.

}
\item{TSprocess}{
Logical, what time series model is being fitted.
Choices are 
\code{'AR', 'MA', 'ARMA'} and  \code{'ARIMA'}.

}
\item{type.w}{
Character. What type of \code{weights} are to be used. 
Default is \code{"prior"}.
These are extracted from the slot \code{@prior.weights} of
\code{object}.

}
\item{to.complete}{
 Use this argument to fill in the first 'p' observations when
 computing the lagged vectors in time series.
 
}
\item{...}{
Additional parameters required by function
\code{\link[VGAMextra:extract.Residuals]{extract.Residuals}}.

}

}
\details{

A set of utility functions in \pkg{VGAMextra} for different
purposes.

Specially for time series family functions in \pkg{VGAMextra} which 
involve specific checks on the majority of arguments entered by 
the user.


}
\value{

 \code{Is.Numeric()} returns a logical vector (or value)
 (\code{TRUE} or \code{FALSE}), after verifying whether quantiles
 \code{x} satisfies all conditions entered.
 
 For \code{is.FormulaAR()}, this function returns a logical value,
 after verifying whether the expression entered for the \code{Model}
 argument in \code{\link[VGAMextra:cm.ARMA]{cm.ARMA}} is an 
 object of class `\code{formula}'. 
 
 Particularly, \code{cross.gammas()} computes either the single lagged 
 covariance(s) from quantiles given in \code{x} or the 
 lagged cross-covariance(s) from values given in \code{x} and \code{y}. 
 
 \code{extract.Residuals()} extracts the residuals of the process from
 slot \code{@residuals}, whilst

 \code{fittedVGAMextra} and \code{weightsVGAMextra} return the
 fitted values and the weights from the \code{vglm} object,
 correspondingly.
 
 \code{isNA} and \code{inspectVGAMextra} are essentially required when
 implementing link functions in \pkg{VGAMextra}.
  
}

%\references{


  %Kleiber, C. and Kotz, S. (2003) 
  %\emph{Statistical Size Distributions in Economics and Actuarial Sciences}.
  %Wiley Series in Probability and Statistics. Hoboken, New Jersey, USA.  
  %Abramowitz, M. and Stegun, I. A. (1972) 
  %\emph{Handbook of Mathematical Functions with Formulas, Graphs, and 
  %Mathematical Tables}, ch.6, p.255.
  %Dover, New York, USA.
  
%}
\author{V. Miranda and T. W. Yee.}
\seealso{

    \code{\link[VGAMextra:cm.ARMA]{cm.ARMA}}.
    %\code{\link[VGAMextra:ARMAff]{ARMAff}}.
  
}
\examples{

# Example 1.
myModel1 <- ~ x1 + x2
is.FormulaAR(myModel1)       # TRUE

test <- list( cbind(y1, y2) ~ x1, ~ x2 - 1)
is.FormulaAR(test)          # FALSE
is.FormulaAR(test[[1]], 2)  # TRUE

# Example 2.

x1 <- c(1:3, 4.5, -Inf)
Is.Numeric(x1)                                        # TRUE
Is.Numeric(x1, length.arg = 5)                        # TRUE
Is.Numeric(x1, length.arg = 5, isInteger = TRUE)      # FALSE
Is.Numeric(x1, length.arg = 5, Nnegative = TRUE)      # FALSE

# Example 3. 
# Here, 'cross.gammas' computes Cov(x, y_{t - 1}), Cov(x, y_{t - 2}) and
# Cov(x, y_{t - 3}). 

x <- runif(50)
y <- runif(50)
cross.gammas(x, y, lags = 3)


}


%\keyword{distribution}
