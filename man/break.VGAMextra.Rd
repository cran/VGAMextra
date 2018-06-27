\name{break.VGAMextra}
\alias{break.VGAMextra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Names/Value of linear predictors/parameters in time series 
family functions. }
\description{

  Splitting out the names of linear predictors or Numeric values for 
  parameters in time series family functions in \pkg{VGAMextra}.

}
\usage{
      break.VGAMextra(eta      = NULL,
                      M1       = NULL,
                      noInter  = NULL,
                      bOrder   = NULL,
                      NOS      = NULL,
                      lInter   = "identitylink",
                      lvar     = "loge",
                      lsd      = "loge",
                      lcoeff1  = "rhobit",
                      lcoeff2  = "rhobit",
                      typeTS   = "AR",
                      namesLP  = FALSE,
                      Complete = FALSE,
                      varArg   = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{eta}{
  A matrix of dimensions \code{c(n, M)} storing the linear predictors 
  values coming from the \code{vglm} fit. Here, \eqn{M} is the number
  of parameters. See \code{warning} below for further information.
  %By default \code{order = NULL}. Thus, a valid value must be entered
  %initially.
  
  }
  \item{M1}{
  Number of parameters involved in the \code{vglm} fit.
  
  }
  \item{noInter}{
   Logical. To determine whether the intercept is 
   estimated. If \code{'TRUE'},
   the intercept is not estimated and set to 0.
  
  }
  \item{bOrder}{
  A vector. 
  The order of the linear process fitted. Either a single number (if one 
  response), or a vector (if multiple responses).
  
  }
  \item{NOS}{
  Integer. Number of respones set in the \code{vglm} call.
  
  }
  \item{lInter, lvar, lsd, lcoeff1, lcoeff2}{
  Link functions applied to parameters. Same as in 
  \code{\link[VGAMextra:ARXff]{ARXff}}, or 
  \code{\link[VGAMextra:MAXff]{MAXff}}.
  
  }
  \item{typeTS}{
  Character. Currently, options \code{"AR"} for Autoregressive, and 
  \code{"MA"} for Moving Average processes are handled.
  
  }
  \item{namesLP}{
  Logical. This function returns either the names of linear the
  predictors/parameters ( if \code{namesLP = TRUE} ) or parameter
  values (default) broken down from the \code{eta} matrix.
  
  }
  \item{Complete}{
  Logical. If \code{TRUE}, columns of zeros are incorporated into the
  matrix \code{eta}. See below for further details.
  
  }
  \item{varArg}{
  Sames as in \code{\link[VGAMextra:ARXff]{ARXff}} or 
  \code{\link[VGAMextra:MAXff]{MAXff}}
  
  }
}

\details{

Time series family functions in \pkg{VGAMextra} currently recycle 
the order set in the \code{\link[VGAM:vglm]{vglm}}. 
Particularly, it occurs when the number of responses is fewer than 
the specified order. For instance, if the order set in \code{vglm}
is \eqn{c(1, 3)}, and 5 responses are managed, then the new order
becomes \eqn{c(1, 3, 1, 3, 1)}. 

Due to such flexibility, time series family functions require 
specific functions to unload the amount of code within each one. 

Moreover, when the order is recycled, the matrix \code{eta} is 
\emph{completed}, as if the order was the same for each response.
This feature is enabled when \code{Complete = TRUE}.
This `common' order turns out to be the maximum order established in
the vector \code{order}. This trick makes the family function to work
properly. To return to the riginal `order', \code{eta} is reduced 
in the same number of colums initially added.
%by calling the 
%function \code{\link[VGAMextra:reduceCol.VGAMextra]{reduceCol.VGAMextra}}.

\code{break.VGAMextra} works in this context. It may return either 
the names of the linear predictors/parameters, or the parameter values 
splitted out as a list. Thus, link functions entered in the \code{vglm} 
call must be passed down to this functions.
For further details on link functions refer to 
\code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.



}

\value{

A list containing either the names of the linear predictors or the 
parameters values (not linear predictors) unwrapped from tje
\code{eta} matrix, as follows:

a) If \code{namesLP = FALSE} (default), value of parameters are 
returned in this order: the intercept (1), standard deviation and
variance of the white noise (2, 3), and the coefficients (4).

b) If \code{namesLP = TRUE}, \bold{names of linear predictors} are
returned in the first entry, whereas \bold{parameter names} are
allocated to the second entry.

Yee and Wild (1996) provide more detailed information about 
the relationship between linear predictors and parameters within 
the VGLM statistical framework.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{


Yee, T. W. and Wild, C. J. (1996)
Vector Generalized Additive Models.
\emph{Journal of the Royal Statistical Society}, 
Series B, Methodological, \bold{58(3)}, 481--493.

}
\section{Warning}{

Note that library \pkg{VGAM} is definitely required.

}

\author{

 Victor Miranda and T. W. Yee

}

\section{Warning}{

Be aware of the dimensions of matrix \code{eta}. It is \code{c(n, M)},
where \eqn{n} is the sample size, and \eqn{M} is the number of 
parameters. If multiple responses, then \eqn{M} equals the summation
of parameters individually.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

  %\code{\link[VGAMextra:pre2.wz]{pre2.wz}},
  \code{\link[VGAMextra:ARXff]{ARXff}},
  \code{\link[VGAMextra:MAXff]{MAXff}},
  %\code{\link[VGAMextra:reduceCol.VGAMextra]{reduceCol.VGAMextra}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[VGAM:vglm]{vglm}}.
  %\code{\link[VGAMextra:cm.ARMA]{cm.ARMA}},
  %\code{\link[VGAM:vglm.control]{vglm.control}},
  %\code{\link[VGAMextra:MAXff]{MAXff}},
  %\code{\link[VGAMextra:cm.ARMA]{cm.ARMA}},


}

\examples{
\donttest{

library(VGAM)

eta     <- matrix(runif(100), nrow = 10, ncol = 10)
M1      <- c(5, 5)
noInter <- FALSE
bOrder  <- c(3, 3)
NOS     <- 2

### ONLY LINEAR PREDICTORS/PARAMETERS NAMES!
### RETURNED OBJECT IS A LIST !

break.VGAMextra(M1      = M1, 
                noInter = noInter, 
                bOrder  = bOrder, 
                NOS     = NOS, 
                typeTS  = "AR", 
                namesLP = TRUE, 
                varArg  = TRUE)

### PARAMETER VALUEs... "UNWRAPPED". Inverse link functions are applied.
###  Note that namesLP must be set to FALSE

break.VGAMextra(eta     = eta,
                M1      = M1, 
                noInter = noInter, 
                bOrder  = bOrder, 
                NOS     = NOS, 
                typeTS  = "AR", 
                namesLP = FALSE, 
                varArg  = TRUE)



}

}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
