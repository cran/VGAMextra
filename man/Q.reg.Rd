\name{Q.reg}
\alias{Q.reg}
\title{
   %
   Conditional quantile regression with \pkg{VGAM}
   %
}
\description{
    Use this function to adequately confer the \code{formula} 
    in \pkg{VGAM} when fitting quantile regression models.
   %\code{q.reg} must be used to specify the \code{formula} when
   %fiiting quantile regression models.

}
\usage{
    Q.reg(y, pvector = NULL, length.arg = NULL)
}
\arguments{
  \item{y}{
   Numeric, a vector or a matrix. It is the response or dependent 
   variable in the \code{formula} of the model to be fit, as in 
   \code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}}.
   See below for further details.
  
  }
  \item{pvector}{
  A prototype vector. Entries are the conditional \eqn{p}--quantiles
  in the fitting process.
  
  }
  \item{length.arg}{
  A length--1 positive integer. It is the number of \eqn{p}--quantiles
  to be modelled.
  
  }
}
\details{
    
    Conditional quantile regression can be carried out using family
    functions in \pkg{VGAM} and \pkg{VGAMextra}.
    The \code{formula} must be set up using this function, \code{Q.reg}.
    Here, the \eqn{p}--quantiles of interest may be entered via
    \code{pvector}. Alternatively, use argument \code{length.arg} by
    establishing the length of \code{pvector}.
    
    Besides, the corresponding link must be entered.
    For example, \code{\link[VGAMextra:gamma1Qlink]{gamma1Qlink}}
    is the proper link to fit models of conditional quantiles for 
    data distributed as Gamma via the family function
    \code{\link[VGAM:gamma1]{gamma1}}.
    
    See examples for further details.
    
}
\value{
   A matrix, each column adequately arranged for regression on conditional
   quantiles, conforming with \pkg{VGAM}. 
   
   Indeed, this is equivalent to \code{cbind(y, y, ...)}, where the 
   total number of columns is, either the length of \code{pvector}, or
   \code{length.arg}.
   
}
\author{
  V. Miranda and T. W. Yee.
  
}
\note{
  Link functions for quantile regression in \pkg{VGAM} require the
  vector of \eqn{p}--quantiles of interest via the argument \code{p}.
  See \code{\link[VGAMextra:normal1sdQlink]{normal1sdQlink}} or
  \code{\link[VGAMextra:maxwellQlink]{maxwellQlink}} for instance.
  
  Therefore, the integer entered at \code{length.arg} in this function, 
  if utilized, must match the length of the vector \code{p}. Else,
  it will be recycled.
 
}
\seealso{
  \code{\link[VGAMextra:normal1sdQlink]{normal1sdQlink}},
  \code{\link[VGAMextra:maxwellQlink]{maxwellQlink}},
  \code{\link[VGAMextra:gamma1Qlink]{gamma1Qlink}},
  \code{\link[VGAM:gamma1]{gamma1}},
  \code{\link[VGAM:vglm]{vglm}},
  \code{\link[VGAM:vgam]{vgam}}

}
\examples{


  ###  Quantile regression with data distributed as Maxwell(s)  ###
   set.seed(12073)
   x2 <- seq(0, 100,length.out = 100)       # independent variable
   b0 <- 0.5                                # true intercept
   b1 <- 0.25                               # true slope
   b2 <- 0.02                               # true second order coef.
   alpha <- b0 + b1 * x2 + b2 * x2^2        # Quadratically modelling the parameters
   nn <- 100                                # Sample size
  
   # The data as a data frame. #
   mdata <- data.frame(y = rmaxwell(n = nn, rate = alpha), x2 = x2, x3 = x2^2)
   
   # Quantile regression using our link function maxwellQlink(). #
   # Quantiles 25\%, 50\%, 75\% are of interest #
   my.p <- c(0.25, 0.50, 0.75)
   

   fit <- vglm(Q.reg(y, pvector = my.p) ~ x2 + x3, 
   
  # OPTIONALLY Q.reg(y, length = length(my.p)) ~ x2 + x3
   
              maxwell(link = "maxwellQlink(p = my.p)"), 
              data = mdata,  trace = TRUE, crit = "coef")

   coef(fit, matrix = TRUE) 
   summary(fit)
   head(predict(fit))
   constraints(fit)
   
}
