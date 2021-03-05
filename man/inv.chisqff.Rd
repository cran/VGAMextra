\name{inv.chisqff}
\alias{inv.chisqff}

\title{
   Inverse Chi--squared Distribution.

}
\description{
  Maximum likelihood estimation of the degrees of freedom for an
   inverse chi--squared distribution using Fisher scoring.

}
\usage{
              inv.chisqff(link = "loglink", zero = NULL)
}
\arguments{
  \item{link, zero}{
  \code{link} is the link function applied to the degrees of freedom,
  leading to the unique linear predictor in this family function.
  By default, the link is \code{\link[VGAM:loglink]{loglink}}.
   
  \code{zero}  allows to model the single linear predictor as
   intercept--only.  
   
   For further details, see 
   \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  
  
  }
}
\details{
  
  The inverse chi--squared distribution 
  with \eqn{df = \nu \geq 0}{df = \nu \ge 0} degrees of 
  freedom implemented here has density
  
   %
  \deqn{f(x; \nu) = \frac{ 2^{-\nu / 2} x^{-\nu/2 - 1} 
               e^{-1 / (2x)} }{ \Gamma(\nu / 2) }, }{
        f(x; \nu) = 2^(-\nu / 2) x^(-\nu/2 - 1) e^(-1 / (2x)) /
                 \Gamma(\nu / 2),
               }
  %
  where \eqn{x > 0}{x > 0}, and
  \eqn{\Gamma} is the \code{\link[base:gamma]{gamma}} function.
  The mean of \eqn{Y} is \eqn{1 / (\nu - 2)} (returned as the fitted
  values), provided \eqn{\nu > 2}.
  
  That is, while the expected information matrices used here are 
  valid in all regions of the parameter space, the regularity conditions 
  for maximum likelihood estimation are satisfied only if \eqn{\nu > 2}.
  To enforce this condition, choose 
  \code{link = logoff(offset = -2)}.
  
  As with, \code{\link[VGAM:chisq]{chisq}}, the degrees of freedom are
  treated as a parameter to be estimated using (by default) the
  link \code{\link[VGAM:loglink]{loglink}}. However, the mean can also
  be modelled with this family function. 
  See \code{\link[VGAMextra:inv.chisqMlink]{inv.chisqMlink}} 
  for specific details about this.
  
  This family \pkg{VGAM} function handles multiple responses.
  
}
\value{
  An object of class \code{"vglmff"}.
  See \code{\link[VGAM:vglmff-class]{vglmff-class}} for further details.

}
\author{
    V. Miranda.

}
\note{

  As with \code{\link[VGAM:chisq]{chisq}} or 
  \code{\link[stats:Chisquare]{Chisquare}}, the degrees of freedom are
  non--negative but allowed to be non--integer.
    
}
\section{Warning}{
   By default, the single linear/additive predictor in this family 
   function, say \eqn{\eta = \log dof}{\eta = log (dof)}, 
   can be modeled in terms of covariates, 
   i.e., \code{zero = NULL}. 
   To model \eqn{\eta} as intercept--only set \code{zero = "dof"}.
   
   See \code{\link[VGAM:zero]{zero}} for more details about this.
  
}
\seealso{
   \code{\link[VGAM:loglink]{loglink}},
   \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
   \code{\link[VGAMextra:inv.chisqMlink]{inv.chisqMlink}},
   \code{\link[VGAM:zero]{zero}}.
   

}
\examples{
   set.seed(17010504)
   dof   <- 2.5 
   yy    <- rinv.chisq(100, df = dof)     
   ics.d <- data.frame(y = yy)             # The data.
 
   
   fit.inv <- vglm(cbind(y, y) ~ 1, inv.chisqff, 
                   data = ics.d, trace = TRUE, crit = "coef")
   Coef(fit.inv) 
   summary(fit.inv)
    
}
