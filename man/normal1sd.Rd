\name{normal1sdff}
\alias{normal1sdff}
\alias{uninormalsd}

\title{
   Estimation and Inference for Conditional Quantiles of
   a 1--parameter Univariate Normal Distribution.

}
\description{
  Maximum likelihood estimation of the standard deviation, including 
  inference for conditional quantiles, of a univariate normal 
  distribution.

}
\usage{
            normal1sdff(zero = NULL, link = "loge",
                        fixed.mean = 0, p.quant = NULL,
                        var.arg = FALSE)
}
\arguments{
  
  \item{zero}{
  Allows to model the single linear predictor in this family function
  as intercept--only. See below for important details about this.
  
  
  }
  \item{link}{
  This is the link function applied to the standard deviation.
  If \code{var.arg} is \code{TRUE}, then \code{link} is applied
  to the variance.
  The default is \code{\link[VGAM:loge]{loge}}. For 
  inference on conditional quantiles entered at \code{p.quant},
  however, it must be manually changed to
  \code{\link[VGAMextra:normal1sdQlink]{normal1sdQlink}}.
  See below for further details.
  
  
  }
  \item{fixed.mean}{
  Numeric, a vector or a matrix. It allocates the (fixed) mean of the 
  response in the fitting process. See below for further details.
  

  }
  \item{p.quant}{
   Numeric. A prototype vector of probabilities indicating the quantiles 
   of interest, when quantile regression is to be performed.
  
  }
  \item{var.arg}{
  If \code{TRUE}, then the variance is estimated, else the
  standard deviation is used.
  
  }
  
}
\details{

  This family function is a variant of 
  \code{\link[VGAM:uninormal]{uninormal}} 
  to estimate the standard deviation of a Normal distribution 
  with \emph{known} mean. The estimated values are returned as
  the fitted values, unlike some other family functions where the mean
  is returned as \emph{fitted values}. However, here the mean is 
  assumed to be known.
  
  By default, the response is supposedly \emph{centered} on its mean, 
  that is, \code{fixed.mean}\eqn{ = 0}. Change this accordingly:
  For a single response or multiple responses, \code{fixed.mean} must
  be a numeric vector where each entry is the mean of each response, 
  only \emph{if} the mean is \emph{fixed}. When the mean is not constant, 
  \code{fixed.mean} must be matrix with the number of columns matching 
  the number of responses.
  
  

  \emph{Quantile regression:}
  The (single) linear/additive predictor by default is the \code{log} 
  of the standard deviation. However, if quantile regression is of 
  primary interest, then the response must be entered using the function
  \code{\link[VGAMextra:Q.reg]{Q.reg}}, and the corresponding 
  \eqn{p}--quantiles through \code{p.quant} in the
  \code{\link[VGAM:vglm]{vglm}} or \code{\link[VGAM:vgam]{vgam}} call.
  Additionally, set \code{\link[VGAMextra:normal1sdQlink]{normalsdQlink}} 
  as the link function via the argument \code{link}. 
  
  This family \pkg{VGAM} function handles multiple responses.
  
}
\value{
  An object of class \code{"vglmff"}.
  See \code{\link[VGAM:vglmff-class]{vglmff-class}} for further details.

}
\author{
    V. Miranda.

}
\section{Warning}{
   Be aware of the argument \code{\link[VGAM:zero]{zero}}:
   by default, the single linear/additive predictor in this family 
   function, say \eqn{\eta}, can be modeled in terms of covariates, 
   i.e., \code{zero = NULL}. 
   To model \eqn{\eta}  as intercept--only, set \code{zero = "sd"}.
   
   See \code{\link[VGAM:zero]{zero}} for more details about this.
  
}
\seealso{
   \code{\link[VGAMextra:normal1sdQlink]{normal1sdQlink}},
   \code{\link[VGAM:loge]{loge}},
   \code{\link[VGAM:uninormal]{uninormal}},
   \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
   \code{\link[VGAM:zero]{zero}},
   \code{\link[VGAM:vgam]{vgam}},
   \code{\link[VGAM:vglm]{vglm}}.
   

}
\examples{
   set.seed(121216)
   my.mean <- -1      #  Mean (CONSTANT)
   my.sd   <- 2.5
   y <- rnorm(100, mean = my.mean, sd = 2.0)      # Generate some data.
   normdat <- data.frame(y = y)                   # Setting up our data.
 
   \donttest{ 
   # Plotting the data
     plot(y, main = c("Y ~ Normal ( mean(known), sd = 2.5 ). "),
          ylab = "The data", pch = 20, 
          xlim = c(0, 100), ylim = c(-7, 7), col = "blue")
     abline(h = 0, v = 0, lwd = 2, col = "black")
   }

   ### EXAMPLE 1. Estimate the SD with two responses. The mean is fixed. ###
   
   fit1 <- vglm(cbind(y, y) ~ 1, family = normal1sdff(fixed.mean = my.mean), 
               data = normdat, trace = TRUE, crit = "coef")
   Coef(fit1) 
   summary(fit1)
    
    
   ### EXAMPLE 2. Quantile regression. The link normal1sdQlink() is used. ###
  
   my.p <- c(25, 50, 75) / 100  # Quantiles 25\%, 50\% and 75\% are of interest.
   
   fit2 <- vglm(Q.reg(y, length.arg = 3) ~ 1, 
                family = normal1sdff(fixed.mean = my.mean, p.quant = my.p,
                                   link = normal1sdQlink), 
                data = normdat, trace = TRUE, crit = "coef")
    summary(fit2)
    head(predict(fit2))
    constraints(fit2)


   ### EXAMPLE 3. Complete the plot. Quantiles matching. ###
   
   \donttest{
   ( my.c3Q <- coef(fit2, matrix = TRUE) )
   with(normdat, lines(rep(my.c3Q[1], 100), col = "tan"   , lty = "dotted", lwd = 2))
   with(normdat, lines(rep(my.c3Q[2], 100), col = "orange", lty = "dotted", lwd = 2))
   with(normdat, lines(rep(my.c3Q[3], 100), col = "brown1", lty = "dotted", lwd = 2))
   legend(20, 7.0, c("Percentil 75", "Percentil 50", "Percentil 25"),
          col = c("brown1", "orange", "tan"),
          lty = rep("dotted", 3), lwd = rep(2, 3), cex = 0.75)
   }
}
