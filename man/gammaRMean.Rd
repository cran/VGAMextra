\name{gammaRMean}
\alias{gammaRMean}

\title{
 
  2--parameter Gamma Distribution
    
}
\description{
   %
   Estimates the 2--parameter gamma distribution by maximum likelihood
   estimation by modelling its mean.
   %
}
\usage{
        gammaRMean(zero = "rate", lmu = "gammaRMeanlink",
                   lrate = "loge", lshape = NULL,
                   irate = NULL,   ishape = NULL,
                   lss = TRUE)
                       
}
\arguments{
    \item{zero}{
     Specifies the parameters to be modelled as intercept--only.
     
     See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}
     for further information.
    
    }
    \item{lmu}{
    The link function applied to model the mean of this distribution:
    \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}}.
  
    
    }
    \item{lrate, lshape, irate, ishape, lss}{
    Same as \code{\link[VGAM:gammaR]{gammaR}}.
    
    }
}
\details{
  
  This family function slightly enlarges the functionalities of 
  \code{\link[VGAM:gammaR]{gammaR}} allowing the extent of the 
  linear predictor to model the mean of the 2--parameter
  gamma distribution, with density
  %
  \deqn{f(y; \alpha, \beta) = \frac{ \beta^\alpha }{ \Gamma(\alpha) } 
         e^{-\beta y} y^{\alpha - 1},}{
         f(y; alpha, beta) = rate^(shape)  exp(-rate * y) 
                y^(shape - 1) / \Gamma(shape),
         }
  %
  Here, \eqn{\alpha}{shape} and \eqn{\beta}{rate} are positive shape
  and rate parameters as in \code{\link[VGAM:gammaR]{gammaR}}.
  The mean is given by
  \eqn{\mu = \alpha / \beta}{ \mu = shape / rate}.
  
  
  Unlike \code{\link[VGAM:gammaR]{gammaR}} where
  the default linear predictors are
  \eqn{\eta1 = \log \alpha}{\eta1 = log (shape)} and
  \eqn{\eta2 = \log \beta}{\eta2 = log (rate)}, this family function
  re--defines this structure by setting up
  \eqn{\eta1 = {\tt{gammaRMeanlink}}(\alpha; \beta)}{
  \eta1 = gammaRMeanlink(shape; rate)} and
  \eqn{\eta2 = \log \beta}{\eta2 = log (rate)} by default, where
  \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} is 
  the link function for the mean of \eqn{Y}.

  
  To mimic the work of \code{\link[VGAM:gammaR]{gammaR}},
  set \code{lmu = NULL} and \code{lshape = "loge"}.
  Particularly, \code{lss} works exactly as with
  \code{\link[VGAM:gammaR]{gammaR}}
  


  For further choices on link functions for \eqn{\eta_2}{\eta2},
  see \code{\link[VGAM:Links]{Links}}.  

  This family function also differs from
  \code{\link[VGAM:gamma2]{gamma2}}. The latter is re-parametrization
  of the gamma distribution to estimate \eqn{\mu}{\mu} and
  \eqn{\alpha}{shape}. That is, the density is re--expressed
  in terms of \eqn{\mu}{\mu} and \eqn{\alpha}{shape} specifically.
  

  Notice, to model the mean of the gamma distribution with this 
  family function, the link
  \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} 
  must be necessarily used via \code{lmu}. Here, \code{lmu}
  overrides the work of \code{lshape}. Then, the transformed mean
  \deqn{{\tt{gammaRMeanlink}}(\alpha; \beta)}{gammaRMeanlink(shape; rate)}
  are returned as the fitted values, for estimated
  \eqn{\alpha}{shape} and \eqn{\beta}{rate}.
  
  
}
\value{
  An object of class \code{"vglm"}.
  See \code{\link[VGAM:vglm-class]{vglm-class}} for full details.
  
}
\section{Note}{
  
  The fitted values returned are \code{gammaRMeanlink}--transformed,
  provided the mean is modelled via \code{lmu}
  
  The parameters \eqn{\alpha}{shape} and \eqn{\beta}{rate}
  match with the arguments \code{shape} and \eqn{rate} of
  \code{\link[stats:rgamma]{rgamma}}.
  
  Multiple responses are handled.
  
}
%\section{Warning}{
% 
% Conforming with the 2--parameter gamma distribution above, 
% both, the \eqn{\alpha}{shape} and \eqn{\beta}{rate} parameters
% must be positive. However, in order for the 
% \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} link
% to be positive, the ratio
 %
% \deqn{ \alpha_i / \beta_i,~~i = 1, \ldots, n}{
%         shape_i / rate_i,  i = 1, \ldots, n}

% must be strictly greater than 1.0.
% Failure on metting this requirement at any IRLS iteration 
% including initial values, will cause numerical instability and,
% eventually, the algorithm will stop.
%
%}
\author{
  V. Miranda and Thomas W. Yee.

}
\references{
Yee, T. W. (2015)
\emph{Vector Generalized Linear and Additive Models: 
With an Implementation in R}. 
Springer, New York, USA.

}
\seealso{
    \code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}},
    \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
    \code{\link[VGAM:gammaR]{gammaR}},
    \code{\link[VGAM:gamma2]{gamma2}},
    \code{\link[VGAM:Links]{Links}}.
     
}
\examples{
 
  ### Modelling the mean in terms of x2, two responses.
  
    set.seed(2017022101)
    nn <- 80
    x2 <- runif(nn)
    mu <- exp(1 + 0.5 * x2)
  
  # Shape and rate parameters in terms of 'mu'
    rate <- rep(exp(1), nn)
    shape <- gammaRMeanlink(theta = mu, rate = rate,
                            inverse = TRUE, deriv = 0)
  
  # Generate some random data
    y1 <- rgamma(n = nn, shape = shape, rate =  rate)
    gdata <- data.frame(x2 = x2, y1 = y1)
    rm(y1)

  # lmu = "gammaRMeanlink" replaces lshape, whilst lrate = "loge"
    fit1 <- vglm(cbind(y1, y1) ~ x2,
                 gammaRMean(lmu = "gammaRMeanlink", lss = FALSE, zero = "rate"),
                 data = gdata, trace = TRUE, crit = "log")
     coef(fit1, matrix = TRUE)
     summary(fit1)
    
  # Compare fitted values with true values.
    compare1 <- cbind(fitted.values(fit1)[, 1, drop = FALSE], mu)
    colnames(compare1) <- c("Fitted.vM1", "mu")
    head(compare1)
 
  
  ### Mimicking gammaR. Notice lmu = NULL.
    fit2 <- vglm(y1 ~ x2, gammaRMean(lmu = NULL, lshape = "loge",
                                     lrate = "loge", lss = FALSE, zero = "rate"),
                 data = gdata, trace = TRUE, crit = "log")
 
  # Compare fitted values with true values.
    compare2 <- with(gdata, cbind(fitted.values(fit2), y1, mu))
    colnames(compare2) <- c("Fitted.vM2", "y", "mu")
    head(compare2)
 
    
  ### Fitted values -- Model1 vs Fitted values -- Model2
    fit1vsfit2 <- cbind(fitted.values(fit1)[, 1, drop = FALSE], 
                        fitted.values(fit2))
    colnames(fit1vsfit2) <- c("Fitted.vM1", "Fitted.vM2")
    head(fit1vsfit2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
