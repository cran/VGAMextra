\name{gammaRff}
\alias{gammaRff}
%\alias{gammaRMeanff}
%\alias{gammaRMean}

\title{
 
  2--parameter Gamma Distribution
    
}
\description{
   %
   Estimates the 2--parameter gamma distribution by maximum likelihood.
   One linear predictor models the mean.
   %
}
\usage{
        gammaRff(zero = "shape", lmu = "gammaRMlink",
                 lrate = NULL, lshape = "loglink",
                  irate = NULL,   ishape = NULL, lss = TRUE)
                       
}
\arguments{
    \item{zero}{
     Specifies the parameters to be modelled as intercept--only.
     
     See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
    
    }
    \item{lmu}{
    The link function applied to the gamma distribution mean, i.e.,
    \code{\link[VGAMextra:gammaRMlink]{gammaRMlink}}.
  
    
    }
    \item{lrate, lshape, irate, ishape, lss}{
    Same as \code{\link[VGAM:gammaR]{gammaR}}.
    
    }
}
\details{
  
  This family function slightly enlarges the functionalities of 
  \code{\link[VGAM:gammaR]{gammaR}} by directly modelling the mean
  of the gamma distribution. It performs very much like
  \code{\link[VGAM:gamma2]{gamma2}}, but involves the ordinary 
  (not reparametrized) density, given by
  %
  \deqn{f(y; \alpha, \beta) = \frac{ \beta^\alpha }{ \Gamma(\alpha) } 
         e^{-\beta y} y^{\alpha - 1},}{
         f(y; alpha, beta) = rate^(shape)  exp(-rate * y) 
                y^(shape - 1) / \Gamma(shape),
         }
  %
  Here, \eqn{\alpha}{'shape'} and \eqn{\beta}{'rate'} are positive shape
  and rate parameters as in \code{\link[VGAM:gammaR]{gammaR}}.
  %The mean is
  %\eqn{\mu = \alpha / \beta}{ \mu = shape / rate}.
  The default linear predictors are 
  \eqn{\eta1 = {\tt{gammaRMlink}}(\alpha; \beta) = \log \mu =
  \log (\alpha / \beta)}{
  \eta1 = gammaRMlink(rate; shape) = log \mu =
  log (rate / shape)}, and
  \eqn{\eta2 = \log \alpha}{\eta2 = log (shape)},
  unlike \eqn{\eta1 = \log \beta}{\eta1 = log (rate)} and
  \eqn{\eta2 = \log \alpha}{\eta2 = log (shape)}
  from \code{\link[VGAM:gammaR]{gammaR}}.
  %the default linear predictors are
  %\eqn{\eta1 = \log \alpha}{\eta1 = log (shape)} and
  %\eqn{\eta2 = \log \beta}{\eta2 = log (rate)}, this family function
  %re--defines this structure by setting up
  % by default, where
  %\code{\link[VGAMextra:gammaRMlink]{gammaRMlink}} is 
  %the link function for the mean of \eqn{Y}.
  
  \code{lmu} overrides \code{lrate} and no link other than
  \code{\link[VGAMextra:gammaRMlink]{gammaRMlink}} is a valid entry
  ({\code{lmu}}). To
  mimic \code{\link[VGAM:gammaR]{gammaR}} simply
  set \code{lmu = NULL} and \code{lrate = "loglink"}.
  The mean (\eqn{\mu}) is returned as the fitted values.
  
  %Particularly, \code{lss} works exactly as in
  %\code{\link[VGAM:gammaR]{gammaR}}.
  %For further choices on link functions for \eqn{\eta_2}{\eta2},
  %see \code{\link[VGAM:Links]{Links}}.  

  \code{gammaRff} differs from
  \code{\link[VGAM:gamma2]{gamma2}}. The latter estimates a
  re-parametrization
  of the gamma distribution in terms \eqn{\mu}{\mu} and
  \eqn{\alpha}{shape}. 
    This \bold{VGAM} family function does not handle censored data.
  % viz. the density is re--expressed
%  in terms of \eqn{\mu}{\mu} and \eqn{\alpha}{shape}.
  %To he link
  %must be necessarily used via \code{lmu}, which overrides
  %overrides \code{lrate}. 
  %Then, the transformed mean
  %\deqn{{\tt{gammaRMlink}}(\alpha; \beta)}{gammaRMlink(shape; rate)}
  %are returned as the fitted values, for estimated
  %\eqn{\alpha}{shape} and \eqn{\beta}{rate}.
  
  
}
\value{
  An object of class \code{"vglm"}.
  See \code{\link[VGAM:vglm-class]{vglm-class}} for full details.
  
}
\section{Note}{
  
  %The fitted values returned are \code{gammaRMlink}--transformed,
  %provided the mean is modelled via \code{lmu}
  
  The parameters \eqn{\alpha}{shape} and \eqn{\beta}{rate}
  match the arguments \code{shape} and \eqn{rate} of
  \code{\link[stats:rgamma]{rgamma}}.
  
  Multiple responses are handled.
  
}
%\section{Warning}{
% 
% Conforming with the 2--parameter gamma distribution above, 
% both, the \eqn{\alpha}{shape} and \eqn{\beta}{rate} parameters
% must be positive. However, in order for the 
% \code{\link[VGAMextra:gammaRMlink]{gammaRMlink}} link
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
    \code{\link[VGAMextra:gammaRMlink]{gammaRMlink}},
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
    mu <- exp(2 + 0.5 * x2)
  
  # Shape and rate parameters in terms of 'mu'
    shape <- rep(exp(1), nn)
    rate <- gammaRMlink(theta = log(mu), shape = shape,
                            inverse = TRUE, deriv = 0)
  
  # Generating some random data
    y1 <- rgamma(n = nn, shape = shape, rate =  rate)
    gdata <- data.frame(x2 = x2, y1 = y1)
    rm(y1)

  # lmu = "gammaRMlink" replaces lshape, whilst lrate = "loglink"
    fit1 <- vglm(cbind(y1, y1) ~ x2,
                 gammaRff(lmu = "gammaRMlink", lss = TRUE, zero = "shape"),
                 data = gdata, trace = TRUE, crit = "log")
     coef(fit1, matrix = TRUE)
     summary(fit1)
    
  # Comparing fitted values with true values.
    compare1 <- cbind(fitted.values(fit1)[, 1, drop = FALSE], mu)
    colnames(compare1) <- c("Fitted.vM1", "mu")
    head(compare1)
 
  
  ### Mimicking gammaR. Note that lmu = NULL.
    fit2 <- vglm(y1 ~ x2, gammaRff(lmu = NULL, lrate = "loglink",
                            lshape = "loglink", lss = FALSE, zero = "shape"),
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

  ### Use gamma2()
     fit3 <- vglm(y1 ~ x2, gamma2,
                 data = gdata, trace = TRUE, crit = "log")
    fit1.fit3 <- cbind(fitted.values(fit1)[, 1, drop = FALSE], 
                        fitted.values(fit2), fitted.values(fit3))
    colnames(fit1.fit3) <- c("Fitted.vM1", "Fitted.vM2", "Fitted.vM3")
    head(fit1.fit3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
