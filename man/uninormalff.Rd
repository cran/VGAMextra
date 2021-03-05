\name{uninormalff}
\alias{uninormalff}


\title{
 
     Normal (distribution--specified) quantile regression
    
}
\description{
   %
   Distribution--specified quantile regression.
   An extension of \code{uninormal} from \pkg{VGAM}.
   It handles effectively \code{uninormalQlink} via the first
   linear predictor.

   %
}
\usage{
    %
     uninormalff(link1 = "identitylink", lsd = "loglink",
                 percentile = 50,
                 imethod = 1, isd = NULL, parallel = FALSE,
                 smallno = 1.0e-5, zero = "sd")
             
}
\arguments{
    \item{link1}{
    Link function for the first linear predictor.
    By default
    \code{link1 = "identitylink"},
    same as \code{lmean} from
   \code{\link[VGAM:uninormal]{uninormal}}.
   Set \code{link1 = "uninormalQlink"} for normal
   quantile regression.
   See details below.
   
    }
    \item{percentile}{
      Numeric. A vector with the percentiles of interest, between
      0 and 100.
    Used only when \code{link1 = "uninormalQlink"}.
    
    
    }
    \item{lsd, imethod, isd, parallel, smallno, zero}{
    
    Same as in
    \code{\link[VGAM:uninormal]{uninormal}}, except
    that \code{"sd"} is the only accepted value for \code{zero}. 
    
    
    }
    %\item{lmu, lmedian}{
    %Link functions applied to the (positive) mean or median, i.e.,
    %\code{\link[VGAMextra:weibullRMeanlink]{weibullRMeanlink}} or
    %%\code{\link[VGAMextra:weibullRMedianlink]{weibullRMedianlink}}.
  %
   % 
  %  }
  %
  %   \item{imu, imedian}{
  %  Optional initial values for \code{mu} and the \code{median} in the
  % fisrt linear predictor. 
  %}
}
\details{
  An extension of 
  \code{\link[VGAM:uninormal]{uninormal}} adapted to handle
  \code{\link[VGAMextra:uninormalQlink]{uninormalQlink}}, 
  for normal quantile regression (QR) via the first
  linear predictor. 
  
  The standard deviation only can be estimated.
  The second linear predictor is fixed to
  \eqn{\eta_2 = \log \sigma}{\eta_2 = log \sigma}, and
  \code{var.arg} is set internally to \code{FALSE}.
  
   Unlike usual QR where the distribution of \eqn{Y|X} is unspecified,
  \code{uninormalff()} estimates normal distributions at different
  quantiles (as entered in \code{percentile})
  of the \eqn{Y|X}. For this, set
  \code{link1 = uninormaQlink}. To mimic
   \code{\link[VGAM:uninormal]{uninormal}} set
   \code{link1 = "identitylink"} (default).
 
  
  
  Initial developments of this work are in 
  \emph{Miranda & Yee (2019)}.
  See, e.g., 
  \code{\link[VGAMextra:weibullRff]{weibullRff}},
   for another example on distribution specified
   quantile regression with the two--parameter Weibull
   distribution.
  
  %\code{\link[VGAM:gammaR]{gammaR}}.
  %For further choices on link functions for \eqn{\eta_2}{\eta2},
  %see \code{\link[VGAM:Links]{Links}}.  

  %\code{gammaRMeanff} differs from
  %\code{\link[VGAM:gamma2]{gamma2}}. The latter estimates a
  %re-parametrization
  %of the gamma distribution in terms \eqn{\mu}{\mu} and
  %\eqn{\alpha}{shape}. % viz. the density is re--expressed
%  in terms of \eqn{\mu}{\mu} and \eqn{\alpha}{shape}.
  %To he link
  %must be necessarily used via \code{lmu}, which overrides
  %overrides \code{lrate}. 
  %Then, the transformed mean
  %\deqn{{\tt{gammaRMeanlink}}(\alpha; \beta)}{gammaRMeanlink(shape; rate)}
  %are returned as the fitted values, for estimated
  %\eqn{\alpha}{shape} and \eqn{\beta}{rate}.
  
  
}
\value{
  An object of class \code{"vglm"}.
  See \code{\link[VGAM:vglm-class]{vglm-class}} for full details.
  
}
\section{Note}{

  \code{\link[VGAMextra:Q.reg]{Q.reg}} must be used in
  the \code{vglm()} or \code{vgam()} to enter the response.
  See example below.
  
  %The fitted values returned are \code{gammaRMeanlink}--transformed,
  %provided the mean is modelled via \code{lmu}
  
 % The parameters \eqn{\alpha}{shape} and \eqn{\beta}{scale}
%  match the arguments \eqn{shape} and \eqn{scale} from
%  \code{\link[stats:rweibull]{rweibull}}.
  
%  Multiple responses are handled.

  This \bold{VGAM} family function does not handle censored data.
    
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
Miranda & Yee (2019) 
\emph{New Link Functions for Distribution--Specific
Quantile Regression Based on Vector Generalized
Linear and Additive Models}.
Journal of Probability and Statistics, 
Article ID 3493628.

Miranda & Yee (2021) 
\emph{Two--Parameter Link Functions,
With Application to Negative Binomial, Weibull and
Quantile Regression}. In preparation.

%Yee, T. W. (2015)
%\emph{Vector Generalized Linear and Additive Models: 
%With an Implementation in R}. 
%Springer, New York, USA.

}
\seealso{

\code{\link[VGAMextra:uninormalQlink]{uninormalQlink}},
    \code{\link[VGAM:uninormal]{uninormal}},
    \code{\link[VGAMextra:Q.reg]{Q.reg}},
    \code{\link[VGAMextra:weibullQlink]{weibullQlink}},
    \code{\link[VGAMextra:weibullRff]{weibullRff}},
    %\code{\link[base:gamma]{gamma}},
    %\code{\link[VGAMextra:weibullRMedianlink]{weibullRMedianlink}},,
    \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
    %\code{\link[VGAM:weibullR]{weibullR}},
    %\code{\link[VGAM:gamma2]{gamma2}},
    %\code{\link[VGAM:Links]{Links}}.
     
}
\examples{
\dontrun{

x2 <- seq(0,10,length.out = 100)             # independent variable
sig <- exp(0.5 + 0.15*x2)                    # non-constant variance
b_0 <- 10                                    # true intercept
b_1 <- 0.15                                  # true slope
set.seed(17221)                              # make the next line reproducible
e <- rnorm(100,mean = 0, sd = sig)           # normal random error with non-constant variance
y <- b_0 + b_1*x2 + e                        # dependent variable

## Data
ndata <- data.frame(y = y, x2 = x2)

## Some percentiles of interest
percentile <- c(10, 25, 50, 90)

## Normal quantile regression, zero = NULL
fit1 <- vglm(Q.reg(y, length.arg = 4) ~ x2, 
             uninormalff(link1 = "uninormalQlink", percentile = percentile, zero = NULL), 
             data = ndata, trace = TRUE)
#summary(fit1)
( my.coef3Q <- coef(fit1, mat = TRUE) )

## Plots - percentile curves.
plot(y ~ x2, pch = 19, ylim = c(-1, 25), 
main =" Normal quantile regression")
abline(h = -3:25, v = 0, col = "gray", lty = "dashed")
with(ndata, lines(x2, my.coef3Q[1, 1] + my.coef3Q[2, 1] * x2, 
                  col = "red", lty = "dotted", lwd = 4))
with(ndata, lines(x2, my.coef3Q[1, 3] + my.coef3Q[2, 3] * x2, 
                  col = "orange", lty = "dotted", lwd = 4))
with(ndata, lines(x2, my.coef3Q[1, 5] + my.coef3Q[2, 5] * x2, 
                  col = "blue", lty = "dotted", lwd = 4))
with(ndata, lines(x2, my.coef3Q[1, 7] + my.coef3Q[2, 7] * x2, 
                  col = "brown", lty = "dotted", lwd = 4))
legend("topleft", c("90th", "50th", "25th", "10th"),
col = c("brown", "blue", "orange", "red"), lty = rep("dotted", 4), lwd = rep(4, 4))


## Mimicking 'VGAM:uninormal'
fit2 <- vglm(y ~ x2,  uninormalff(link1 = "identitylink", percentile = NULL, zero = NULL), 
             data = ndata, trace = TRUE)

     }
}

