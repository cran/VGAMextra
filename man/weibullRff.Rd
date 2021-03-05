\name{weibullRff}
\alias{weibullRff}
%\alias{weibullRMeanff}
%\alias{weibullRMedianff}
%\alias{weibullRMean}
%\alias{weibullRMedian}

\title{
 
     Distribution--specified quantile regression: 
     2--parameter Weibull Distribution
    
}
\description{
   %
     
   Estimates the 2--parameter Weibull distribution by maximum likelihood.
   An extension of \code{weibullR} from \pkg{VGAM}.
   Weibull quantile regression and Weibull--mean
   modelling are also handled via the first
   linear predictor.
   %
}
\usage{
    %
     weibullRff(link1 = c("loglink", "weibullMlink", "weibullQlink")[1],
                lshape = "loglink", percentile = 50,
                imu = NULL, iscale = NULL, ishape = NULL,
                lss = TRUE, nrfs = 1, probs.y = c(0.2, 0.5, 0.8),
                imethod = 1, zero = "shape")
             
}
\arguments{
    \item{link1}{
    Link function for the first linear predictor. 
    Default
    is \code{link1 = "loglink"}, mimicking 
    \code{\link[VGAM:weibullR]{weibullR}}. The other options
    are the 2--parameter
    \code{\link[VGAMextra:weibullQlink]{weibullQlink}}, applied
    to the Weibull quantile function, and the 2--parameter
    \code{\link[VGAMextra:weibullMlink]{weibullMlink}}, applied to 
    the Weibull mean function. See below for more details.
    
    }
    \item{percentile}{
    
    Numeric. A vector with the percentiles of interest, between 0
    and 100.
    Used only in Weibull quantile regression, that is,
    when \code{link1 = "weibullQlink"}.
    
    }
    \item{lshape,imu, iscale, ishape, lss,
    nrfs, probs.y, imethod}{
    
    Same as \code{\link[VGAM:weibullR]{weibullR}}.
    }
    \item{zero}{
     Specifies the parameters to be modelled as intercept--only.
     Further details below.
     
     See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
    
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
  
  \code{weibullRff} is a modified version of
  \code{\link[VGAM:weibullR]{weibullR}} adapted to handle
  \code{\link[VGAMextra:weibullQlink]{weibullQlink}}
  and \code{\link[VGAMextra:weibullMlink]{weibullMlink}},
  two 2-parameter linear predictors that model
  the Weibull mean and quantiles respectively.
  %The densit\code{\link[VGAM:weibullR]{weibullR}}
  %with \eqn{\beta}{scale} 
  The underlying density is the ordinary 
  \code{scale}(\eqn{\beta}) & \code{shape}(\eqn{\alpha})
  Weibull density (see \code{\link[VGAM:weibullR]{weibullR}}).
  
  The second linear predictor is always 
  \eqn{\eta_2 = \log~\alpha}{\eta2 = log shape}.
  The argument \code{link1} handles the first linear predictor.
  
  \pkg{** Mimicking \code{weibullR} **}
  
  The default is \code{link1 = "loglink"}, i.e., 
  \eqn{\eta_1 = \log~\beta = \log~scale}{\eta1 = log scale}, and 
  \eqn{\eta_2 = \log~\alpha = \log~shape}{\eta2 = log shape},
  as with \code{\link[VGAM:weibullR]{weibullR}}.
  The mean (\eqn{\mu}) is returned as the fitted value.
  
  \pkg{** Weibull quantile regression **}
  
  For Weibull quantile regression set \code{link1 = "weibullQlink"}
  and enter a numeric vector of percentiles of interest
  via \code{percentile}.  
  See examples.
  
  NOTE: Enter the response using 
  \code{\link[VGAMextra:Q.reg]{Q.reg}}. See example below.
    The Weibull quantiles are returned as the fitted values.
  
  \pkg{** Weibull-mean modelling **}
  
  For Weibull-mean modelling (viz. mean time to failure)
  set \code{link1 = "weibullMlink"}. 
  The mean (\eqn{\mu}) is returned as the fitted value.


  %Particularly, \code{lss} works exactly as in
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
  
  %The fitted values returned are \code{gammaRMeanlink}--transformed,
  %provided the mean is modelled via \code{lmu}
  
  The parameters \eqn{\alpha}{shape} and \eqn{\beta}{scale}
  match the arguments \eqn{shape} and \eqn{scale} from
  \code{\link[stats:rweibull]{rweibull}}.
  
  Multiple responses are handled.

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
    \code{\link[VGAMextra:Q.reg]{Q.reg}},
    \code{\link[VGAMextra:weibullQlink]{weibullQlink}},
    \code{\link[VGAMextra:weibullMlink]{weibullMlink}},
    \code{\link[VGAM:weibullR]{weibullR}},
    \code{\link[base:gamma]{gamma}},
    %\code{\link[VGAMextra:weibullRMedianlink]{weibullRMedianlink}},,
    \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
    %\code{\link[VGAM:weibullR]{weibullR}},
    %\code{\link[VGAM:gamma2]{gamma2}},
    %\code{\link[VGAM:Links]{Links}}.
     
}
\examples{
\dontrun{
set.seed(18121)
nn <- 300
x2 <- sort(runif(nn, 0, 3))  # Predictor/covariate.
bb <- exp(1.1 + 0.2 * x2)    # Scale parameter as function of x2.
aa <- exp(1.0 - 0.35 * x2)     # Shape parameter as function of x2.
mymu <- bb * gamma(1 + 1/aa)  # The Weibull mean.

## Use weibullMlink to generate appropriate scale parameter.
newbb <- weibullMlink(theta = log(mymu), shape = aa, inverse = TRUE, deriv = 0)

## A single random response
wdata <- data.frame(y1 = rweibull(nn, shape = aa, scale = newbb), x2 = x2)

# Plotting the data / Histogram
plot(y1  ~ x2, xlim = c(0, 3.1), ylim = c(-1, 35),
     pch = 20, data = wdata, col = "black", 
     main = "Weibull Quantile regression~ x2")
abline(h = 0, v = 0, col = "grey", lty = "dashed")
with(wdata, hist(y1, col = "red", breaks = 15))

## Weibull regression - percentile = c(25, 50, 75)
## Note the use of Q.reg.
fit1 <- vglm(Q.reg(y1, length.arg = 3) ~ x2, 
             weibullRff(link1 = "weibullQlink", zero = NULL,
                                 percentile = c(25, 50, 75)), 
             trace = TRUE, data = wdata)
head(fitted(fit1))
summary(fit1)
my.coef3Q <- coef(fit1, mat = TRUE)

### Proportion of data below the estimated 25\% Quantile line.
100 * (1 - (sum(wdat$y1 >= fitted(fit2)[, 1]) / nn))  # Around 25\%
### Proportion of data below the estimated 50\% Quantile line.
100 * (1 - (sum(wdat$y1 >= fitted(fit2)[, 2]) / nn))   # Around 50\%
### Proportion of data below the estimated 75\% Quantile line.
100 * (1 - ( sum(wdat$y1 >= fitted(fit2)[, 3]) / nn ))   # Around 75\%

## The quantile plots ##
my.coef3Q <- coef(fit2, matrix = TRUE)
with(wdat, lines(x2, exp(my.coef3Q[1, 1] + my.coef3Q[2, 1] * x2), 
                    col = "red", lty = "dotted", lwd = 4))
with(wdat, lines(x2, exp(my.coef3Q[1, 3] + my.coef3Q[2, 3] * x2), 
                 col = "orange", lty = "dotted", lwd = 4))
with(wdat, lines(x2, exp(my.coef3Q[1, 5] + my.coef3Q[2, 5] * x2), 
                 col = "blue", lty = "dotted", lwd = 4))

## Adding the 'mean' or expected Weibull regression line.
fit2 <- vglm(y1 ~ x2, 
             weibullRff(link1 = "weibullMlink", zero = NULL), 
             trace = TRUE, data= wdat)
my.coef3Q <- coef(fit2, mat = TRUE)
with(wdat, lines(x2, exp(my.coef3Q[1, 1] + my.coef3Q[2, 1] * x2), 
                 col = "yellow", lty = "dashed", lwd = 3))


legend("topleft", c("25h Perc", "50th Perc", "Mean", "75th Perc"),
       col = c("red", "orange", "cyan", "blue"),
       lty = c("dashed", "dashed", "dashed", "dashed"), lwd = rep(4, 4))
     }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\dontrun{
%  The defaul linear predictors are:
%
%  a) \eqn{\eta_1 = {\tt{weibullRMeanlink}}(\alpha; \beta) = \log \mu =
%  \log (\beta * \Gamma(1 + 1 / \alpha))}{
%  \eta1 = weibullRMeanlink(scale; shape) = log \mu =
%  log (scale * \Gamma(1 + 1/\alpha))}, and
%  \eqn{\eta_2 = \log \alpha}{\eta2 = log (shape)}, for
%  \code{weibullRMeanff}, and
%  
%  b) \eqn{\eta_1 = {\tt{weibullRMedianlink}}(\alpha; \beta) = \log (Median) =
%  \log (\beta * (\log 2)^(1/\alpha))}{
%  \eta1 = weibullRMedianlink(scale; shape) = log (Median) =
%  log (scale * (log 2)^(1 / \alpha)) }, and
%  \eqn{\eta_2 = \log \alpha}{\eta2 = log (shape)}, for
%  \code{weibullRMedianff},
%  
%  unlike \eqn{\eta_1 = \log \beta}{\eta1 = log (scale)} and
%  \eqn{\eta_2 = \log \alpha}{\eta2 = log (shape)}
%  from \code{\link[VGAM:weibullR]{weibullR}}.
%  Here, \eqn{\beta}{'scale'} and 
%  \eqn{\alpha}{'shape'} are scale and shape parameters.
%  
%  \code{lmu}, \code{lmedian} override \code{lscale} and can only be assigned 
%  to either \code{\link[VGAMextra:weibullRMeanlink]{weibullRMeanlink}} and
%  \code{\link[VGAMextra:weibullRMedianlink]{weibullRMedianlink}}.
  %\code{\link[VGAMextra:gammaRMeanlink]{gammaRMeanlink}} is a valid entry
  %({\code{lmu}}). 
%  To mimic \code{\link[VGAM:weibullR]{weibullR}} 
%  set \code{lmu = NULL} and \code{lscale = "loglink"}.
%  The mean (\eqn{\mu}) and the median are returned as the fitted values
%  by \code{weibullRMeanff} and \code{weibullRMedianff} 
%  respectively.
%  }
