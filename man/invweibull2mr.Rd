\name{invweibull2mr}
\alias{invweibull2mr}
%\alias{weibull3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2- parameter Inverse Weibull Distribution }
\description{

  Maximum likelihood estimation of the 2-parameter Inverse Weibull distribution.
  No observations should be censored.
  
  
}
\usage{

  invweibull2mr(lscale  = loglink, 
                lshape  = logofflink(offset = -2),
                iscale  = NULL, 
                ishape  = NULL, 
                imethod = 2, 
                lss     = TRUE, 
                gscale  = exp(-4:4), 
                gshape  = exp(-4:4),
                probs.y = c(0.25, 0.50, 0.75),
                zero    = "shape")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape}{
  Parameter link functions applied to the 
  (positive) \code{shape} parameter (called \eqn{a}{a} below) and
  (positive) \code{scale} parameter (called \eqn{b}{b} below). 
  Given that the \code{shape} parameter must be greater than \eqn{2},
  \code{lshape = logofflink(offset = -2)} by default.
  See \code{\link[VGAM:Links]{Links}} for more choices.
  
  
  }
  \item{iscale, ishape}{
  Optional initial values for the \code{shape} and \code{scale} parameters.


  }
  \item{gscale, gshape}{
  See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  
  
  }
  \item{lss, probs.y}{
  Details at 
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  
  
  }
  \item{imethod}{
  Initializing method internally implemented. 
  Currently only the values 1 and 2 are allowed and NO observations 
  should be censored.


  }
  \item{zero}{
  Numeric or character vector. The position(s) of the name(s) of 
  the parameters/linear predictors to be
  modeled as intercept--only. Default is \code{"shape"}.
  %If given, their absolutes values must be either 
  %1 or 2. By default, the shape parameter (after \code{lshape} is applied) 
  %is modelled (\code{zero = -2}) as the single unknown number that is estimated. 
  Details: 
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  
  
  }
 % \item{zero, probs.y, lss}{
 %  Details at \code{\link{CommonVGAMffArguments}}.
 % }
}

\details{

  The Weibull distribution and the Inverse Weibull distributions are related
  as follows:
  
  Let \eqn{X}{X} be a Weibull random variable with paramaters scale =\eqn{b}{b}
  and shape =\eqn{a}{a}. Then, the random variable \eqn{Y = 1/X}{Y = 1/X} has 
  the Inverse Weibull density with parameters scale = \eqn{1/b}{1/b} and 
  shape = \eqn{a}{a}. 

  The  Inverse weibull density for a response \eqn{Y} is given by
  \deqn{f(y;a,b) = a (b^a) y^{-a-1} \exp[-(y/b)^(-a)] }{%
        f(y;a,b) = a (b^a) y^(-a-1) * exp(-(y/b)^(-a))}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{y > 0}.
%  The cumulative distribution function is 
%  \deqn{F(y;a,b) = 1 - \exp[-(y/b)^a].}{%
%        F(y;a,b) = 1 - exp(-(y/b)^a).}
%  The mean of \eqn{Y} is \eqn{b \, \Gamma(1 - 1/a)}{b * gamma(1 - 1/a)}
 The mean, that is returned as the fitted values, (if \eqn{a > 1}{a > 1}) and 
 the variance (if \eqn{a > 2}{a >2}) are
  \deqn{E[Y] = b \ \Gamma(1 - 1/a); \ \ \ Var[Y] = b^{2} \ 
  [\Gamma(1 - 2/a) - (\Gamma(1 - 1/a))^2].}{%
        E[Y] = b * \Gamma(1 - 1/a),  Var[Y] = (b^2) * 
        [\Gamma(1 - 2/a) - (\Gamma(1 - 1/a))^2].}
  
  Fisher scoring is used to estimate both parameters.
  Although the expected information matrices used are valid in all regions of 
  the parameter space, the regularity conditions for maximum likelihood 
  estimation (MLE) are satisfied only if \eqn{a>2}{a>2} 
  (according to Kleiber and Kotz (2003)). If this is violated then a warning 
  message is issued. To enforce \eqn{a > 2}{a > 2}, it has been set by default 
  that \code{lshape = logofflink(offset = 2)}.
  
  As a result of the math relation between the Weibull and the Inverse Weibull
  distributions, regularity conditions for inference for the latter, are the
  following:
  %Harper et a. (2011) summarizes the MLE's conditions for inference: 
  if \eqn{a \le 1} then the MLE's are not consisten, 
  if \eqn{1 < a < 2}{1 < a < 2} then MLEs exist but are not assymptotically 
  normal, if \eqn{a = 2}{a = 2}, the MLE's exist and are normal and 
  asymptotically efficient but the convergence rate is slower compared when 
  \eqn{a > 2}{a > 2}. If \eqn{a > 2}, then the MLE's have classical
  asymptotic properties.
  %choosing \code{lshape = logoff(offset = -2)}.
  %Common values of the shape parameter lie between 0.5 and 3.5.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOT SURE... CHECK THIS UP. IT'S VERY IMPORTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Summarized in Harper et al. (2011),
  %for inference, there are 4 cases to consider.
  %If \eqn{a \le 1} then the MLEs are not consistent
  %(and the smallest observation becomes a hyperefficient
  %solution for the location parameter in the 3-parameter case).
  %If \eqn{1 < a < 2} then MLEs exist but are not asymptotically normal.
  %If \eqn{a = 2} then the MLEs exist and are normal and asymptotically
  %efficient but with a slower convergence rate than when \eqn{a > 2}.
  %If \eqn{a > 2} then MLEs have classical asymptotic properties.
  %The 3-parameter (location is the third parameter) Weibull can
  %be estimated by maximizing a profile log-likelihood (see,
  %e.g., Harper et al. (2011) and Lawless (2003)), else try
  %\code{\link{gev}} which is a better parameterization.

}

\value{

  An object of class \code{"vglmff"} 
  (see \code{\link[VGAM:vglmff-class]{vglmff-class}}).
  The object is used to model special models such as 
  \code{\link[VGAM:vglm]{vglm}} and \code{\link[VGAM:vgam]{vgam}}.


}

\references{

Harper, W. V., Eschenbach, T. G. and James, T. R. (2011)
Concerns about Maximum Likelihood Estimation for the Three-Parameter
Weibull Distribution: Case Study of Statistical Software.
\emph{The American Statistician}, \bold{65(1)}, 44-54.

Kleiber, C. and Kotz, S. (2003)
\emph{Statistical Size Distributions in Economics and Actuarial Sciences},
Hoboken, NJ, USA: Wiley-Interscience.


Johnson, N. L. and Kotz, S. and Balakrishnan, N. (1994)
\emph{Continuous Univariate Distributions},
2nd edition, Volume 1, New York: Wiley.

%Lawless, J. F. (2003)
%\emph{Statistical Models and Methods for Lifetime Data},
%2nd ed.
%{Hoboken, NJ, USA: John Wiley & Sons}.
%Rinne, Horst. (2009)
%\emph{The Weibull Distribution: A Handbook}.
%Boca Raton, FL, USA: CRC Press.
%Gupta, R. D. and Kundu, D. (2006)
%On the comparison of Fisher information of the
%Weibull and GE distributions,
%\emph{Journal of Statistical Planning and Inference},
%\bold{136},
%3130--3144.
%Harper, W. V. and Eschenbach, T. G. and James, T. R. (2011)
%Concerns about Maximum Likelihood Estimation for
%              the Three-Parameter {W}eibull Distribution:
%              Case Study of Statistical Software,
%\emph{The American Statistician},
%\bold{65(1)},
%{44--54}.
% Smith, R. L. (1985)
%  Maximum likelihood estimation in a class of nonregular cases.
%  \emph{Biometrika}, \bold{72}, 67--90.
% Smith, R. L. and Naylor, J. C. (1987)
%  A comparison of maximum likelihood and Bayesian estimators
%  for the three-parameter Weibull distribution.
%  \emph{Applied Statistics}, \bold{36}, 358--369.


}

\author{ 

Victor Miranda and T. W. Yee.

}

\note{

  By default, the shape paramater is modeled as intercept only.
  
  Successful convergence depends on having reasonably good initial
  values. If the initial values chosen by this function are not good,
  make use the two initial value arguments, \code{iscale} and \code{ishape}.

  This \pkg{VGAM} family function currently handles \emph{multiple responses} 
  however, it does not handle censored data. This feature will be considered in 
  a later version of the package.
  
  The Inverse Weibull distribution, which is that of \eqn{Y = 1/X}{Y = 1/X} 
  where \eqn{X}{X} has the Weibull density, is known as the 
  log-Gompertz distribution.
  %This \pkg{VGAM} family function handles multiple responses.
  %The Weibull distribution is often an alternative to the lognormal
  %distribution. The inverse Weibull distribution, which is that of
  %\eqn{1/Y} where \eqn{Y} has a Weibull(\eqn{a,b}) distribution, is
  %known as the log-Gompertz distribution.
%  There are problems implementing the three-parameter Weibull
%  distribution. These are because
%  the classical regularity conditions for the
%  asymptotic properties of the MLEs are not satisfied because the
%  support of the distribution depends on one of the parameters.

% Other related distributions are the Maxwell and Rayleigh
% distributions.


}

\section{Warning}{

  Note that \code{zero} can be a \bold{numeric} or a \bold{character} 
  vector specifying the position of the names 
  (partially or not) of the linear predictor modeled as intercept only.
  In this family function these names are
  
  \code{c("scale", "shape")}.
  
  Numeric values can be entered as usual.
  See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}} 
  for further details.
  For simplicity, the second choice is recommended.

  %This function is under development to handle censoring situations.
  %The version of this function which will handle censored data will be
  %called \code{cenweibull()}. It is currently being written and will use
  %\code{\link{SurvS4}} as input. 
  %It should be released in later versions of \pkg{VGAMextra}.


  If the shape parameter is less than two (i.e. less than 
  \code{exp(0.69315)}), then misleading inference may
  result ! (see above the regularity condition for the 'variance'), 
  e.g., in the \code{summary} and \code{vcov} of the object.
  
  However, the larger the shape parameter is (for instance, greater than 
  \code{exp(2.5)}, plus reasonable scale), the more unstable the algorithm
  may become. The reason is that inverse weibull densities under such 
  conditions are highly peaked and left skewed. Thus, density values are 
  too close to zero (or values represented as zero in computer arithmetic).


}

\seealso{

    \code{\link[VGAMextra:invweibullDist]{invweibullDist}},
    \code{\link[VGAM:weibullR]{weibullR}}.
    %\code{\link[stats:Weibull]{dweibull}},
    %\code{\link{truncweibull}},
    %\code{\link{gev}},
    %\code{\link{lognormal}},
    %\code{\link{expexpff}},
    %\code{\link{maxwell}},
    %\code{\link{rayleigh}},
    %\code{\link{gumbelII}}.


}

\examples{
#-----------------------------------------------------------------------#
# Here, covariate 'x2' affects the scale parameter.
# See how data is generated.

set.seed(102)
wdata <- data.frame(x2 = runif(nn <- 1000))  # Complete data
wdata <- transform(wdata,
            y1 = rinvweibull(nn, scale = exp(2.5 - (0.5) * x2), 
                             shape = exp(1.5) ),
                             
            y2 = rinvweibull(nn, scale = exp(1.5 + x2), 
                             shape = exp(1.25) ))
            
#------------------------------------------------------------------------#
# Fitting the Inverse Weibull distribution accordingly.
# Note that multiple responses are handled. 

fit1 <- vglm(cbind(y1, y2) ~ x2, 
             invweibull2mr(zero = "shape",
                           # OPTIONAL INITIAL VALUE. Be carefull here when
                           # entered initial value. Sensitive distribution
                             ishape = exp(1.2),
                          lss = TRUE),                           
             data = wdata, trace = TRUE, crit = "log")

coef(fit1, matrix = TRUE)
vcov(fit1)
summary(fit1)

###   A second option (producing same results!!) might be to use the 
###   constraints argument in the 'vglm()' call. Note that 'x2' affects
###   the scale parameter only.

fit2 <- vglm(y1 ~ x2, 
             invweibull2mr(zero = NULL), 
             data = wdata, trace = TRUE,
             constraints = list(x2 = rbind(1, 0)))
            
coef(fit2, matrix = TRUE)
vcov(fit2)
summary(fit2)
constraints(fit2)



}
\keyword{Inverse Weibull distribution}
\keyword{Fisher scoring}
\keyword{VGAM}
