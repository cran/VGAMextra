\name{invgamma2mr}
\alias{invgamma2mr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2 - parameter Inverse Gamma Distribution }
\description{

  Estimates the 2-parameter Inverse Gamma distribution
  by maximum likelihood estimation.

}
\usage{
  invgamma2mr(lmu      = "loge", 
              lshape   = logoff(offset = -2), 
              parallel = FALSE, 
              ishape   = NULL, 
              imethod  = 1, 
              zero     = "shape") 

}
%- maybe also 'usage' for other objects documented here.
% apply.parint = FALSE,
\arguments{
  \item{lmu, lshape}{
  Link functions applied to the (positives) \emph{mu} and \emph{shape}
  parameters (called \eqn{\mu}{mu} and \eqn{a}{shape} respectively), 
  according to \code{\link[VGAM:gamma2]{gamma2}}.
  See \code{\link{CommonVGAMffArguments}} for further information.


  }
  \item{parallel}{
  Same as \code{\link[VGAM:gamma2]{gamma2}}.
  Details at \code{\link{CommonVGAMffArguments}}.  
  %If \code{parallel = TRUE} then the constraint is not applied 
  %to the intercept.


  }
  \item{ishape}{
  Optional initial value for \emph{shape}, same as 
  \code{\link[VGAM:gamma2]{gamma2}}
  %A \code{NULL} means a value is computed internally.
  %If a failure to converge occurs, try using this argument.
  %This argument is ignored if used within \code{\link{cqo}}; see the
  %\code{iShape} argument of \code{\link{qrrvglm.control}} instead.


  }
  \item{imethod}{
  Same as \code{\link[VGAM:gamma2]{gamma2}}.
  %An integer with value \code{1} or \code{2} which
  %specifies the initialization method for the \eqn{\mu}{mu} parameter.
  %If failure to converge occurs try another value (and/or specify a 
  %value for \code{ishape}).


  }
  \item{zero}{
  Numeric or character vector. Position or name(s) of the 
  parameters/linear predictors to be 
  modeled as intercept--only. Default is \code{"shape"}. Details at 
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  %Specifies which of the two linear/additive predictors 
  %are modelled as an intercept only. By default, the shape parameter 
  %(after \code{lshape} is applied) is modelled as a single unknown 
  %number that is estimated.  It can be modelled as a function of the 
  %explanatory variables by setting \code{zero = NULL}.  A negative
  %value means that the value is recycled, so setting \eqn{-2} means 
  %all shape parameters are intercept only.
  %See \code{\link{CommonVGAMffArguments}} for more information.


  }
}
\details{

 The Gamma distribution and the Inverse Gamma distribution are related 
 as follows:Let \emph{X} be a random variable distributed as 
 \eqn{Gamma (a, \beta)}{Gamma (a , b)},  where \eqn{a > 0}{a > 0} 
 denotes the \emph{shape} parameter and \eqn{\beta > 0}{b > 0} is the 
 \emph{scale} paramater. Then \eqn{Y = 1/X}{Y = 1/X} is an Inverse Gamma 
 random variable with parameters \emph{scale} = \eqn{a}{a} and 
 \emph{shape} = \eqn{1/\beta}{1/b}. 
 
 The Inverse Gamma density function is given by
 
  \deqn{f(y;\mu, a) = \frac{(a - 1)^{a} \mu^{a}}{\Gamma(a)}y^{-a- 1} \
                    e^{-\mu(a - 1)/y},}{%
        f(y; mu, shape) = exp(-mu * (shape - 1)/y) * y^(-shape - 1) * 
                    mu^(shape) * (shape - 1)^(shape) / gamma(shape),}
  
%  \deqn{f(y;\mu,a) = \frac{\exp(-a y / \mu) \times
%               (a y / \mu)^{a-1}
%               \times a}{
%               \mu \times \Gamma(a)}}{%
%   f(y;mu,shape) = exp(-shape * y / mu) y^(shape-1) shape^(shape) / 
%          [mu^(shape) * gamma(shape)]}

  for \eqn{\mu > 0}{mu > 0}, \eqn{a > 0}{shape > 0} and \eqn{y > 0}.
  Here, \eqn{\Gamma(\cdot)}{gamma()} is the gamma function, as in 
  \code{\link[base:Special]{gamma}}. The mean of \emph{Y} is 
  \eqn{\mu=\mu}{mu=mu} (returned as the fitted values) with variance 
  \eqn{\sigma^2 = \mu^2 / (a - 2)}{sigma^2 = mu^2 / (shape - 2)} 
  if \eqn{a > 2}{shape > 2}, else is infinite. Thus, the 
  \emph{link function} for the \emph{shape} parameter is 
  \code{\link[VGAM:loglog]{loglog}}. Then, by default, the two 
  linear/additive predictors are \eqn{\eta_1=\log(\mu)}{eta1=log(mu)}, 
  and \eqn{\eta_2=\log(a)}{eta2=loglog(shape)}, i.e in the VGLM context, 
  \eqn{\eta = (log(\mu), loglog(a)}{eta = (log(mu), loglog(shape)).}
  
  This \pkg{VGAM} family function handles \emph{multiple} reponses by 
  implementing Fisher scoring and unlike 
  \code{\link[VGAM:gamma2]{gamma2}}, the working-weight matrices
  are \emph{not} diagonal.
  % Not 100% sure about weights.
  % On 20141208 I'm 100% sure about weights.  
  The Inverse Gamma distribution is right-skewed and either for small values
  of \eqn{a}{shape} (plus modest \eqn{\mu}{mu}) or very large values of 
  \eqn{\mu}{mu} (plus moderate \eqn{a > 2}{shape > 2}), the density has 
  values too close to zero.
 
}
\value{

  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}}
  and \code{\link{vgam}}.

}
\references{

  McCullagh, P. and Nelder, J. A. (1989)
  \emph{Generalized Linear Models}, 
  2nd ed. London, UK. Chapman & Hall.
  %description in \code{\link[VGAM:gamma2]{gamma2}}
  %McCullagh, P. and Nelder, J. A. (1989)
  %\emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.

}
\section{Warning}{


  Note that \code{zero} can be a \bold{numeric} or a \bold{character} 
  vector specifying the position of the names 
  (partially or not) of the linear predictor modeled as intercept only.
  In this family function such names are
  
  \code{c("mu", "shape")}.
  
  Numeric values can be entered as usual.
  See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}} 
  for further details.




}
\author{ Victor Miranda and T. W. Yee }

\note{

  The response must be strictly positive.
  %A moment estimator for the shape parameter may be implemented in
  %the future.
  
  If \code{mu} and \code{shape} are vectors, then \code{rinvgamma(n = n,
  shape = shape, scale = mu/(shape - 1)} will generate random inverse gamma 
  variates of this parameterization, etc.; 
  see \code{\link[VGAMextra:invgammaDist]{invgammaDist}}.
  %For \code{\link{cqo}} and \code{\link{cao}}, taking the logarithm
  %of the response means (approximately) a \code{\link{gaussianff}} family
  %may be used on the transformed data.

  Given the math relation between the Gamma and the Inverse Gamma
  distributions, the parameterization of this \pkg{VGAM} family function 
  underlies on the parametrization of the 2-parameter gamma distribution 
  described in the monograph

}
\seealso{

  \code{\link[VGAMextra:invgammaDist]{invgammaDist}},
  \code{\link[VGAM:gamma2]{gamma2}} for the 2-parameter gamma distribution,
  \code{\link[stats:GammaDist]{GammaDist}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  %\code{\link[VGAM:loge]{loge}}.
  
}
\examples{
#------------------------------------------------------------------------#
# Essentially fitting a 2-parameter inverse gamma distribution
# with 2 responses.

set.seed(101)
y1 = rinvgamma(n = 500, scale = exp(2.0), shape = exp(2.0))
y2 = rinvgamma(n = 500, scale = exp(2.5), shape = exp(2.5))
gdata <- data.frame(y1, y2)

fit1 <- vglm(cbind(y1, y2) ~ 1, 
            family = invgamma2mr(zero = NULL, 
            
                                 # OPTIONAL INITIAL VALUE
                                 # ishape = exp(2),
                                 
                                 imethod = 1),
            data = gdata, trace = TRUE)

Coef(fit1)
c(Coef(fit1), log(mean(gdata$y1)), log(mean(gdata$y2)))
summary(fit1)
vcov(fit1, untransform = TRUE)

#------------------------------------------------------------------------#
# An example including one covariate.
# Note that the x2 affects the shape parameter, which implies that both,
# 'mu' and 'shape' are affected.
# Consequently, zero must be set as NULL !

x2    <- runif(1000)
gdata <- data.frame(y3 = rinvgamma(n = 1000, 
                                   scale = exp(2.0), 
                                   shape = exp(2.0 + x2)))

fit2 <- vglm(y3 ~ x2, 
            family = invgamma2mr(lshape = "loge", zero = NULL), 
            data = gdata, trace = TRUE)

coef(fit2, matrix = TRUE)
summary(fit2)
vcov(fit2)

}
\keyword{VGAM}
\keyword{Inverse Gamma distribution}
\keyword{Fisher scoring}

