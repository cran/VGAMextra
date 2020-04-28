\name{gen.betaIImr}
\alias{gen.betaIImr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generalized Beta Distribution of the Second Kind family function }
\description{

  Maximum likelihood estimation of the 4-parameter 
  generalized beta II distribution using Fisher scoring.

}
\usage{
 gen.betaIImr(lscale    = "loglink", 
              lshape1.a = "loglink", 
              lshape2.p = "loglink", 
              lshape3.q = "loglink", 
              iscale    = NULL, 
              ishape1.a = NULL,
              ishape2.p = NULL, 
              ishape3.q = NULL, 
              imethod   = 1,
              lss       = TRUE, 
              gscale    = exp(-5:5), 
              gshape1.a = exp(-5:5),
              gshape2.p = exp(-5:5), 
              gshape3.q = exp(-5:5), 
              probs.y   = c(0.25, 0.50, 0.75),
              zero      = "shape" )            
 }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lscale, lshape1.a, lshape2.p, lshape3.q}{
  Parameter link functions applied to the 
  shape parameter \code{a},
  scale parameter \code{scale},
  shape parameter \code{p}, and
  shape parameter \code{q}.
  All four parameters are positive.
  See \code{\link{Links}} for more choices.
  
  
  }
  \item{iscale, ishape1.a, ishape2.p, ishape3.q}{
  Optional initial values for \code{b}, \code{a}, \code{p} and \code{q}.
  Default is \code{NULL} for all of them, meaning initial values are 
  computed internally.
  
  
  }
  \item{imethod}{
  Initializing method to internally compute the initial values.
  Currently, only \code{method = 1} is handled.
  
  
  }
  \item{gscale, gshape1.a, gshape2.p, gshape3.q}{
  Grid search initial values.
  See \code{\link{CommonVGAMffArguments}} for further information.
  
  
  }
  \item{zero}{
  Numeric or Character vector. Position(s) or name(s) of the 
  parameters/linear predictors to be modeled as intercept-only.
  Details at 
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  %The values must be from the set \{1, 2, 3, 4\}.
  %Here, the shape1.a parameter is modelled as intercept-only by default.
  
  
  }
  \item{lss, probs.y}{
  See \code{\link{CommonVGAMffArguments}} for important information.
  
  
  } 
}
\details{

This distribution is most useful for unifying a substantial
number of size distributions. For example, the Singh-Maddala,
Dagum, Fisk (log-logistic), Lomax (Pareto type II),
inverse Lomax, beta distribution of the second kind
distributions are all special cases.
Full details can be found in Kleiber and Kotz (2003), and
Brazauskas (2002). The argument names given here are used 
by other families that are special cases of this family.
Fisher scoring is used here and for the special cases too.

The 4-parameter generalized beta II distribution has density
  \deqn{f(y) = a y^{ap-1} / [b^{ap} B(p,q) \{1 + (y/b)^a\}^{p+q}]}{%
        f(y) = a y^(ap-1) / [b^(ap) B(p,q) (1 + (y/b)^a)^(p+q)]}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{p > 0}, \eqn{q > 0}, 
                           \eqn{y \geq 0}{y >= 0}.
Here \eqn{B} is the beta function, and 
\eqn{b} is the scale parameter \code{scale},
while the others are shape parameters.
The mean is 
  \deqn{E(Y) = b \, \Gamma(p + 1/a) \, \Gamma(q - 1/a) / 
                                (\Gamma(p) \, \Gamma(q))}{%
        E(Y) = b  gamma(p + 1/a)  gamma(q - 1/a) / ( gamma(p)  gamma(q))}
provided \eqn{-ap < 1 < aq}; these are returned as the fitted values.
%The distribution is motivated by the incomplete beta function
%\eqn{B_y(p,q)} which is the integral from 0 to \eqn{y} of the integrand
%\eqn{u^{p-1} (1-u)^{q-1}}{u^(p-1) (1-u)^(q-1)} where \eqn{y>0}.


}
\section{Warning}{

\code{zero} can be a \bold{numeric} or a \bold{character} vector 
specifying the position(s) or the name(s) (partially or not) of the linear 
predictors modeled as intercept--only. Numeric values can be 
entered as usual. If names are used, note that the linear predictors in 
this family function are

\code{c("scale", "shape1.a", "shape2.p", "shape3.q")}.

For simplicity, using names rather than numeric vectors is 
recommended.

}
\value{

  An object of class \code{"vglmff"} (see \code{\link{vglmff-class}}).
  The object is used by modelling functions such as \code{\link{vglm}},
  and \code{\link{vgam}}.


}
\references{

 Brazauskas, V. (2002)
 Fisher information matrix for the Feller-Pareto distribution.
 \emph{Statistics & Probability Letters}, 
 \bold{59}, 159--167.
 
 Kleiber, C. and Kotz, S. (2003) 
 \emph{Statistical Size Distributions in Economics and Actuarial Sciences}.
 Wiley Series in Probability and Statistics. 
 Hoboken, New Jersey, USA.

 McDonald, J. B. and Xu, Y. J. (1995)
 A generalization of the beta distribution with applications.
 \emph{Journal of Econometrics}. 
 \bold{66}, p.133--152.
  
 McDonald, J. B. (1984)
 Some generalized functions for the size distribution of income.
 \emph{Econometrica}, 
 \bold{52}, 647--663.

}

\author{ T. W. Yee and V. Miranda.}

\note{

Paramaters \code{"shape1.a", "shape2.p", "shape3.q"} are modeled 
as intercept only, by default.

If the self-starting initial values fail, try experimenting
with the initial value arguments, \code{iscale}, \code{ishape1.a},
\code{ishape2.p} and \code{ishape3.q} whose default is \code{NULL}.
Also, the constraint \eqn{-ap < 1 < aq} may be violated as the 
iterations progress so it is worth monitoring convergence, e.g., 
set \code{trace = TRUE}.

Successful convergence depends on choosing good initial values. This 
process might be difficult for this distribution, since 4 parameters are 
involved. Presently, only \code{method = 1} is internally handled to set 
initial values. It involves \emph{grid search}, an internal implementation 
of the well-known grid search algorithm for exhaustive searching through 
a manually specified subset of the hyperparameter space.
%Aiming to provide more options to the user, two 
%sefl-initializing methods have been set, i.e. 
%\code{imethod = 1 or 2}. In particular, \code{imethod = 2} involves
%the Gauss's Hypergeometric Function referred to as 
%\eqn{{}_{2}F_{1}}{2F1}, whose algorithm is included in the package
%\code{hypergeo}. Note that packages \code{contfrac} and 
%\code{elliptic} must be additionally installed for suitable operation
%of \code{hypergeo}.

Default value of \code{lss} is \code{TRUE} standing for the 
following order: location (b), shape1.a (a), shape2.p (p), 
shape3.q (q). In order to match the arguments of existing R functions, 
the option \code{lss = FALSE} might be set leading to switch the 
position of location (b) and shape1.a (a), only.


}
\seealso{

    \code{\link{betaff}},
    \code{\link{betaII}},
    \code{\link{dagum}},
    \code{\link{sinmad}},
    \code{\link{fisk}},
    \code{\link{lomax}},
    \code{\link{inv.lomax}},
    \code{\link{paralogistic}},
    \code{\link{inv.paralogistic}},
    %\code{\link{lino}},
    \code{\link{genbetaIIDist}}.
    
}
\examples{

%\dontrun{
#----------------------------------------------------------------------- #
# An example.- In this data set parameters 'shape1.a' and 'shape3.q' are 
# generated in terms of x2.

set.seed(1003)
nn <- 200
gdata1  <- data.frame(x2 = runif(nn))
gdata   <- transform(gdata1,
           y1 = rgen.betaII(nn, scale = exp(1.1), shape1.a = exp(1.2 + x2), 
                            shape2.p = exp(0.7) , shape3.q = exp(2.1 - x2)),
           y2 = rgen.betaII(nn, scale = exp(2.0), shape1.a = exp(1.8 + x2),
                            shape2.p = exp(2.3) , shape3.q = exp(1.9 - x2)), 
           y3 = rgen.betaII(nn, scale = exp(1.5), shape1.a = exp(1.8),
                            shape2.p = exp(2.3) , shape3.q = exp(1.3)))
                            
#------------------------------------------------------------------------#
# A single intercept-only model. No covariates.
# Note the use of (optional) initial values.
fit  <- vglm(y2 ~ 1,   #y3 ~ 1
             gen.betaIImr(lss = TRUE,
                          # OPTIONAL INITIAL VALUES
                           #iscale    = exp(1.5), 
                           #ishape1.a = exp(1.8),
                           #ishape2.p = exp(2.3),
                           #ishape3.q = exp(1.3),
                          
                          imethod = 1),  
             data = gdata, trace = TRUE, crit = "loglik")
             
Coef(fit)
coef(fit, matrix = TRUE)
summary(fit)

#------------------------------------------------------------------------#
# An intercept-only model. Two responses.
fit1 <- vglm(cbind(y2, y2) ~ 1,   # cbind(y1, y2)
             gen.betaIImr(lss = TRUE),  
             data = gdata, trace = TRUE, crit = "loglik")
             
Coef(fit1)
coef(fit1, matrix = TRUE)
summary(fit1)
vcov(fit1, untransform = TRUE)

#------------------------------------------------------------------------#
# An example incorporating one covariate. Constraints are set accordingly.
# x2 affects shape1.a and shape3.q.
# Note that the first option uses 'constraints', whilst in the second
# choice we use the argument 'zero' to 'set' the same constraints.

### Option 1.
c1 <- rbind(0, 1, 0, 0)
c2 <- rbind(0, 0, 0, 1)
mycons <- matrix( c(c1, c2), nc = 2, byrow = FALSE)

fit2 <- vglm(y1 ~ x2, gen.betaIImr(lss = TRUE, zero = NULL),
             data = gdata, trace = TRUE, crit = "loglik",
             constraints = list(x2 = mycons ))

coef(fit2, matrix = TRUE)
summary(fit2)
vcov(fit2)
constraints(fit2)


### Option 2.
fit3 <- vglm(y1 ~ x2, 
             gen.betaIImr(lss = TRUE, 
                          zero = c("scale", "shape2.p")),
             data = gdata, trace = TRUE, crit = "loglik")

coef(fit3, matrix = TRUE)
summary(fit3)
vcov(fit3)
constraints(fit3)

  %}
}
\keyword{Singh-Maddala}
\keyword{Dagum}
\keyword{Inverse Weibull}
