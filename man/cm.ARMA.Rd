\name{cm.ARMA}
\alias{cm.ARMA}
%\alias{ARXff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Constraint matrices for vector generalized time series
        family functions. }
\description{

 Constraint matrices for coefficients of vector genelized time series 
 family functions in \pkg{VGAMextra}.
 
 }


\usage{
       cm.ARMA(Model      = ~ 1, 
               Resp       =  1,
               lags.cm    =  2,
               offset     = -2, 
               whichCoeff =  1,
               factorSeq  =  2) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{
  A symbolic description of the model being fitted. Must match that
  formula specified in the \code{vglm()} call.
  
  
  }
  \item{lags.cm}{
  Vector of POSITIVE integers greater than 1 indicating the \emph{order} 
  for each response. It must match the \emph{orders} entered in the 
  \code{\link[VGAM:vglm]{vglm}} call. Its default
  value is \code{2}, assuming that a TS process of order greater
  than 1 is being fitted. If \code{lags.cm < 2}, then NO constraints
  are required as only one coefficient (AR, MA or ARMA) is being 
  estimated.
  
  
  }
  \item{offset}{
  Vector of integers specifying the position of the ARMA coefficient 
  at which constraints initiate FOR each response.
  If negative, it is recycled and the absolute value is used. The default 
  value is \code{-2}, which refers to the fourth position on the vector 
  parameter, right after the drift or mean, the white noise sd, and
  the first ARMA coefficient.
  
  Particularly, if only one coefficient is being estimated, i.e, an AR
  or MA process of order-1 is being fitted, then NO restrictions over the
  (unique) coefficient are needed. Consequently, \code{abs(offset < 2)} 
  leads to a message error.
  
  
  }
  \item{whichCoeff}{
  Vector of POSITIVE integers strictly less than `\code{abs(offset)}',
  each entry aplies to each response in the \code{vglm(...)} call.
  This argument allows the user to specify the \emph{unrestricted} 
  coefficient to be considered for constraints. 
  For instance, 
  \code{whichCoeff = 2} means
  that \eqn{\theta_2}{theta[2]} is the required coefficient to 
  compute the constraint matrices. By default, \code{whichCoeff = -1} 
  which implies that \eqn{\theta_1}{theta[1]} is used for this purpose.
  
  If \code{whichCoeff} is greater than or equal to
  \code{abs(offset)}, an error message is displayed since 
  constraints must be function of unrestricted
  parameters. 
  
  
  }
  \item{Resp}{
  The number of responses in the model fitted. Must match the number of 
  responses given in \code{formula} in the \code{vglm} call.
  
  
  }
  \item{factorSeq}{
  Vector of POSITIVE integers. Thus far, restrictions handled are 
  \emph{geometric sequences} and \emph{arithmetic progressions}.
  Hence, \code{factorSeq} specifies either the initial \emph{power} or
  \emph{factor} at restrictions. 
  
  See below for further details.
  
  
  }
}

\details{

NOTE: Except for the \code{Model}, all arguments of length 1 are 
recycled when \code{Resp} \eqn{ \geq 2}{\ge 2}.

Time Series family functions in \pkg{VGAMextra} that are derived from 
AR(p) or MA(q) processes include the \emph{drift} term (or mean) and the 
\emph{white noise} standard deviation as the first two elements
in the vector parameter. For an MA(4), for example, it is given 
by
\deqn{(\mu, \sigma_\varepsilon, \phi_1, \phi_2, \phi_3, \phi_4).}{%
(\mu, \sigma[e], \phi[1], \phi[2], \phi[3], \phi[4]).
}

Thus, constraint matrices on coefficients can be stated from the 
\emph{second} coefficient, i.e., from
\eqn{\phi_2}{\phi[2]}. This feature is specified with 
\code{offset = -2} by default. 

%Note that NO constraints can be set for \eqn{\phi_1}{\phi[1]} only.

In other words, \code{offset} indicates the exact position at 
which parameter restrictions commence. For example, \code{offset = -3} 
indicates that \eqn{\phi_3}{\phi[3]} is the first coefficient over
which constraints are applied. Then, in order to successfully utilize 
this argument, it must be greater than or equal to 2 in absolute value.
Otherwise, an error message will be displayed as no single restriction
are amenable with \eqn{\phi_1}{\phi[1]} only.

Furthermore, if \code{lags.cm = 1}, i.e, a AR or MA process of order one 
is being fitted, then NO constraints are required either, as only one 
coefficient is directly considered.

Hence, the miminum absolute value for argument \code{offset} is 
\code{2}
%At any case, the first two parameters in the
%\emph{vector} parameter are the \emph{drift} or \emph{mean}, and the
%\emph{standard deviation} of the white noise involved. The AR (or MA, 
%or ARMA) coefficients are positioned right after both parameters.
%---------
%Furthermore, \code{cm.ARMA()} currently handles restrictions over the
%AR, MA or ARMA coefficients strictly. 
%Both statements above is particularly important given that the 
%\code{offset} argument considers such property. 

%Hence, if only one coefficient is being estimated, 
%then NO restrictions over the (unique) coefficient are needed. 

As for the \code{factorSeq} argument, its defaul value is 2. 
Let \code{factorSeq = 4}, \code{lags.cm = 5}, \code{offset = -3}, and
\code{whichCoeff = 1}. The coefficient restrictions if a 
\emph{geometric progression} is assumed are

\deqn{ \theta_3 = \theta_1^4, }{\theta[3] = \theta[1]^4, }
\deqn{ \theta_4 = \theta_1^5, }{\theta[4] = \theta[1]^5, } 
\deqn{ \theta_5 = \theta_1^6, }{\theta[5] = \theta[1]^6, }

If coefficient restrictions are in \emph{arithmetic sequence}, 
constraints are given by

\deqn{ \theta_3 = 4 * \theta_1, }{\theta[3] = 4 * \theta[1], }
\deqn{ \theta_4 = 5 * \theta_1, }{\theta[4] = 5 * \theta[1], }
\deqn{ \theta_5 = 6 * \theta_1, }{\theta[5] = 6 * \theta[1], }

The difference lies on thelink function used: 
\code{\link[VGAM:loglink]{loglink}} for the first case, and 
\code{\link[VGAM:identitylink]{identitylink}} for the latter.

Note that conditions above are equivalent to test the following 
two Null Hypotheses:
\deqn{Ho: \theta_k = \theta_1^k}{Ho: \theta[k] = \theta[1]^k}
or
\deqn{Ho: \theta_k = j * \theta_1}{Ho: \theta[k] = k * \theta[1] }
for \eqn{k = 3, 4, 5}. 

Simpler hypotheses can be tested by properly setting all arguments
in \code{cm.ARMA()}.
For instance, the default list of constraint matrices returned by 
\code{cm.ARMA()} allows to test
\deqn{Ho: \theta_k = \theta_1^j}{H: \theta[k] = \theta[1]^k}
for \eqn{k = 2}, in a TS model of order-2 with one response.


}

\value{

A list of constraint matrices with specific restrictions over
the AR(\eqn{p}), MA(\eqn{q}) or ARMA (\eqn{p, q}) coefficients. 
Each matrix returned is conformable with the VGAM/VGLM framework.

Paragrpah above means that each constraint matrix returned by 
\code{cm.ARMA()} is full-rank with \emph{M} rows (number of parameters), 
as required by  \pkg{VGAM}. Note that constraint matrices within 
the VGAM/VGLM framework are \emph{M} by \emph{M} identity matrices 
by default. 

Restrictions currently handled by \code{cm.ARMA()} are (increasing)
arithmetic and geometric progressions.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

}

\section{Warning}{

Hypotheses above can be tested by properly applying
\emph{parameter link functions}. If the test
\deqn{Ho: \theta_k = \theta_1^k,}{Ho: \theta[k] = \theta[1]^k,}
arises, then constraint matrices returned by \code{cm.ARMA()} are
conformable to the use of \code{\link[VGAM:loglink]{loglink}}.

On the other hand, the following hypothesis
\deqn{Ho: \theta_k = k * \theta_1,}{Ho: \theta[k] = k * \theta[1] ,}
properly adapts to the link function 
\code{\link[VGAM:identitylink]{identitylink}}. \eqn{k = 2, 3,\ ldots}.

For further details on parameter link functions within \pkg{VGAM}, see
\code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.


}

\references{

Yee, T. W. and Hastie, T. J. (2003) 
Reduced-rank vector generalized linear models. 
\emph{Statistical Modelling}, \bold{3}, 15--41.

%Yee, T. W. and Wild, C. J. (1996)
%Vector Generalized Additive Models.
%\emph{Journal of the Royal Statistical Society},
%Series B, Methodological, \bold{58(3)}, 481--493.

Yee, T. W. (2008) 
The \code{VGAM} Package. 
\emph{R News}, \bold{8}, 28--39.

}

\author{

 Victor Miranda and T. W. Yee

}

\note{

\code{cm.ARMA()} can be utilized to compute constraint matrices 
for many VGLTSM fmaily functions, e.g.,
\code{\link[VGAMextra:ARXff]{ARXff}} and 
\code{\link[VGAMextra:ARXff]{MAXff}} in \pkg{VGAMextra}.

More improvements such as restrictions on the
\emph{drift parameter} and \emph{white noise standard deviation}
will be set later.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

  \code{\link[VGAM:loglink]{loglink}},
  \code{\link[VGAM:rhobitlink]{rhobitlink}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
  %\code{\link[stats:arima]{arima}}.

}

\examples{

#############
# Example 1.
#############
# Constraint matrices for a TS family function (AR or MA) 
# with 6 lagged terms.
# Restriction commences at the third position (theta[3]) powered to
# or multiplied by 4. Intercept-only model.
position   <- -3
numberLags <-  6
myfactor   <-  4
cm.ARMA(offset = position, lags.cm = numberLags, factorSeq = myfactor)

# With one covariate
cm.ARMA(Model =  ~ x2, offset = position, 
        lags.cm = numberLags, factorSeq = myfactor)

\donttest{
# Or 2 responses...
cm.ARMA(offset = position, lags.cm = numberLags, 
        factorSeq = myfactor, Resp = 2)


# The following call causes an ERROR.
# cm.ARMA(offset = -1, lags.cm = 6, factorSeq = 2)
}

##############
# Example 2.
##############

# In this example, the use of constraints via 'cm.ARMA()' is
# included in the 'vglm' call. Here, two AR(2) models are fitted
# in the same call (i.e. two responses), where different constraints
# are set, as follows:
# a) list(ar = c(theta1, theta1^2)) and
# b) list(ar = c(theta2, theta2^2 )).

# 2.0 Generate the data.
set.seed(1001)
nn     <- 100
# A single covariate.
covdata <- data.frame(x2 =  runif(nn)) 

theta1 <- 0.40; theta2 <- 0.55
drift  <- c(0.5, 0.75)
sdAR   <- c(sqrt(2.5), sqrt(2.0))

# Generate AR sequences, TS1 and TS2, considering Gaussian white noise

# Save both in a data.frame object: the data.
tsdata  <- 
  data.frame(covdata, # Not used 
             TS1 = arima.sim(nn, 
                             model = list(ar = c(theta1, theta1^2)), 
                             rand.gen = rnorm, 
                             mean = drift[1], sd = sdAR[1]),
             TS2 = arima.sim(nn, 
                             model = list(ar = c(theta2, theta2^2)), 
                             rand.gen = rnorm, 
                             mean = drift[2], sd = sdAR[2]))

# 2.1 Fitting both time series with 'ARXff'... multiple responses case.
fit1 <- vglm(cbind(TS1, TS2) ~ 1, 
             ARXff(order = c(2, 2), type.EIM = "exact"), 
             data = tsdata,  
             trace = TRUE)

Coef(fit1)                
coef(fit1, matrix = TRUE)
summary(fit1)

## Same length for both vectors, i.e. no constraints.
length(Coef(fit1))
length(coef(fit1, matrix = TRUE))




###2.2 Now, fit the same models with suitable constraints via 'cm.ARMA()'
# Most importantly, "loglink" is used as link function to adequately match 
# the relationship between coefficients and constraints. That is:
# theta2 = theta1^2, then log(theta2) = 2 * log(theta1).

fit2 <- vglm(cbind(TS1, TS2) ~ 1, 
             ARXff(order = c(2, 2), type.EIM = "exact", lARcoeff = "loglink"), 
             constraints = cm.ARMA(Model = ~ 1, 
                                   Resp = 2,
                                   lags.cm = c(2, 2),
                                   offset  = -2),
             data = tsdata,  
             trace = TRUE)
Coef(fit2)
coef(fit2, matrix = TRUE)
summary(fit2)

# NOTE, for model 1, Coeff2 = Coeff1^2, then log(Coeff2) = 2 * log(Coeff1)
( mycoef <- coef(fit2, matrix = TRUE)[c(3, 4)] )
2 * mycoef[1] - mycoef[2]    # SHOULD BE ZERO

# Ditto for model 2: 
( mycoef <- coef(fit2, matrix = TRUE)[c(7, 8)] )
2 * mycoef[1] - mycoef[2]    # SHOULD BE ZERO

## Different lengths, due to constraints
length(Coef(fit2))
length(coef(fit2, matrix = TRUE))

 




}



%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
