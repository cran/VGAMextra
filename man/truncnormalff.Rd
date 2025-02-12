\name{truncnormal}
\alias{truncnormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

Truncated normal Distribution Family Function
}
\description{
Maximum likelihood estimate of the two--parameter normal
distribution with lower/upper truncation.
}
\usage{
 truncnormal(lmean = "identitylink", lsd = "loglink",
             min.support = -Inf, max.support = Inf, zero = "sd") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmean, lsd}{
  Link functions applied to mean and standard deviation/variance. 
}
  \item{min.support, max.support}{
  Vector of lower and upper truncation limits (recycled).
  \code{min.support} enables LHS truncation and
  \code{max.support} enables RHS truncation. 
  The default imply no truncation (mimicks 
  \code{\link[VGAM:uninormal]{uninormal}}).
}
  \item{zero}{
  See \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}} for more information.
}
}
\details{
MLE of the two--parameter (univariate) normal distribution subject to 
lower/upper truncation. All response values are greater
then \code{min.support} and/or lower than \code{max.support}.

The truncated--normal density for a response \eqn{Y} is


\deqn{f(y; \mu, \sigma) =    f(y; \mu, \sigma) /
  [\Phi(\texttt{max.support}, \mu, \sigma) -
  \Phi(\texttt{min.support},\mu, \sigma) ], }{
  f(y; \mu, \sigma) = f_N(y; \mu, \sigma) /
  [\Phi(\texttt{max.support}, \mu, \sigma) -
  \Phi(\texttt{min.support},\mu, \sigma) ],
}

where \eqn{f} is the probability density function of standard normal
distribution and \eqn{\Phi} is the standard normal CDF.


The mean of Y, given by

\deqn{
\mu + [\varphi(\texttt{min.support}) + \varphi(\texttt{max.support})/\Delta \Phi(\mu,\sigma)]\cdot \sigma,
}{
\mu + [\varphi(\texttt{min.support}) + \Phi(\texttt{max.support})/\Delta \Phi(\mu,\sigma)] * \sigma,
}

with \eqn{\Delta \Phi(\mu, \sigma) = 
 \Phi((\texttt{max.support} - \mu)/\sigma )  -
 \Phi( (\texttt{min.support} - \mu)/\sigma ),}{
 \Delta \Phi(\mu, \sigma) = 
 \Phi( (\texttt{max.support} - \mu)/\sigma )  -
 \Phi( (\texttt{min.support} - \mu)/\sigma ),
}
are returned as the fitted values.

}


\value{
  An object of class \code{"vglmff"} 
  (see \code{\link[VGAM:vglmff-class]{vglmff-class}}).
  The object is used by modelling functions such as 
  \code{\link[VGAM:vglm]{vglm}},
  and \code{\link[VGAM:vgam]{vgam}}.
}

\references{
Nadarajah, S. and Kotz, S. (2003).
\code{R} Programs for Computing Truncated Distributions.
\emph{Journal of Statistical Software, Code
Snippets},
\bold{16}(2), 1--8. 


Cohen, A.C. (1991)
\emph{Truncated and Censored Samples:
Theory and Applications},
New York, USA. Marcel Dekker.

}
\author{
Siqi (Vicky) Liu, Victor Miranda, and Thomas W. Yee.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[VGAM:uninormal]{uninormal}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}}.
}
\examples{
nn <- 2000
set.seed(14290909)

## Parameters
mysd <- exp(1.0)   # sd
LL   <- -0.5       # Lower bound
UL   <-  8.0       # Upper bound

## Truncated data
ldata2 <- data.frame(x2 = runif(nn))
ldata2 <- transform(ldata2, y1 = rtruncnorm(nn, 1 + 1.5 * x2, mysd, 
                                min.support = LL, max.support = UL))
# head(ldata2)
# hist(ldata2$y1, breaks = 22, col = "blue", xlim = c(-5, 10))

##############################################################
# Fitting a truncated normal distribution - sd is intercept only
fit1 <- vglm(y1 ~ x2, truncnormal(zero = "sd", min.support = LL, max.support = UL),
             data = ldata2, trace = TRUE)
coef(fit1, matrix = TRUE)
vcov(fit1)
             
##############################################################
# Fitting a truncated lognormal distribution - zero = NULL
fit2 <- vglm(y1 ~ x2, truncnormal(zero = NULL, min.support = LL, max.support = UL),
             data = ldata2, trace = TRUE)
coef(fit2, matrix = TRUE)
vcov(fit2)

##############################################################
# Mimicking uninormal()
fit3 <- vglm(y1 ~ x2, truncnormal(zero = "sd"),
             data = ldata2, trace = TRUE)
coef(fit3, mat = TRUE)

# Same as
fit3bis <- vglm(y1 ~ x2, uninormal(zero = "sd"),
                 data = ldata2, trace = TRUE)
coef(fit3bis, mat = TRUE)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
