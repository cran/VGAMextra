\name{trunclognormal}
\alias{trunclognormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

Truncated Log-normal Distribution Family Function
}
\description{
Maximum likelihood estimate of the two--parameter lognormal
distribution with lower/upper truncation.
}
\usage{
 trunclognormal(lmeanlog = "identitylink", lsdlog = "loglink",
                min.support = 1e-6, max.support = Inf, zero = "sdlog") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmeanlog, lsdlog, zero}{
  Same as \code{\link[VGAM:lognormal]{lognormal}}.
}
  \item{min.support, max.support}{
  Positive lower and upper truncation limits (recycled).
  \code{min.support} enables LHS truncation;
  %(default is\eqn{1e^{-6}})
  \code{max.support} enables RHS truncation (default is 
  none).
}
}
\details{
MLE of the two--parameter (univariate) lognormal distribution subject 
to lower/upper truncation. All response values are greater
than \code{min.support} and lower than \code{max.support}.

Default values of \code{min.support}, \code{max.suppport} should 
effectively reproduce \code{\link[VGAM:lognormal]{lognormal}}.

The truncated--lognormal density for a response \eqn{Y} is


\deqn{f(y; \mu, \sigma) =    f_N(y; \mu, \sigma) /
  [\Phi(\texttt{max.support}, \mu, \sigma) -
  \Phi(\texttt{min.support},\mu, \sigma) ], }{
  f(y; \mu, \sigma) = f_N(y; \mu, \sigma) /
  [\Phi(\texttt{max.support}, \mu, \sigma) -
  \Phi(\texttt{min.support},\mu, \sigma) ],
}

where \eqn{f_N} is the ordinary lognormal density (see
\code{\link[VGAM:lognormal]{lognormal}}) and 
\eqn{\Phi} is the standard normal CDF.


The mean of Y, given by
\deqn{
\exp{(\mu + \sigma^2/2)} \cdot
[\Phi(((\log(\texttt{max.support}) - \mu)/\sigma) - \sigma) -
 \Phi(((\log(\texttt{min.support}) - \mu)/\sigma) - \sigma)  ]  /  
\Delta \Phi(\mu,\sigma),}{
\exp{(\mu + \sigma^2/2)} \cdot
[\Phi(((\log(\code{max.support}) - \mu)/\sigma) - \sigma) -
 \Phi(((\log(\code{min.support}) - \mu)/\sigma) - \sigma)  ]  /  
\Delta \Phi(\mu,\sigma), 
} 

with \eqn{\Delta \Phi(\mu, \sigma) = 
 \Phi( (\log(\texttt{max.support}) - \mu)/\sigma )  -
 \Phi( (\log(\texttt{min.support}) - \mu)/\sigma ),}{
 \Delta \Phi(\mu, \sigma) = 
 \Phi( (\log(\texttt{max.support}) - \mu)/\sigma )  -
 \Phi( (\log(\texttt{min.support}) - \mu)/\sigma ),
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
Victor Miranda, Siqi (Vicky) Liu and Thomas W. Yee.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[VGAM:lognormal]{lognormal}},
  \code{\link[VGAM:uninormal]{uninormal}},
  \code{\link[VGAM:CommonVGAMffArguments]{CommonVGAMffArguments}},
  \code{\link[stats:Lognormal]{Lognormal}}.
}
\examples{
##########
set.seed(10470923)
nn <- 3000

## Parameters
mysdlog   <- exp(-1.5)   # sdlog
LL   <- 3.5              # Lower bound
UL   <- 8.0              # Upper bound

## Truncated data
ldata2 <- data.frame(x2 = runif(nn))
ldata2 <- transform(ldata2, y1 = rtrunclnorm(nn, 1 + 1.5 * x2, mysdlog, 
                                             min.support = LL, max.support = UL))
# head(ldata2)
# hist(ldata2$y1, breaks = 22, col = "blue", xlim = c(0, 10))

##############################################################
# Fitting a truncated lognormal distribution - sd is intercept only
fit1 <- vglm(y1 ~ x2, trunclognormal(zero = "sdlog", min.support = LL, max.support = UL),
             data = ldata2, trace = TRUE)
coef(fit1, matrix = TRUE)
vcov(fit1)
             
##############################################################
# Fitting a truncated lognormal distribution - zero = NULL
fit2 <- vglm(y1 ~ x2, trunclognormal(zero = NULL, min.support = LL, max.support = UL),
             data = ldata2, trace = TRUE)
coef(fit2, matrix = TRUE)
vcov(fit2)

##############################################################
# Mimicking lognormal()
fit3 <- vglm(y1 ~ x2, trunclognormal(zero = "sdlog"),
             data = ldata2, trace = TRUE)
coef(fit3, mat = TRUE)

# Same as
fit3bis <- vglm(y1 ~ x2, lognormal(zero = "sdlog"),
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
