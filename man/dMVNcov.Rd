\name{dmultinorm}
\alias{dmultinorm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density for the multivariate Normal distribution}
\description{
Density for the multivariate Normal distribution

 

}
\usage{
      dmultinorm(vec.x, vec.mean = c(0, 0),
                 mat.cov = c(1, 1, 0),
                 log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec.x}{
   For the \eqn{R}--multivariate Normal, an \eqn{R}--vector of 
   quantiles.
   
  }
  \item{vec.mean}{
   The vector of means. 
   
   }
   \item{mat.cov}{
   The vector of variances and covariances, arranged in that order.
   See below for further details.
   
   }
   \item{log}{
   Logical. If \code{TRUE}, the logged values are returned.
   
   }
  }
  \details{
  
  This implementation of the multivariate (say
  \eqn{R}--dimensional) Normal density 
  handles the variances and covariances, instead of 
  correlation parameters. 
  
  For more than one observation, arrange all entries in
  matrices accordingly.
  
  For each observation, \code{mat.cov} is a vector of
  length \eqn{R \times (R + 1) / 2}{R * (R + 1) / 2},
  where the first \eqn{R} entries are the variances
  \eqn{\sigma^2{i}}{\sigma^2[i]}, \eqn{i = 1, \ldots, R}, and 
  then the covariances arranged as per rows, that is,
  \eqn{cov_{ij}}{cov[ij]} \eqn{i = 1, \ldots, R,
  j = i + 1, \ldots, R}. 
  
  By default, it returns the density of two independent
  standard Normal distributions. 
  
  
}
\value{

 The density of the multivariate Normal distribution.
  

}

\section{Warning}{

For observations whose covariance matrix is not positive definite,
\code{NaN} will be returned.

}

\author{

 Victor Miranda
 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{

  \code{\link[VGAM:binormal]{binormal}}.
  
}

\examples{
###
### Two - dimensional Normal density.
###
set.seed(180228)
nn  <- 25
mean1 <- 1; mean2 <- 1.5; mean3 = 2
var1 <- exp(1.5); var2 <- exp(-1.5); var3 <- exp(1); cov12 = 0.75
dmvndata <- rbinorm(nn, mean1 = 1, mean2 = 1.5, var1 = var1, var2 = var2,
                    cov12 = cov12)

## Using dbinorm() from VGAM.
d2norm.data <- dbinorm(x1 = dmvndata[, 1], x2 = dmvndata[, 2],
                        mean1 = mean1, mean2 = mean2, var1 = var1, var2 = var2,
                        cov12 = cov12)
## Using dmultinorm().
d2norm.data2 <- dmultinorm(vec.x = dmvndata, vec.mean = c(mean1, mean2),
                        mat.cov = c(var1, var2, cov12))
summary(d2norm.data)
summary(d2norm.data2)
##
## 3--dimensional Normal.
##
dmvndata <- cbind(dmvndata, rnorm(nn, mean3, sqrt(var3)))

d2norm.data3 <- dmultinorm(dmvndata, vec.mean = c(mean1, mean2, mean3),
                       mat.cov = c(var1, var2, var3, cov12, 0, 0))

hist(d2norm.data3)
summary(d2norm.data3)

}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
