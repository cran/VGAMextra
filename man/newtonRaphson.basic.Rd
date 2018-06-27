\name{newtonRaphson.basic}
\alias{newtonRaphson.basic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Newton--Raphson algorithm
  
}
\description{

  Newton--Raphson algorithm to approximate the roots of
  univariate real--valued functions.
  
  This function is vectorized.
}
\usage{
%
newtonRaphson.basic(f, fprime, a, b, 
                    tol = 1e-8, n.Seq = 20,
                    nmax = 15, ...)
%
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
   A univariate function whose root(s) are approximated.
   This is the target function. Must return a vector.

  }
  \item{fprime}{
   A function. The first derivative of \code{f}. Must return a vector.
   
  }
  \item{a, b}{
    Numeric vectors.
    Upper and lower real limits of the open interval \eqn{(a, b)} 
    where the root(s) of  \code{f} will be searched. Notice, entries 
    \code{Inf}, \code{-Inf}, \code{NA} and \code{NaN} are not handled.
    
    These vectors are subject to be recycled if \code{a} and \code{b}
    lenghts differ.
  
  }
  \item{tol}{
   Numeric. A number close to zero to test whether the 
   approximate roots from iterations \eqn{k} and \eqn{(k + 1)}
   are close enough to stop the algorithm.
  
  }
  \item{n.Seq}{
  Numeric. The number of equally spaced initial points within
  the interval (\code{a}, \code{b}) to internally set up initial 
  values for the algorithm.
  
  }
  \item{nmax}{
   Maximum number of iterations. Default is \eqn{15}.
   
  }
  \item{...}{
  Any other argument passed down to functions \code{f} and
  \code{fprime}.  
  
  }
}


\details{

 This is an implementation of the well--known Newton--Raphson 
 algorithm to find a real root, \eqn{r}{r}, \eqn{a < r < b},
 of the function \eqn{f}.
 
 Initial values, \eqn{r_0}{r[0]} say, for the algorithm are 
 internally computed by drawing `\code{n.Seq}' equally spaced points 
 in \eqn{(a, b)}. Then, the function \code{f} is evaluated at this 
 sequence. Finally, \eqn{r_0}{r[0]} results from the closest image to 
 the horizontal axis.
 
 At iteration \eqn{k}, the \eqn{(k + 1)^{th}}{(k + 1)^th} approximation 
 given by
 %
 \deqn{r^{(k + 1)} = r^{(k)} - 
           {\tt{f}}(r^{(k), ...)} / {\tt{fprime}}(r^{(k)}, ...)}{
       r[k + 1]    = r[k]  - f(r[k], ... ) / fprime(r[k], ...)
 }
 %
 is computed, unless the approximate root from step \eqn{k} is the 
 desired one.
 
 \code{newtonRaphson.basic} approximates this root up to
 a relative error less than \code{tol}. That is, at each iteration, 
 the relative error between the estimated roots from iterations 
 \eqn{k} and \eqn{k + 1} is calculated and then compared to \code{tol}.
 The algorithm stops when this condition is met.
 
 Instead of being single real values, arguments \code{a} and \code{b} 
 can be entered as vectors of length \eqn{n}, say
 \eqn{{\tt{a}} = c(a_1, a_2, \ldots, a_n)}{
                 a = c(a[1], a[2], \ldots, a[n])
    } and
    \eqn{{\tt{b}} = c(b_1, b_2,\ldots, b_n)}{
         b = c(b[1], b[2],\ldots, b[n])
    }. 
 In such cases, this function approaches the (supposed) root(s)
 at each interval \eqn{(a_j, b_j)}{(a[ j ], b[ j ])}, 
 \eqn{j = 1, \ldots, n}. Here, initial values are searched 
 for each interval \eqn{(a_j, b_j)}{(a[ j ], b[ j ])}.
 

}
\value{
 The approximate roots in the intervals  
 \eqn{(a_j, b_j)}{(a[ j ], b[ j ])}.
 When \eqn{j = 1}, then a single estimated root is returned, if any.

}
\author{
   V. Miranda.
   
}
\note{

 The explicit forms of the target function \code{f} and its 
 first derivative \code{fprime} must be available for the algorithm. 
 
 \code{\link[VGAMextra:newtonRaphson.basic]{newtonRaphson.basic}} 
 does not handle yet numerically approximated derivatives.
 
 A warning is displayed if no roots are found, or if more than one 
 root might be lying in
 \eqn{(a_j, b_j)}{(a[ j ], b[ j ])}, for any \eqn{j = 1, \ldots, n}.
 
 If \code{a} and \code{b} lengths differ, then the recyling rule 
 is applied. Specifically, the vector with minimum length
 will be extended up to match the maximum length by repeating
 its values. 
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[VGAM:bisection.basic]{bisection.basic}}
}
\examples{
# Find the roots in c(-0.5, 0.8), c(0.6, 1.2) and c(1.3, 4.1) for the
# f(x) = x * (x - 1) * (x - 2). Roots: r1 = 0, and r2 = 1, r3 = 2.

f <- function(x) x * (x - 1) * (x - 2)
fprime <- function(x) 3 * x^2 - 6 * x + 2

# Three roots.
newtonRaphson.basic(f = f, fprime  = fprime, 
                    a = c(-0.5, 0.6, 1.3), 
                    b = c(0.8, 1.2, 4.1))              ## 0.0, 1.0 and 2.0
                    
# Recycling rule. Intervals analysed are (-0.5, 1.2) and (0.6, 1.2)
newtonRaphson.basic(f = f, fprime  = fprime, 
                    a = c(-0.5, 0.6), b = c(1.2)) 

## Warning: There is more than one root in (-0.5, 1.2)!
}