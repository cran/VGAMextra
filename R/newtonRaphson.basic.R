##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

# newtonRaphson.basic() vectorized, 20160909.
newtonRaphson.basic <- function(f, fprime, a, b, 
                                tol = 1e-8,
                                n.Seq = 20,
                                nmax = 15, ...) {
  
  if (!all(sign(f(a, ...)) * sign(f(b, ...)) <= 0)) {
    warning("Either no roots or more than one root between the \n",
            "  given limits 'a' and 'b'.")
  }
  
  my.grid <- apply(cbind(a, b), 1, 
                   function(x) seq(x[1], x[2], length.out = n.Seq))
  if (ncol(my.grid) == 1) {
    my.init <- apply(my.grid, 2, function(x) abs(f(x, ...)))
    my.init <- which.min(my.init)
  } else {
    my.init <- apply(my.grid, 1, function(x) abs(f(x, ...)))
    my.init <- apply(t(my.init), 2, which.min)
  }
  x.current <- if (length(my.init) == 1) my.grid[my.init, ] else 
                                          diag(my.grid[my.init, ])
  if (all( f(x.current, ...) == 0)) {
    return(x.current)
  }
  
  N    <- 1
  flag <- FALSE
  
  while(N <= nmax) {
    
    if ( !(all(x.current > a) & all (x.current < b)) ) {
      flag <- TRUE
      x.current[which( x.current < a )] <- a[which( x.current < a )] + tol
      x.current[which( x.current > b )] <- b[which( x.current > b )] - tol
    }
    
    fprime.x <- fprime(x.current, ...)
    x.next   <- x.current - f(x.current, ...) / fprime.x
    xcu.xne <- t(cbind(x.current, x.next))
     
    if (all(dist(xcu.xne) /(tol + sqrt(sum(x.current^2))) < tol)) {  
      return(x.next)
    } 
    
    N <- N + 1
    x.current <- x.next
  }
  
  warning("Newton-Raphson did not converged after ", nmax, " iterations.",
          " Returning the approximated root from the latest iteration.\n")
  return(x.current)
  
}

