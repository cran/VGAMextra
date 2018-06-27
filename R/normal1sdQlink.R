##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.



### Quantile link for the normal1sd family function ###
# 2016/12/08

normal1sdQlink <- function(theta,  
                           mean = stop("Please, enter the fixed 'mean'."),
                           p = stop(" Please, enter argument 'p'."),
                           bvalue  = NULL, inverse = FALSE, deriv = 0,
                           short   = TRUE, tag = FALSE) {
  internal.mean <- mean
  
  if (!is.Numeric(deriv, length.arg = 1, 
                             integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv'out of range.")
    
  if (!is.Numeric(mean))
    stop("Invalid 'mean'.")
  dim(internal.mean) <- NULL #  Manage it as a vector.
  
  
  if (length(mean) != 1) # Else, recycling rule over the mean
    mean <- "mean"       # dimnames(fv) <- list(yn, predictors.names) 
  
  if (!is.Numeric(p, positive = TRUE) || p >= 1)
    stop("Invalid value for 'p'.")
    
  if (is.character(theta)){
    cha.theta <- 
      if (short) paste("normal1sdQlink(", theta, "; ", 
                       p, ", ", mean, ")",sep = "") else
        paste(as.char.expression(mean), " + ", 
              as.char.expression(theta),
              " * sqrt(2) * erf^(-1) (2 * ", p, " - 1)", sep = "")
    if (tag)
      cha.theta <- paste("1-parameter (sd) Normal Quantile Link: ",
                         cha.theta, sep = "")
    return(cha.theta)
  }
  
  if (any(p == 0.5))  # Invlink doesn't exist for p = 0.5
    p[p == 0.5] <- 0.5009
  
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE) 
    } else {
      p <- p[1]
      warning("Taking only the first entry of 'p'. Extend this by enter",
              " 'theta' as a matrix.")
    }
  
  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  #if (length(bvalue) & any(theta <= 0))
  #  theta[theta <= 0] <- bvalue
    
  int.const <- sqrt(2) * erf(x = 2 * p - 1, inverse = TRUE)

  if (inverse) {
    switch(deriv + 1, 
           (theta - internal.mean) / int.const,
           1 / int.const, 0)
  } else {
    switch(deriv + 1,
           internal.mean + theta * int.const,
           int.const, 0)
  }
}