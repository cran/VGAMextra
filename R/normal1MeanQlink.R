##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.


### Quantile link for the normal1sd family function ###
# 2018/11/19

normal1MeanQlink <- function(theta,  
                           sd.cons = stop("Please, enter the fixed 'sd'."),
                           p = stop(" Please, enter argument 'p'."),
                           bvalue  = NULL, inverse = FALSE, deriv = 0,
                           short   = TRUE, tag = FALSE) {
  internal.sd <- sd.cons
  
  if (!is.Numeric(deriv, length.arg = 1, 
                             integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv'out of range.")
    
  if (!is.Numeric(sd.cons))
    stop("Invalid 'sd.cons'.")
  dim(internal.sd) <- NULL #  Manage it as a vector.
  
  
  if (length(sd.cons) != 1) # Else, recycling rule over the mean
    sd.cons <- "sd.cons"       # dimnames(fv) <- list(yn, predictors.names) 
  
  if (!is.Numeric(p, positive = TRUE) || any(p >= 1))
    stop("Invalid value for 'p'.")
    
  if (is.character(theta)){
    cha.theta <- 
      if (short) paste("normal1MeanQlink(", theta, "; ", "p = ",
                       p, ", fixed-sd = ", sd.cons, ")",sep = "") else
        paste(as.char.expression(sd.cons), " + ", 
              as.char.expression(theta),
              " * sqrt(2) * erf^(-1) (2 * ", p, " - 1)", sep = "")
    if (tag)
      cha.theta <- paste("1-parameter Normal Quantile Link (models the mean): ",
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
           theta - internal.sd *  int.const,
           1 , 0)
  } else {
    switch(deriv + 1,
           internal.sd + theta * int.const,
           1, 0)
  }
}