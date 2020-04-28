##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0


expQlink <- function(theta, p = stop("Argument 'p' must be entered."),
                     bvalue = NULL, inverse = FALSE,
                     deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for argument 'p'.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    e.string <- if (short) paste("expQlink(", theta, "; ", p,
                                 ")", sep = "") else
      paste("logloglink[(1 - ", p, ")^(-1 / ",
            as.char.expression(theta), 
            ")]", sep = "")
    
    if (tag) 
      e.string <- paste("Exponential quantile link:", e.string)
    
    return(e.string)
  } 
  
  if (!inverse)
    theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE) 
    } else {
      p <- p[1]
      warning("Taking only the first entry of 'p'. Extend this by enter",
              " 'theta' as a matrix.")
    }
  
  if (inverse) {
    
    switch(deriv + 1,
           -log1p(-p) / exp(theta),
           -theta,
           theta)
    
  } else{
    
    switch(deriv + 1,
           logloglink((1 - p)^(-1/theta)),
           -1 / theta,
            1 / theta^2)
  }
}
