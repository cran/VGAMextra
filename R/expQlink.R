##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.



expQlink <- function(theta, p = stop("Argument 'p' must be entered."),
                     bvalue = NULL, inverse = FALSE,
                     deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for argument 'p'.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    
    e.string <- if (short) paste("expQlink(", theta, ")", sep = "") else
      paste("log(1 - ", p, ")^(-1 / ",
            as.char.expression(theta), 
            ")", sep = "")
    
    if (tag) 
      e.string <- paste("Exponential quantile link:", e.string)
    
    return(e.string)
  } 
  
  if (length(bvalue))
    theta[theta <= 0] <- bvalue
  
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
           -log1p(-p) / theta ,
           theta^2 / log1p(-p),
           theta^3 / log1p(-p)^2)
    
  } else{
    
    switch(deriv + 1,
           -log1p(-p) / theta,
           log1p(-p) / theta^2,
           -log1p(-p) / theta^3)
  }
}
