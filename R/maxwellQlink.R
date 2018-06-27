# These functions are
# Copyright (C) 2014-2018 V. Miranda and T. W. Yee, University of Auckland

maxwellQlink <- function(theta,
                         p = stop("Argument 'p' must be specified."),
                         bvalue = NULL, inverse = FALSE,
                         deriv = 0, short = TRUE, tag = FALSE) {
  
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for argument 'p'. Must lie between 0 and 1.")
  
  if (!is.Numeric(deriv, length.arg = 1, 
                  integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (is.character(theta)) {
    m.string <- 
      if (short) paste("maxwellQlink(", 
                       theta, "; ", p, ")", sep = "") else
        paste("(2 * qgamma(", p, ", 1.5) / ", 
              as.char.expression(theta), 
              ")^0.5",  sep = "")
    
    if (tag)
      m.string <- paste("Maxwell quantile link: ", m.string, sep = "")
    
    return(m.string)
  }
  
  if (length(bvalue))
    theta[theta <= 0] <- bvalue
  
  if (any(theta <= 0)) 
    theta[theta <= 0] <- NaN
  
  #if (length(p) > 1)
  #  p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE)
  
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
           2 * qgamma(p = p, shape = 1.5, log.p = FALSE) / theta^2,
           -2 * theta^(1.5) / 
             sqrt(2 * qgamma(p = p, shape = 1.5, log.p = FALSE)),
           3 * theta^2 / qgamma(p = p, shape = 1.5, log.p = FALSE))
    
  } else {
    
    switch(deriv + 1, 
           sqrt(2 * qgamma(p = p, shape = 1.5)) / theta^(0.5),
           -sqrt(2 * qgamma(p = p, shape = 1.5)) / (2 * theta^(1.5)),
           3 * sqrt(2 * qgamma(p = p, shape = 1.5)) / (4 * theta^(2.5)))
    
  }
}