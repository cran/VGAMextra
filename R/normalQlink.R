#############################################################################
# These functions are Copyright (C) 2014--2020
# V. Miranda-Soberanis, Auckland University of Technology
# T. Yee, University of Auckland

uninormalQlink <- function(tau.arg = stop("Argument 'tau.arg' must be entered."),
                           theta, sd = NULL,
                           inverse = FALSE, wrt.param = NULL, deriv = 0, 
                           bvalue = NULL, short = TRUE, tag = FALSE) {
  
  p <- tau.arg; rm(tau.arg)
  if (length(p) & (!is.Numeric(p, positive = TRUE) || any(p >= 1) ))
    stop("Invalid input for argument 'p'.")
  
  if (!is.Numeric(deriv, length.arg = 1, integer.valued = TRUE) || deriv > 2)
    stop("Argument 'deriv' unmatched.")
  
  if (length(wrt.param) && (!(wrt.param %in% 1:2) ||
                            !is.Numeric(wrt.param, length.arg = 1)) )
    stop("argument 'wrt.param' should be 1 or 2")
  
  
  if (length(theta) && length(sd)) {
    sd.mat <- cbind(sd)
    theta.mat  <- cbind(theta)
    
    if (any(sd.mat < 0))
      stop("Negative values for the 'sd' argument not allowed.")
    
    if (nrow(sd.mat) != nrow(theta.mat))
      stop("Number of rows of 'theta' and 'sd' do not match.")
  }
  
  
  if (is.character(theta)) {
    n.string <- if (short) {
      paste("uninormalQlink(", theta, ", sigma; ", p, ")", sep = "") 
      }  else {
        paste(as.char.expression(theta), 
              " + sd * sqrt(2) * inv-erf( 2 * ",p, "- 1 )", sep = "")
      }
    if (tag)
      n.string <- paste("Normal quantile link:", n.string)
    
    return(n.string)
  }
  
  if (length(p) > 1)
    if (is.matrix(theta)) {
      p <- matrix(p, nrow = nrow(theta), ncol = ncol(theta), byrow = TRUE)
    } else {
      p <- p[1]
      warning("Unmatched dimensions. Taking the first entry of 'p'
              since 'theta' is not a matrix ")
    }
  
  
  ## For this link, theta is 'mu' ideally spanning R.
  #if (!inverse)
  #  theta[theta <= 0] <- if (length(bvalue)) bvalue else NaN
  
  if (inverse) {
    
    theta.ret <- switch(deriv  + 1,
                        theta - sd * sqrt(2) * erf(2 * p - 1, inverse = TRUE),
                        
                        # wrt = 1 -> wrt eta1 ELSE  wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          
                          #  c(d mu / d eta1 , d sd / d eta1 )
                          myvec <- c(rep(1, NCOL(sd.mat)), 
                                     rep(0, NCOL(sd.mat)))
                          matrix(myvec, nrow = nrow(sd.mat),
                                 ncol = 2 * ncol(sd.mat), byrow = TRUE)
                          
                        } else {
                          
                          #  c(d mu / d eta2 , d sd / d eta2 )
                          # Strictly: d sd/deta2 not needed...computed by
                          # the current link function
                          # d sd / deta2 not required. Computed by deta.dtheta
                        cbind(matrix(0, nrow = nrow(sd.mat), 
                                     ncol = ncol(sd.mat)), NULL)
                          
                          #Should be
                          #cbind(matrix(0, nrow = nrow(sd.mat), 
                          #             ncol = nrow(sd.mat)), theta)
                          
                        },
                        
                        # wrt = 1 -> wrt eta1 ELSE  wrt = 2 -> wrt eta2
                        if (wrt.param == 1) {
                          
                          #  c(d2mu / d2eta1 , d2sd / d2eta1 )
                          matrix(0, nrow = nrow(sd.mat), ncol = 2)
                          
                        } else {
                          
                          #  c(d2mu / d2eta2 , d2sd / d2eta2 )
                          # Strictly: d2sd/d2eta2 not needed...computed by
                          # the current link function
                          cbind(matrix(0, nrow = nrow(sd.mat), 
                                       ncol = ncol(sd.mat)), NULL)
                          #Should be
                          # cbind(matrix(0, nrow = nrow(sd.mat), 
                          #             ncol = nrow(sd.mat)), theta)
                        })
    
    eta.ret <- theta.ret
  } else {
    ## Here, eta
    ##       d eta / d thetaj, j = wrt = 1, 2 (mu or sd)
    ##       d2eta / dthetaj2, j = wrt = 1, 2.
    
    mymat <- matrix(c(1, 0), nrow = nrow(sd.mat), ncol = 2,
                    byrow = TRUE)
    
    eta.ret <- switch(deriv + 1,
                      theta + sd * sqrt(2) * erf(2 * p - 1, inverse = TRUE),
                      #  wrt = 1 -> wrt mu (= theta)
                      #  wrt = 2 -> wrt sd
                      if (wrt.param == 1) mymat[, 1, drop = FALSE] else
                                                 mymat[, 2, drop = FALSE],
                      mymat[, 2, drop = FALSE])
  }
  
  eta.ret

}

