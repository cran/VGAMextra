# These functions are
# Copyright (C) 2014-2018 V. Miranda and T. W. Yee, University of Auckland

isNA <- function(x) is.na(x) & !is.nan(x)

inspectVGAMextra <- function(x, inverse = FALSE, the.NR = NULL,
                             na.s = TRUE, nan.s  = TRUE, 
                             b.valueG = NULL, b.valueL = NULL, 
                             inf.s = NULL, neginf = NULL, 
                             extra = list(NULL)) {
  
  if (!is.list(extra))
    stop("Argument 'extra' must be a list.")
  
  if (length(b.valueG) && !is.Numeric(b.valueG))
     stop("Invalid input for argument 'b.valueG'")
  
  if (length(b.valueL) && !is.Numeric(b.valueL))
    stop("Invalid input for argument 'b.valL'")
  
  if (!inverse) {
    #if (length(b.valueG))
    thx.Gb  <- which(x >= b.valueG)
    
    #if (length(b.valueL))
    thx.Lb  <- which(x <= b.valueL)
    
    #if (na.s)
    thxNA  <- which(isNA(x))
    #if (nan.s)
    thxNaN <- which(is.nan(x))
    
    #if (inf.s)
    thxinf <- which(is.infinite(x) & (x > 0) )
    thxNif <- which(is.infinite(x) & (x < 0) )
    if (length(c(thx.Gb, thx.Lb, thxNA, thxNaN, thxinf, thxNif))) 
      x <- x[-c(thx.Gb, thx.Lb, thxNA, thxNaN, thxinf, thxNif)]
    
    return(list( x = x, list(x.bvalG = thx.Gb, x.bvalL = thx.Lb, 
                 thxNA = thxNA, thxNaN = thxNaN,
                 thxinf = thxinf, thxNif = thxNif,
                 extra = NULL)))
    
  } else {
    
    if (length(extra$x.bvalG) && !length(b.valueG))
      warning("Argument 'b.valueG' might be NULL")
    
    if (length(extra$x.bvalL) && !length(b.valueL))
      warning("Argument 'b.valueL' might be NULL")
    
    if (length(extra$thxinf) && !length(inf.s))
      warning("Argument 'inf.s' for replacement might be NULL")
    
    if (length(extra$thxNif) && !length(neginf))
      warning("Argument 'neginf' for replacement might be NULL")
    
    my.dim <- dim(x); dim(x) <- NULL
    
    if(length(extra$x.bvalG))
      x[extra$x.bvalG] <- b.valueG 
    
    if(length(extra$x.bvalL))
      x[extra$x.bvalG] <- b.valueL
    
    if(length(extra$thxNA))
      x[extra$thxNA] <- NA
    
    if(length(extra$thxNaN))
      x[extra$thxNaN] <- NaN
    
    if(length(extra$thxinf))
      x[extra$thxinf] <- inf.s
    
    if(length(extra$thxNif))
      x[extra$thxNif] <- neginf
    
    
    if (length(unlist(extra))) {
      big.v <- unlist(extra)
        x[-big.v] <- if (length(the.NR)) the.NR else x[-big.v]
    } else {
       x <- the.NR
    }

    dim(x) <- my.dim
    
    x
  }
}