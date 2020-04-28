##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

cm.ARMA <- function(Model   = ~ 1,
                    Resp    =  1,
                    lags.cm =  2,     # It's the order
                    offset  = -2,     # Default on theta2
                    whichCoeff = 1,
                    factorSeq  = 2) { 
  
  if (!Is.Numeric(Resp, 
                  isInteger = TRUE,
                  Nnegative = TRUE))
    stop("Entries for argument 'Resp' must be positive integers.")
  
  if(!is.FormulaAR(Model, Resp = Resp))
    stop("Wrong input for argument 'Model'.")
  
  Lag <- lags.cm; rm(lags.cm)
  if (!Is.Numeric(Lag, 
                  isInteger = TRUE,
                  Nnegative = TRUE))
    stop("Entries for argument 'lags.cm' must be positive integers.")
  
  if (any(Lag < 2))
    stop("Value of 'Lag' is less than 2. No constrained", 
         " matrices to be considered." )
  
  if (!Is.Numeric(offset, 
                  isInteger = TRUE))

  if(abs(max(offset)) < 2 || all(rep(abs(max(offset)), length(Lag)) > Lag))
    stop("Wrong input for argument 'offset'. abs(offset) must be greater",
         "\n",
         " than 1 and less than or equal to 'lags.cm' for each response.")
  
  if(!Is.Numeric(factorSeq, 
                 isInteger = TRUE, 
                 Nnegative = TRUE))
    stop("Entries for argument 'factorSeq' must be positive integers.")
  
  if (!Is.Numeric(whichCoeff, 
                  isInteger = TRUE,
                  Nnegative = TRUE))
    stop("Entries for argument 'whichCoeff' must ne positive integers.")

  if (max(whichCoeff)  >= abs(max(offset)))
    stop("Wrong input for argument 'whichCoeff'. Must be less",
         "\n",
         " than 'abs(offset)' strictly for each response.")

  if (Resp > 1 && (length(Lag) != Resp)) {
    Lag <- rep(Lag, Resp)
    warning("Values of 'lags.cm' (and others) are recycled. Its",
            "\n",
            " lenght differs from number of responses.")
  }
    
  # Offset starts at -2... + mean + sigma positions.
  myMaxp <- which(Lag == max(Lag))[1]
  myMax  <- Lag[myMaxp]
  
  # Matching arguments #
  Lag    <- rep(Lag, Resp)[1:Resp]
  offset <- rep(offset, Resp)[1:Resp]
  whichCoeff <- rep(whichCoeff, Resp)[1:Resp]
  factorSeq  <- rep(factorSeq, Resp)[1:Resp]
  
  # Temporal variables.
  tM1        <- myMax + 2
  twhichCoe  <- whichCoeff[myMaxp]
  toffset    <- abs(offset[myMaxp])
  tfactorSeq <- factorSeq[myMaxp]
  
  # Delete first two columns plus myoffset.
  H <- diag(tM1); H <- H[, -((toffset + 2):tM1)]
  H[, twhichCoe + 2] <- 
    c( rep(0, times = twhichCoe + 1 ),  1 , 
       rep(0, ((toffset + 2) - 1) - (twhichCoe + 2) ), 
       c(tfactorSeq + 0:(tM1 - (toffset + 2)) ) )
  
  diagMat <- diag(Resp)
  Hstar <- kronecker(diagMat, H)
  M1    <- Lag + 2
  maxM1 <- max(M1)
  Haux  <- Hstar
  auxMy <- 0
  
  # Deleting unnecessary ROWS.
  for (ii in seq(length(M1))) {
    if ( M1[ii] < maxM1 )
      Haux <- Haux[-((auxMy + M1[ii] + 1):(auxMy + maxM1)), ] 
    auxMy <- auxMy + M1[ii]
  }
  
  vars.only <- attr(terms(Model), "term.labels")
  interTest <- attr(terms(Model), "intercept")
  lvars     <- length(vars.only) + interTest
  myList    <- vector("list", lvars)
  names(myList) <- 
    c( if (as.logical(interTest)) "(Intercept)" else NULL, vars.only )
  for (jj in 1:lvars) 
    myList[[jj]] <- Haux
  
  myList
}
