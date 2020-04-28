##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

checkTS.VGAMextra <- function(thetaEst = NULL, # A vector
                              tsclass  = c("AR", "MA"), 
                              chOrder  = 1,
                              NofS     = 1, 
                              retmod   = TRUE,
                              pRoots   = TRUE) {
  
  if ( !length(thetaEst) )
    stop("Please, enter vector 'thetaEst'.")
  
  if(length(thetaEst) && !is.vector(thetaEst))
    stop("Wrong input for argument 'thetaEst'.",
         " It must be a vector.")
  
  if(length(NofS) && !is.Numeric(NofS, 
                                 integer.valued = TRUE, 
                                 positive = TRUE))
    stop("Bad input for argument 'NofS'.")
  
  if(length(chOrder) && !is.Numeric(chOrder, length.arg = 1,
                                    integer.valued = TRUE, 
                                    positive = TRUE))
    stop("Bad input for argument 'chOrder'.")
  
  tsclass <- match.arg(tsclass, c("AR", "MA"))[1]
  MM1 <- chOrder
  myRoots <- matrix(0.0, nrow = MM1, ncol = NofS)
  # Arranging parameters. One column per response.
  myPars <- matrix( thetaEst, ncol = NofS, nrow = MM1,
                    byrow = FALSE)  # ( * * )

  if (tsclass == "AR")
    for (jj in 1:NofS)  {
    checkAux <- polyroot( c( 1, -myPars[ , jj] ) )
    myRoots[seq(length(checkAux)), jj] <- checkAux
    }
  
  if (tsclass == "MA")  
    for (jj in 1:NofS) {
      checkAux <- polyroot( c( 1, myPars[ , jj] ) )
      myRoots[seq(length(checkAux)), jj] <- checkAux
    }

  if (pRoots) {
    # --- A fancy way to display the roots --- #
    dimnames(myRoots) <- 
      list(paste("Root" , 1:MM1 , sep = ""),
           paste("Model", 1:NofS, sep = ""))
    catHelp <- if (tsclass == "AR") 
             "AR component" else "MA component"
        
    cat("\n", 
        "Polynomial roots of the",catHelp, "computed from the given",
        "\n coefficients:",
        "(Examining stationarity/invertibility)",
        "\n", "\n")
    
    # Insert blank spaces and printing out #
    myRootsNB <- myRoots
    myRootsNB[which(Mod(myRoots) == 0 )] <- NA
    show.vanova(x  =  data.frame(Mod(myRootsNB)), 
                digits = .Options$digits)
  } else {
    # Compute and return the Modulus of each root (A matrix).
    if ( retmod  ) Mod(myRoots) else invisible(Mod(myRoots))
  }
}
