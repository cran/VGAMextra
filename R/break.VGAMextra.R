##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.


break.VGAMextra <- function( eta      = NULL,
                             M1       = NULL,
                             noInter  = NULL,
                             bOrder   = NULL,  
                             NOS      = NULL,
                             lInter   = "identitylink",
                             lvar     = "loge",
                             lsd      = "loge",
                             lcoeff1  = "rhobit",   # For Odd positions.
                             lcoeff2  = "rhobit",   # For Even positions.
                             typeTS   = "AR",
                             namesLP  = FALSE,      # If TRUE returns names
                             Complete = FALSE,
                             varArg   = NULL) {
  
  if ( !is.logical ( varArg ))
  stop(" Invalid entry for 'varArg'.")
  
  if ( length(M1) != NOS )
    stop("Conflicting number of responses.")
  
  lInter <- as.list(substitute(lInter))
  eInter <- link2list(lInter)
  lInter <- attr(eInter, "function.name")
  
  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")
  
  lvar <- as.list(substitute(lvar))
  evar <- link2list(lvar)
  lvar <- attr(evar, "function.name")
 
  
  lcoeff <- as.list(substitute(lcoeff1))
  ecoeff <- link2list(lcoeff)
  lcoeff <- attr(ecoeff ,"function.name")
  
  #-----------------------------------------------------------#
  if ( !namesLP ) {
    
    justAsum <- matrix(0.0, nrow = nrow(eta) , ncol = NOS)
    myauX    <- 0; ts.drMean <- ts.var <- ts.sd <- numeric(0)
    ts.The   <- vector("list", length = NOS)
    
    for ( kk in 1:NOS ) {
      ts.drMean <- cbind( ts.drMean , 
                          eta2theta(eta[, myauX + M1[kk] - bOrder[kk] - 1 , 
                                        drop = FALSE],
                                    link = lInter , earg = eInter ))
      
      if ( varArg ) 
        ts.var <- cbind( ts.var , 
                         eta2theta(eta[, myauX + M1[kk] - bOrder[kk], 
                                       drop = FALSE],
                                   link = lvar , earg = evar ) ) else 
        ts.sd  <- cbind( ts.sd , 
                         eta2theta(eta[, myauX + M1[kk] - bOrder[kk], 
                                       drop = FALSE], 
                                   link = lsd , earg = esd ) )
                                   
      for ( jj in 1:bOrder[kk] ) {
        auxLinkinv <- eta2theta(eta[, myauX + M1[kk] - bOrder[kk] + jj , 
                                    drop = FALSE],  lcoeff , earg = ecoeff) 
        ts.The[[kk]] <- cbind( ts.The[[kk]] , auxLinkinv)
        justAsum[, kk] <- justAsum[, kk] + auxLinkinv
      }
                                   
      myauX <- myauX +  M1[kk]
    }
    
    
    if ( varArg )
      ts.sd <- sqrt(ts.var) else ts.var <- ts.sd^2
    
    if ( Complete ) {
      aux4 <- NULL
      max2 <- max( bOrder )[1]
      for ( jj in 1:NOS ) {
        if ( ncol( ts.The[[jj]] ) < max2 ) {
          aux3 <- matrix(0.0, nrow = nrow(ts.The[[jj]]) , 
                         ncol = max2 - ncol( ts.The[[jj]] ))
          ts.The[[jj]] <- cbind( ts.The[[jj]], aux3)
        }  
        aux4 <- cbind( aux4, ts.The[[jj]] )
      }
    }
    
    # if ( Complete ) then ts.The is a matrix,
    # else, it's a list: entry/response!! 
    # i.e. if ( Complete ) the 4th entry is a matrix else 
    # 4th entry is splitted in several entries.
    if (Complete) {
      ts.The <- aux4 
    }
    
  } else {
    
    namesCO    <- vector( "list", NOS )
    mean.names <- vector( "character", NOS)
    my.names   <- vector( "character", NOS)
    OnlyNames  <- vector("character")
    predictors.names <- NULL
    
    if ( typeTS != "AR" && typeTS != "MA" )
      stop("Only names for AR and MA coefficients are",
           " currently handled.")
    
    if (typeTS == "AR")
      auxNam <- "drift"
    if (typeTS == "MA")
      auxNam <- "dev.Mean"
    
    for ( jj in 1:NOS ) {
      mean.names[jj] <- 
        if (NOS == 1) paste(typeTS, auxNam, sep = "") else 
                  paste(paste(typeTS, auxNam, sep = ""), jj , sep = "")
      my.names[jj] <- 
        if ( varArg && NOS == 1) "noiseVar" else
          if ( varArg )
            paste("noiseVar", jj, sep = "") else
              if (NOS == 1) "noiseSD" else
                paste("noiseSD", jj, sep = "") 
      OnlyNames <- c(OnlyNames, 
                     if ( noInter ) NULL else mean.names[jj], my.names[jj])
      
      predictors.names <- 
        c(predictors.names,
          if ( noInter ) NULL else
            namesof(mean.names[jj], link = lInter , 
                    earg = eInter , tag = FALSE),
          namesof(my.names[jj], 
                  if ( varArg ) lvar else lsd , 
                  earg = ifelse( varArg , evar , esd ),
                  tag = FALSE ))    
      
      nameAux <- paste(typeTS, "coeff", sep = "")
      namesCO[[jj]] <- if (NOS == 1)
        paste(nameAux, 1:bOrder[jj], sep = "") else                 
              paste(paste(nameAux,1:bOrder[jj], sep = ""), jj, sep = "")
      OnlyNames <- c(OnlyNames, namesCO[[jj]])
      
      for ( kk in 1:bOrder[jj]) {
        predictors.names <-  c(predictors.names, 
                    namesof(namesCO[[jj]][kk], 
                            link = lcoeff , earg = ecoeff , tag = FALSE))
      }
    }
  }    # End of if else (!namesLP)
  
  if (namesLP)
    list(predictors.names, mean.names, my.names, namesCO, OnlyNames) else
  list( if ( noInter ) NA else ts.drMean, 
        ts.sd, ts.var , ts.The , justAsum)
}
