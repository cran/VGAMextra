##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0
# Supports the Wald, score, and lrt tests (20180209)

## ARIMAXff control function.
ARIMAXff.control <- function(save.weights = TRUE,
                             summary.HDEtest = FALSE,
                             ...) { 
  list(save.weights = save.weights,
       summary.HDEtest = summary.HDEtest,...)
}


ARIMAXff <- 
  function(order     = c(1, 1, 0),
           zero      = c("ARcoeff", "MAcoeff"),
           diffCovs  = TRUE,
           xLag      = 0,
           include.current = FALSE,
           type.EIM  = c("exact", "approximate")[1],
           var.arg   = TRUE,
           nodrift   = FALSE,
           noChecks  = FALSE,
           ldrift    = "identitylink",
           lsd       = "loglink",
           lvar      = "loglink",
           lARcoeff  = "identitylink",
           lMAcoeff  = "identitylink",
           idrift    = NULL,        # Must be a vector of length NOS
           isd       = NULL,        # Must be a vector of length NOS
           ivar      = NULL,        # Must be a vector of length NOS
           iARcoeff  = NULL,        # Must be a vectors
           iMAcoeff  = NULL) {
    
    
    if (xLag < 0)
      stop("Wrong input for argument 'xLag'")
    
    if ( !Is.Numeric(order, length.arg = 3)  )
      stop("Invalid 'order'. Should be a vector (p, d, q), length 3.")
    
    ARord2 <- order[1]
    dOrd   <- order[2]
    MAord2 <- order[3]; rm(order)
    flagAR <- (ARord2 == 0)
    flag.2 <- ((ARord2 == 0) && (MAord2 == 0)) 
    
    if (flag.2)
      stop("Seems like a random walk is modelled.",
              " Refer to uninormal() for this. ")
    
   # if ((ARord2 == 0) && (MAord2 == 0))
   #    stop("Wrong input for order. Enter the AR order and/or MA order.")

    ARord  <- if (ARord2 == 0) 1 else ARord2
    MAord  <- if (MAord2 == 0) 1 else MAord2
    
    if (ARord2 == 0) {
      rem.lat <- grep("ARcoeff", zero)
      zero <- zero[-rem.lat]
    }
    
    if (MAord2 == 0) {
      rem.lat <- grep("MAcoeff", zero)
      zero <- zero[-rem.lat]
    }
    
    if (!is.logical(diffCovs))
      stop("Bad input for argument 'diffCovs'")
 
    if (!is.logical(include.current))
      stop("Bad input for argument 'include.current'")
    
    if (!dOrd)
      stop("Refer to ARMAXff(), ARXff() or MAXff(). ",
           "Only the ARIMA--class handled", "\n  ",
           "by this family function")
    
    nOrder <- ARord2 + MAord2
    realM1 <- 2 + ARord2 + MAord2 - nodrift
    
    if (!is.logical(noChecks) && !is.logical (var.arg) &&
          !is.logical(nodrift))
      stop("'noChecks', 'nodrift', and 'var.arg' must be logical")
    
    type.likelihood <- "exact"
    type.EIM <- if (flag.2) "approximate" else
       match.arg(type.EIM, c("exact", "approximate"))[1]
    
    ldrift <- as.list(substitute(ldrift))
    edrift <- link2list(ldrift)
    ldrift <- attr(edrift, "function.name")
    
    lvar <- as.list(substitute(lvar))
    evar <- link2list(lvar)
    lvar <- attr(evar, "function.name")
    
    lsd <- as.list(substitute(lsd))
    esd <- link2list(lsd)
    lsd <- attr(esd, "function.name")
      
    lARcoeff <- as.list(substitute(lARcoeff))
    eARcoeff <- link2list(lARcoeff)
    lARcoeff <- attr(eARcoeff, "function.name")
    
    lMAcoeff <- as.list(substitute(lMAcoeff))
    eMAcoeff <- link2list(lMAcoeff)
    lMAcoeff <- attr(eMAcoeff, "function.name")
    
    pre.blurb1 <- paste("ARcoeff", 1:ARord, sep = "")
    pre.blurb2 <- paste("MAcoeff", 1:MAord, sep = "")
    
    blurb.vec <- 
      c(if (nodrift) NULL else 
        namesof("drift", 
                link = ldrift, earg = edrift, tag = FALSE),
        if (var.arg)
          namesof("noiseVar", link = lvar, earg = evar, tag = FALSE) 
        else 
          namesof("noiseSD", link = lsd, earg = esd, tag = FALSE) )
    
    for (jj in 1:ARord)  
      blurb.vec <-  c( blurb.vec, namesof(pre.blurb1[jj],
                                link = lARcoeff, earg = eARcoeff))
    
    for (kk in 1:MAord)  # tag = FALSE by default.. not needed
      blurb.vec <-  c( blurb.vec, namesof(pre.blurb2[kk],
                                 link = lMAcoeff, earg = eMAcoeff))
    
    blurb.vec2 <- character(0)
    for (kk in 1:realM1) {
      if (kk != realM1)
        aux <- paste(blurb.vec[kk], ", ", sep = "")
      else 
        aux <- paste(blurb.vec[kk], ".", sep = "")
      blurb.vec2 <- c(blurb.vec2, aux)
    }
    
    preAR <- if (ARord == 1) "(ARcoeff1) * Y[t - 1] + "
    preMA <- if (MAord == 1) "(MAcoeff1) * e[t - 1]  "
    
    alo <- 
      new("vgltsmff",
          blurb = c("VGLTSMs: ", nOrder + 1 + !nodrift,
                    "-parameter ARIMAX model ",
                    "of order",
                    "(", ARord2, ",", dOrd, ",", MAord2, ").",
                    "\n", "\n",
                    "LINKS: ",
                    blurb.vec2, "\n",
                "Model: ", "\n" , " Y[t] = ", if (nodrift) NULL else 
                      if (!ARord2) "mu + " else "drift + ",
                    " B^T * X_t + ",
                    if (ARord2) {
                     if (ARord2 == 1) "(ARcoeff1) * Y[t - 1] + " else
                c("(ARcoeff1) * Y[t - 1] +... + (ARcoeff", ARord2,")",
                          " * Y[t -", ARord2,"] + ")
                    },
                    if (MAord2) {
                      if (MAord2 == 1) "(MAcoeff1) * e[t - 1]" else
               c("(MAcoeff1) * e[t - 1] + ... + (MAcoeff", MAord,")",
                          " * e[t - ", MAord,"]") 
                    },
                    if ( !MAord2 ) " e[t] " else 
                    " + e[t].", "\n"),
                    #if ( !ARord2 ) "E[Yt] = mu" else
                   #"E[Yt] = drift / (1 - SUM(ARcoeff))"),
          
          
          
          
   constraints = eval(substitute(expression({
      M1 <- .realM1
      constraints <- cm.zero.VGAM(constraints, x = x, 
                                  zero = .zero , M = M,
                                  # Old : predictors.names...
                                  predictors.names = parameters.names,
                                  M1 = M1)
    }), 
    list( .zero = zero , .realM1 = realM1 ) )),
   
   
   
   
   
    infos = eval(substitute(function(...) {
      
      list(M1        = .realM1 ,
           Order     = .nOrder ,
           Q1        = 1, 
           ldrift    = .ldrift ,
           lsd       = .lsd  ,
           lvar      = .lvar ,
           lARcoeff  = .lARcoeff ,
           lMAcoeff  = .lMAcoeff ,
           edrift    = .edrift ,
           esd       = .esd ,
           evar      = .evar , 
           eARcoeff  = .eARcoeff ,
           eMAcoeff  = .eMAcoeff ,
           ARord     = .ARord ,
           MAord     = .MAord ,
           dOrd      = .dOrd ,
           zero      = .zero ,
           nodrift   = .nodrift ,
           type.EIM  = .type.EIM ,
           expected  = TRUE, 
           multipleResponse = TRUE )
      }, list( .realM1 = realM1 , .nOrder = nOrder,
               .ldrift = ldrift , .lsd = lsd , .lvar = lvar ,
               .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
               .edrift = edrift , .esd = esd , .evar = evar ,
               .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
               .ARord = ARord , .MAord = MAord , .dOrd = dOrd ,
               .nodrift = nodrift , .type.EIM = type.EIM ,
               .zero = zero ))),
   
   
   
   
   
    first = eval(substitute(expression({
      
      extra$y.r <- y
      y   <- cbind(y)
      n   <- NROW(y)
      NOS <- NCOL(y)
      mat.save <- matrix(0, n, .dOrd)
      
      if (NOS != 1)
        stop("Currently this VGLSTM family function handles",
             "univariate responses only.")
      
      if ((NCOL(x) == 1) && 
          (colnames(x) == "(Intercept)") && ( .xLag > 0))
        stop("Set xLag = 0. No covariates entered.")
      
      y.diff <- diff(y, differences = .dOrd )
      
      for (ii in 1:NOS) {
        yaux1 <- y[, ii , drop = FALSE]
         for (jj in 1:( .dOrd )) {
           mat.save[c(1:(n - jj)), jj] <-  diff(yaux1)
           yaux1 <- mat.save[c(1:(n - jj)), jj, drop = FALSE] 
         }
      }
      
      y <- extra$y <- y.diff
      extra$mat.save <- mat.save
      
      difxs <- .diffCovs
      if ((NCOL(x) == 1) && (colnames(x) == "(Intercept)") && difxs) {
        difxs <- FALSE
        warning("No covariates entered. 'diffCovs' set to FALSE")
      }
      
      x.int <- x[, -1, drop = FALSE]
      
      if (difxs) {
        x.diff <-  apply(cbind(x.int), 2, function(x) {
          diff(x, differences = .dOrd )
        })
        colnames(x.diff) <- paste("Diff", colnames(x.int), sep = "")
      } else{
        x.diff <- x.int[ c(1:NROW(y.diff)), , drop = FALSE]
      }
      
      x.matrix <- cbind(1, x.diff)
      colnames(x.matrix) <- c("(Intercept)", colnames(x.diff))
      
      x.temp <- temp3 <- NULL
      if ( .xLag ) {
        for (jj in 1:NCOL(x.diff)) {
          temp1 <-  WN.lags(cbind(x.diff[, jj]), .xLag )
          temp2 <- paste(colnames(x.diff)[jj], "Lag", sep = "")
          colnames(temp1) <- paste(temp2, 1:( .xLag ), sep = "")
          x.temp <- cbind(x.temp, temp1)
        }
        rm(temp1, temp2)
      }
      
      if ( .inc.c ) {
        temp3 <- x[1:NROW(x.matrix), -1, drop = FALSE]
      } 
      x.matrix <- cbind(x.matrix, x.temp, temp3)
      
      list.names <- vector("list", NCOL(x.matrix))
      names(list.names) <- colnames(x.matrix)
      for (ii in 1:(NCOL(x.matrix))) 
        list.names[[ii]] <- ii
      attr(x.matrix, "assign") <- list.names
      
      x <- x.matrix
      w <- w[1:NROW(y)]
      
    }), list( .dOrd = dOrd , .xLag = xLag , .diffCovs = diffCovs ,
              .inc.c = include.current ))),
   
   
   
   
   
    initialize = eval(substitute(expression({
      
      M1 <- .realM1 
      mycheck <- w.y.check(w = w, y = y,
                           Is.positive = FALSE, 
                           ncol.w.max  = Inf,
                           ncol.y.max  = Inf,
                           out.wy      = TRUE,
                           colsyperw   = 1,
                           maximize    = TRUE)
      w   <- mycheck$w
      y   <- extra$y <- mycheck$y
      NOS <- ncol(mycheck$y)
      M   <- NOS * M1; n <- NROW(y)
      nOrder <- .nOrder
      y.r    <- extra$y
      dOrd   <- .dOrd
      
      if (length( .iARcoeff ) || length( .iMAcoeff ))
        quick.check.coeffs(arc = .iARcoeff, mac = .iMAcoeff,
                           NOS = NOS, arOrd = .ARord, maOrd = .MAord)
      
      myIdrift <- rep( .idrift , NOS)[1:NOS]
      myIvar   <- rep( .ivar   , NOS)[1:NOS]
      myIsd    <- rep( .isd    , NOS)[1:NOS]
      iniar    <- if (!length( .iARcoeff )) NULL else
                     matrix( .iARcoeff , NOS, .ARord , byrow = TRUE)
      inima    <- if (!length( .iMAcoeff )) NULL else 
                    matrix( .iMAcoeff , NOS, .MAord , byrow = TRUE)
      
      drift.n <- if (NOS == 1) "drift.mean" else 
                                paste("drift.mean", 1:NOS, sep = "")
      varsd.n <- if ( .var.arg ) if (NOS == 1) "noiseVar" else
                     paste("noiseVar", 1:NOS, sep="")  else
                         if (NOS == 1) "noiseSD" else
                            paste("noiseSD", 1:NOS, sep = "")
      
      ar.n <- ma.n <- character(0)
  for (jj in 1:( .ARord ))
    ar.n <- if (NOS == 1) c(ar.n, paste("ARcoeff", jj, sep = "")) else
       c(ar.n, paste(paste("ARcoeff", jj, sep = ""), 1:NOS, sep = ""))
      
   for (jj in 1:( .MAord )) 
    ma.n <- if (NOS == 1) c(ma.n, paste("MAcoeff", jj, sep = "")) else 
       c(ma.n, paste(paste("MAcoeff", jj, sep = ""), 1:NOS, sep = ""))
      
      parameters.names <- c(if ( .nodrift ) NULL else drift.n, varsd.n,
                            if (!( .ARord2 )) NULL else ar.n,
                            if (!( .MAord2 )) NULL else ma.n)
      
   predictors.names <-  c( if( .nodrift ) NULL else 
       namesof(drift.n, link = .ldrift , earg = .edrift, tag = FALSE),
          namesof(varsd.n, link = if ( .var.arg ) .lvar  else .lsd ,
                           earg = if ( .var.arg ) .evar  else .esd ),
    if (!( .ARord2 )) NULL else
     namesof(ar.n, link = .lARcoeff , earg = .eARcoeff , tag = FALSE),
    if (!( .MAord2 )) NULL else
     namesof(ma.n, link = .lMAcoeff , earg = .eMAcoeff , tag = FALSE))
      
     predictors.names <- 
            predictors.names[interleave.VGAM(M1*NOS, M1 = M1)]
     
      if (!length(etastart)) {
        init.dr  <- matrix(0.0, nrow = n, ncol = NOS)
        init.sig <- matrix(0.0, nrow = n, ncol = NOS)
        init.AR  <- array( 0.0, dim = c( n, NOS , .ARord ))
        init.MA  <- array( 0.0, dim = c( n, NOS , .MAord ))
        
        
        y.sc <- scale(y, center = TRUE, scale = FALSE)  # y 180519
        extra$res <- matrix(NA_real_, nrow = n, ncol = NOS)
        res2.mat  <- matrix(NA_real_, nrow = n, ncol = NOS)
        
        for (rsp in 1:NOS) {
          to.fit <- cbind(y.sc[, rsp, drop = FALSE],
                          WN.lags(y = y.sc[, rsp, drop = FALSE],
                          lags = .ARord + .MAord + 4))
          
          
          to.fit2 <- lm( to.fit[, 1, drop = FALSE] ~
                           to.fit[,  -1, drop = FALSE] - 1)
       
          
          extra$res <- res2.mat[, rsp] <-  cbind(residuals(to.fit2))
          
          if (TRUE) { 
            to.fit2 <- cbind(y.sc[, rsp, drop = FALSE],
                       if ( .ARord2 )
                         WN.lags(y = cbind(y.sc[, rsp, drop = FALSE]),
                          lags = .ARord ) else NULL,
                           if ( .MAord2 )
                          WN.lags(y = res2.mat[, rsp, drop = FALSE],
                                  lags = .MAord ) else NULL)
            
          to.fit2 <- lm( to.fit2[, 1, drop = FALSE] ~
                         to.fit2[,  -1, drop = FALSE] - 1)
     }
   
    initCoe <- initCoebis <- coef(to.fit2)
     extra$res <- res2.mat[, rsp] <-  cbind(residuals(to.fit2))
          
          
      if ( .ARord2 != 0 ) {
        
        just.hlp <-matrix(if (length( .iARcoeff )) iniar[ rsp, ] else
          initCoe[ 1:( .ARord )] , n, .ARord , byrow = TRUE)
        init.AR[, rsp, ] <- just.hlp
        initCoe <- initCoe[ - ( 1:( .ARord ) ) ] 
      } 
       
       
       if ( .MAord2 != 0 ) {
        
        just.hlp <-matrix(if (length( .iMAcoeff )) inima[ rsp, ] else
          initCoe, n, .MAord , byrow = TRUE)
        init.MA[, rsp, ] <- just.hlp
        
        
      }
       
       print(head(init.AR))
       print(head(init.MA))
       
       #stop("")
       
       #print( .ARord2 )     
       #print(  .ARord )
       
       #print( .MAord2 )     
       #print(  .MAord )
          
      #if ( .MAord2 != 0 ) {
        
      #  just.hlp <-matrix(if (length( .iMAcoeff )) inima[ rsp, ] else
      #    initCoe[-(1:( .ARord ))], n, .MAord , byrow = TRUE)
      #  init.MA[, rsp, ] <- just.hlp
        
      #}
        
       init.dr[, rsp] <- ifelse( length( .idrift ) ,
                                 myIdrift[rsp], 
                                 ifelse( .ARord2 == 0 , 
                                         mean(y[, rsp]),
              mean(y[, rsp]) * (1 - sum(initCoebis[1:( .ARord )]))                            ))
              
          
        init.sig[, rsp] <- if( .var.arg ) var(res2.mat[, rsp]) else
                                                sd(res2.mat[, rsp])

  }  # End of for (1:NOS)
        
      
         etastart  <- cbind(if ( .nodrift ) NULL else init.dr ,
                           init.sig, 
                           matrix(init.AR, n, NOS * ( .ARord )), 
                           matrix(init.MA, n, NOS * ( .MAord )))
        
        etastart <- cbind(if ( .nodrift ) NULL else
                theta2eta(init.dr, .ldrift , .edrift ),
               theta2eta(init.sig, if ( .var.arg ) .lvar else .lsd ,
                              if ( .var.arg ) .evar else .esd ),
              if (!( .ARord2 )) NULL else
                    theta2eta(matrix(init.AR, n, NOS * ( .ARord )),
                             .lARcoeff , .eARcoeff ), 
              if (!( .MAord2 )) NULL else
                    theta2eta(matrix(init.MA, n, NOS * ( .MAord ) ),
                              .lMAcoeff , .eMAcoeff))
        etastart <- etastart[, interleave.VGAM(M1 * NOS, M1 = M1)]
        
        etastart
        
      } # End of !etastart
    }), list( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
              .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .evar = evar , .esd = esd ,
              .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff , 
              .idrift = idrift , .ivar = ivar , .isd = isd ,
              .iARcoeff = iARcoeff , .iMAcoeff = iMAcoeff ,
              .var.arg = var.arg , .realM1 = realM1 ,
              .ARord = ARord , .MAord = MAord , .nOrder = nOrder ,
              .nodrift = nodrift , .ARord2 = ARord2 ,
              .MAord2 = MAord2 , .dOrd = dOrd ))),
   
   
   
   
   
    linkinv = eval(substitute(function(eta, extra = NULL){
      
      M1  <- .realM1; n <- nrow(eta)
      NOS <- ncol(eta)/M1
      mat.save <- extra$mat.save
      dOrd < .dOrd; NNro <- 1
      nn <- NROW(mat.save)
      
      
      drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
             eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                          link =  .ldrift , earg = .edrift )
      
      sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder, drop = FALSE],
                         link = if ( .var.arg ) .lvar else .lsd ,
                         earg = if ( .var.arg ) .evar else .esd )
      
      ars <- array(0, dim = c(n, NOS, .ARord ))
      mas <- array(0, dim = c(n, NOS, .MAord ))
      
      if ( .ARord2 )
      for (ai in 1:( .ARord )) 
        ars[, , ai] <- eta2theta(eta[, M1 * (1:NOS) -
                            ( .ARord + .MAord ) + ai, drop = FALSE], 
                                 link = .lARcoeff , earg = .eARcoeff )
      if ( .MAord2 )
      for ( mi in 1:( .MAord )) 
        mas[, , mi] <- eta2theta(eta[, M1 * (1:NOS) -
                           ( .MAord ) + mi , drop = FALSE],
                                 link = .lMAcoeff , earg = .eMAcoeff )
      y.est <- matrix(NA_real_, n , NOS)

      for (ii in 1:NOS) 
         y.est[, ii] <- y.tt <- drifts[, ii] +  
          rowSums(cbind(ars[, ii, ] * 
                  WN.lags(y = cbind(extra$y[, ii]), lags = .ARord ),
                    mas[, ii, ] * WN.lags(y = cbind(extra$res[, ii]), 
                                lags = .MAord )))# + extra$res[, ii] 
      
      extra$res <- extra$y - y.est
      sum.mat <- cbind(extra$y.r, 
                          extra$mat.save[, -c(dOrd), drop = FALSE],
                       c(y.est, rep(0, dOrd )))
      
      t.rev <- rev(1:NCOL(sum.mat))
      for (ii in 1:dOrd) {
        yaux1 <- sum.mat[c(1:(nn - dOrd + ii - 1)), 
                         t.rev[ii], drop = FALSE]
        yaux2 <- sum.mat[c(1:(nn - dOrd + ii - 1)),
                         t.rev[ii] - 1, drop = FALSE]
        yaux3 <- yaux1 + yaux2
        yaux3 <- rbind(sum.mat[1, t.rev[ii] - 1, drop = FALSE], yaux3)
        sum.mat[ , t.rev[ii] - 1] <- c(yaux3, rep(0, dOrd - ii))
      }
      
      y.est <- sum.mat[, 1]
      y.est[c(1:NNro)] <- (NNro + 0.10)*y.est[c(2:(NNro + 1))]
      y.est
      
    }, list( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
             .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
             .edrift = edrift , .evar = evar , .esd = esd ,
             .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff , 
             .idrift = idrift , .ivar = ivar , .isd = isd ,
             .iARcoeff = iARcoeff , .iMAcoeff = iMAcoeff ,
             .var.arg = var.arg , .realM1 = realM1 ,
             .nodrift = nodrift , .nOrder = nOrder ,
             .ARord  = ARord , .MAord = MAord , .dOrd = dOrd , 
             .ARord2 = ARord2 , .MAord2 = MAord2 ))),
   
   
   
   
   
    last = eval(substitute(expression({
      
      M1 <- .realM1
      drifts <- if ( !(.nodrift) ) 
         eta2theta(eta[, M1*(1:NOS) - .nOrder - 1 , drop = FALSE], 
                  link = .ldrift , earg = .edrift ) else 
                    if (.nodrift ) matrix(0.0, nrow = n, ncol = NOS)
      
      sigWN <- if ( .var.arg )
        eta2theta(eta[, M1*(1:NOS) - .nOrder  , drop = FALSE], 
                  link = .lvar , earg = .evar ) else 
          eta2theta(eta[ , M1*(1:NOS) - .nOrder , drop = FALSE], 
                    link = .lsd , earg = .esd )
      
      ars <- mas <- NULL
      prov.names <- c(if ( .nodrift ) NULL else drift.n , varsd.n)
      if ( .ARord2 ) {
        ars <- array(NA, dim = c(n, NOS, .ARord ))
        for (ai in 1:( .ARord ))
          ars[, , ai ] <- 
            eta2theta(eta[, M1*(1:NOS) - (.MAord + .ARord) + ai, 
                          drop = FALSE], 
                      link = .lARcoeff , earg = .eARcoeff )
        for (kk in 1:( .ARord ))
          prov.names <- c(prov.names, ar.n[[kk]])
      }
      
      if ( .MAord2 ) {
        mas <- array(NA, dim = c(n, NOS, .MAord ))
        for ( mi in 1:( .MAord ))
          mas[, , mi ] <-
            eta2theta(eta[, M1*(1:NOS) - ( .MAord ) + mi , 
                           drop = FALSE],
                      link = .lMAcoeff , earg = .eMAcoeff )
        for (kk in 1:( .MAord ))
          prov.names <- c(prov.names, ma.n[[kk]])
      }
      
      misc$link <- rep( .ldrift , times = M1 * NOS )
      names(misc$link) <- 
        prov.names[interleave.VGAM(M1 * NOS , M1 = M1)]
      misc$earg <- vector("list", length = M1 * NOS)
      names(misc$earg) <- 
        prov.names[interleave.VGAM(M1 * NOS , M1 = M1)]
      
  for (jj in 1:NOS) {
        
    if (!( .nodrift ))
      misc$link[ M1 * jj - .nOrder - 1 ] <- .ldrift
      misc$link[M1 * jj - .nOrder] <- if ( .var.arg ) .lvar else .lsd 
        
    if (!( .nodrift ))
      misc$earg[[ M1 * jj - .nOrder - 1 ]] <- .edrift 
    misc$earg[[M1 * jj - .nOrder]] <- if ( .var.arg ) .evar else .esd
        
        if ( .ARord2 )
        for (kk in 1:( .ARord )) {
          misc$link[  M1 * jj - .nOrder + kk ]  <- .lARcoeff 
          misc$earg[[ M1 * jj - .nOrder + kk ]] <- .eARcoeff 
        }
        
        if ( .MAord2 )
        for( ll in 1:( .MAord )) {
          misc$link[  M1 * jj - .MAord + ll ]  <- .lMAcoeff 
          misc$earg[[ M1 * jj - .MAord + ll ]] <- .eMAcoeff
        }  
      }
      
      misc$noChecks  <- .noChecks
      misc$expected  <- TRUE
      misc$process   <- "ARMA"
      misc$var.arg   <- .var.arg 
      misc$nomean    <- .nodrift
      misc$Order     <- c( .ARord , .MAord )
      misc$M1        <- .realM1
      misc$NOS       <- NOS
      misc$ARord     <- .ARord
      misc$MAord     <- .MAord
      misc$flagAR    <- .flagAR
      misc$residuals <- cbind(c(rep(0, .dOrd ), extra$res))
      misc$theta.names <- parameters.names
      misc$multipleResponses <- TRUE
      misc$control$wzepsilon <- control$wzepsilon
      fit$prior.weights <- w #20180315
      y <- c(extra$y.r) 
      
      # Checks removed on March/2/2016. Now, refer to 'summary'.
      if (!( .noChecks ) && !( .flag.2 )) {
        
        M <- M1 * NOS
        FlagArma <- grep("(Intercept)", names(fit$coefficients))
        
        if (length(FlagArma) == M) {
          coefArma <- fit$coefficients
          coefMa   <- coefAr  <- numeric(0)
          for ( kk in 1:NOS ) {
            myAux2 <- coefArma[ M1 * kk - .nOrder + 1:( .ARord )]
            myAux3 <- coefArma[ M1 * kk - .MAord + 1:( .MAord )]
            coefAr <- c(coefAr, myAux2)
            coefMa <- c(coefMa, myAux3)
          }
      coefAr <- eta2theta(coefAr, link = .lARcoeff , earg = .eARcoeff)
      coefMa <- eta2theta(coefMa, link = .lMAcoeff , earg = .eMAcoeff)
          
          # 2016/March/01. Roots printed out by 'summary(vglm)'. '.
          checkAR <- checkTS.ffs(thetaEst = coefAr,
                                 tsclass  = "AR",
                                 NofS     = NOS, 
                                 chOrder  = .ARord ,
                                 pRoots   = FALSE, 
                                 retmod   = TRUE)
          
          checkMA <- checkTS.ffs(thetaEst = coefMa, 
                                 tsclass  = "MA",
                                 NofS     = NOS, 
                                 chOrder  = .MAord ,
                                 pRoots   = FALSE,
                                 retmod   = TRUE)
          
          flag1 <- (all(checkAR > 0) && any(checkAR < 1 + 5e-3))
          flag1 <- 
              ((all(checkMA > 0) && any(checkMA < 1 + 5e-3)) || flag1)
        
        if (!flag1)
        cat("\nChecks on stationarity / invertibility successfully",
            "performed. \nNo roots liying inside the unit circle.",
             "\nFurther details within the 'summary' output.\n") else
              
        cat("\nChecks on stationarity / invertibility successfully",
                "performed. \nNOTE: Some roots lie inside the unit ",
                  "circle. \nFurther details within the 'summary' ", 
                  "output.\n")
        } else {
          cat("\nApparently, some constraints have been set.",
              "\nDetails on stationarity / invertibility",
              "checks \nwhitin the 'summary' output.")
        } 
      } # End of !( .noChecks )
      
    }), list( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
              .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .evar = evar , .esd = esd ,
              .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff , 
              .var.arg = var.arg , .realM1 = realM1 ,
              .ARord = ARord , .MAord = MAord , 
              .nodrift = nodrift , .noChecks = noChecks ,
              .nOrder = nOrder, .flagAR = flagAR ,
              .type.likelihood = type.likelihood ,
              .ARord2 = ARord2 , .MAord2 = MAord2,
              .flag.2 = flag.2 , .dOrd = dOrd ))),
   
   
   
   
   
    loglikelihood = 
        eval(substitute(function(mu, y, w, residuals = FALSE,
                                 eta, extra = NULL, summation = TRUE) {
          
      y   <- extra$y
      M1  <- .realM1; n <- nrow(eta)
      NOS <- ncol(eta)/M1
      eim.ap <- ( .type.EIM == "approximate")
      
      drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
        eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                  link =  .ldrift , earg = .edrift )
      y.est <- drifts + log(1) * ifelse(M1 > 1, 0.15, 0.15) * extra$res
      ars <- mas <- 0
      
      print(head(drifts))
      
      if ( .var.arg ) {
        sigWN <- sqrt(eta2theta(eta[, M1*(1:NOS) - .nOrder,
                                    drop = FALSE], .lvar , .evar ))
      } else {
        sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder,
                               drop = FALSE], .lsd , .esd )
      }
      
      if ( .ARord2 ) {
        #ars <- matrix(NA_real_, n, .ARord )
        ars <- eta2theta(eta[, (M1  - ( .ARord2 + .MAord2 ) 
                                + 1):(M1 - .MAord2 ), drop = FALSE], 
                         link = .lARcoeff , earg = .eARcoeff )
        
     y.est <- y.est + cbind(rowSums(ars * WN.lags(y = cbind(extra$y), 
                                                   lags = .ARord2 )))
        
      } 
      
      if ( .MAord2 ) {
        #ars <- matrix(NA_real_, n, .ARord )
        mas <- eta2theta(eta[, (M1  - ( .MAord2 ) + 1):M1,
                             drop = FALSE], 
                         link = .lMAcoeff , earg = .eMAcoeff )
        
    y.est <- y.est + cbind(rowSums(mas * WN.lags(y = cbind(extra$res), 
                                                   lags = .MAord2 )))
      } 
      
      
       if (residuals) {
         stop("Loglikelihood not implemented yet to",
              " handle residuals.") 
       } else {
         loglik.terms <- c(w) * dARMA(x = y, mean = y.est ,
                                      sd = sigWN, log = TRUE)
       }
      
       if ( summation ) 
         loglik.terms <- sum(loglik.terms) 
        #if (!eim.ap) sum(loglik.terms) else
                                        #round(sum(loglik.terms), 4)
      
      loglik.terms
        
    }, list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar ,
             .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
             .edrift = edrift , .esd = esd , .evar = evar ,
             .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
             .nOrder = nOrder , .ARord = ARord , .MAord = MAord ,
             .type.likelihood = type.likelihood ,
             .nodrift = nodrift , .type.EIM = type.EIM ,
             .var.arg = var.arg , .realM1 = realM1 ,
             .ARord2 = ARord2 , .MAord2 = MAord2 ))),
    
   
   
   
   
    validparams = eval(substitute(function(eta, y, extra = NULL) {
      
     M1  <- .realM1; n <- nrow(eta)
     NOS <- ncol(eta)/M1
     
     drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else
       eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                 link =  .ldrift , earg = .edrift )
     
     sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder, drop = FALSE],
                        link = if ( .var.arg ) .lvar else .lsd ,
                        earg = if ( .var.arg ) .evar else .esd )
     
     ars <- array(0, dim = c(n, NOS, .ARord ))
     mas <- array(0, dim = c(n, NOS, .MAord ))
     
     for (ai in 1:( .ARord )) 
       ars[, , ai] <- eta2theta(eta[, M1 * (1:NOS) -
                            ( .ARord + .MAord ) + ai, drop = FALSE],
                                link = .lARcoeff , earg = .eARcoeff )
     
     for ( mi in 1:( .MAord )) 
       mas[, , mi] <- eta2theta(eta[, M1 * (1:NOS) -
                                   ( .MAord ) + mi , drop = FALSE],
                                link = .lMAcoeff , earg = .eMAcoeff )
     
     okay1 <- (all(is.finite(sigWN)) && all(0 < sigWN) &&
                 all(is.finite(ars)) && all(is.finite(mas)) &&
                 all(is.finite(drifts)))
     okay1
     
   }, list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar ,
            .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
            .edrift = edrift , .esd = esd , .evar = evar ,
            .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
            .nOrder = nOrder , .ARord = ARord , .MAord = MAord ,
            .type.likelihood = type.likelihood ,
            .nodrift = nodrift ,
            .var.arg = var.arg , .realM1 = realM1 ))),
   
   
   
    #vfamily = c("vgtsff", "vgltsff-class", "ARIMAff"),
    vfamily  = c("ARIMAXff", "ARMAvgltsmff"),
   
   
    
   
    deriv = eval(substitute(expression({
      M1  <- .realM1; n <- nrow(eta)
      NOS <- ncol(eta)/M1
      M   <- M1 * NOS
      y   <- extra$y
      
      drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
        eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                  link =  .ldrift , earg = .edrift )
      y.est <- drifts
      ars <- mas <- 0
      
      if ( .var.arg ) {
        sigWN <- sqrt(eta2theta(eta[, M1*(1:NOS) - .nOrder,
                                    drop = FALSE], .lvar , .evar ))
      } else {
        sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder,
                               drop = FALSE], .lsd , .esd )
      }
      
      if ( .ARord2 !=0 ) {
        #ars <- matrix(NA_real_, n, .ARord )
        ars <- eta2theta(eta[, (M1  - ( .ARord2 + .MAord2 ) 
                                + 1):(M1 - .MAord2 ), drop = FALSE], 
                         link = .lARcoeff , earg = .eARcoeff )
        
        y.est <- y.est + 
                    cbind(rowSums(ars * WN.lags(y = cbind(extra$y), 
                                                 lags = .ARord2 )))
        
      } 
      
      if ( .MAord2 !=0 ) {
        #ars <- matrix(NA_real_, n, .ARord )
        mas <- eta2theta(eta[, (M1  - ( .MAord2 ) + 1):M1,
                             drop = FALSE], 
                         link = .lMAcoeff , earg = .eMAcoeff )
        
        y.est <- y.est + 
           cbind(rowSums(mas * WN.lags(y = cbind(extra$res), 
                                                  lags = .MAord2 )))
    
      } 
      
      dl.drfmean <- (y - y.est) / sigWN^2
      if ( .var.arg ) {
        dl.dvar <- (y - y.est)^2 / (2 * sigWN^4) - 1 / (2 * sigWN^2)
      } else {
        dl.dsd  <- (y - y.est)^2 / sigWN^3 - 1 / sigWN
      }
      
      dl.dThe <- array(NA_real_, dim = c(n, NOS, .ARord ))
      dl.dPhi <- array(NA_real_, dim = c(n, NOS, .MAord ))
      
   for (ii in 1:NOS) 
    dl.dThe[, ii, ] <- WN.lags(y = cbind(y[, ii]), lags = .ARord ) *
          matrix((y[, ii]- y.est[, ii]), n , .ARord ) / sigWN[, ii]^2
     
   for (ii in 1:NOS) 
      dl.dPhi[, ii, ] <- WN.lags(y = cbind(extra$res[, ii]),
                                   lags = .MAord ) *
         matrix((y[, ii] - y.est[, ii]), n , .MAord ) / sigWN[, ii]^2
      
      ddrif.deta <- dtheta.deta(drifts , .ldrift , .edrift)
      dsewn.deta <- if (.var.arg) dtheta.deta(sigWN^2, 
                                              .lvar , .evar) else
                                      dtheta.deta(sigWN, .lsd , .esd)
      
      dThe.deta  <- dtheta.deta(ars , .lARcoeff , .eARcoeff)
      dPhi.deta  <- dtheta.deta(mas , .lMAcoeff , .eMAcoeff)
      
      myderiv <- c(w) * cbind(if ( .nodrift ) NULL else
                                         dl.drfmean * ddrif.deta,
                            if ( .var.arg ) dl.dvar * dsewn.deta else
                                               dl.dsd * dsewn.deta,
        if (!( .ARord2 )) NULL else
          matrix( dl.dThe[, 1, ] * dThe.deta, n, NOS * ( .ARord )),
        if (!( .MAord2 )) NULL else
            matrix( dl.dPhi[, 1, ] * dPhi.deta, n, NOS * ( .MAord )))
      
      myderiv <- 
        myderiv[, interleave.VGAM(( .realM1 )*NOS , M1 = .realM1)]
      
      colnames(myderiv) <- NULL
      
      # Required at @weights
      ars <- array(ars, dim = c(n, NOS, .ARord2 ))
      mas <- array(mas, dim = c(n, NOS, .MAord2 ))
      dThe.deta  <- dtheta.deta(ars , .lARcoeff , .eARcoeff)
      dPhi.deta  <- dtheta.deta(mas , .lMAcoeff , .eMAcoeff)
      
      myderiv
      
    }), list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
              .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .esd = esd , .evar = evar ,
              .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
              .ARord = ARord , .MAord = MAord , .nOrder = nOrder ,
              .var.arg = var.arg , .realM1 = realM1,
              .nodrift = nodrift , .ARord2 = ARord2 ,
              .MAord2 = MAord2 ))),
   
   
   
   
   
    weight = eval(substitute(expression({
      
      if (( .type.EIM ) == "exact" ) {
        
        final.eim <- numeric(0)
        for (ii in 1:NOS) {
          if ( ( .ARord2 ) && !( .MAord2 ) )
              exact.eim <-  ARpEIM.G2(y  = cbind(y[, ii]),
                                 drift   = NULL,
                                 sdError = cbind(sigWN[, ii]),
                                 ARpcoeff = cbind(ars[, ii, ]),
                                 var.arg = .var.arg ,
                                 order   = .ARord2 ,
                                 nodrift = .nodrift )
          if ( !( .ARord2 ) && ( .MAord2 ) )
            exact.eim <-  MAqEIM.G2(y  = cbind(y[, ii]),
                              mean  = NULL,
                              sdError  = cbind(sigWN[, ii]),
                              MAqcoeff = cbind(mas[, ii, ]),
                              var.arg  = .var.arg ,
                              order    = .MAord2 ,
                              nomean   = .nodrift )
          
          if ( ( .ARord2 ) && ( .MAord2 ) )
          exact.eim <- ARMA.EIM.G2(y = cbind(y[, ii]),
                              arCoe   = cbind(ars[, ii, ]),
                              maCoe   = cbind(mas[, ii, ]),
                              estRes  = cbind(extra$res[, ii]),
                              sdError = cbind(sigWN[, ii]),
                              var.arg = .var.arg ,
                              order   = c(.ARord , .MAord ),
                              nodrift = .nodrift )
          comb.wz <- combVGAMextra(1:M1)
          a.bind  <- cbind(if ( .nodrift ) NULL else ddrif.deta[, ii],
                           dsewn.deta[, ii],
                           dThe.deta[, ii, ],
                           dPhi.deta[, ii, ])
          
          dthdeta <- apply(comb.wz, 1, function(x) {
            a.bind[, x[1]] * a.bind[, x[2]]
          })
          
          final.eim <- cbind(final.eim, 
                             c(w[, ii]) * exact.eim * dthdeta)
        }
        
         final.eim <-
                 final.eim[, interleave.VGAM( NOS * M1 * (M1 + 1) / 2,  
                                                       M1 = NOS)]
         final.eim <- array(c(final.eim), 
                             dim = c(n, NOS, M1 * (M1 + 1)/ 2))
         wz <- arwz2wz(final.eim, M = M1 * NOS, M1 = M1)
        
      } else {
        
        ned2l.dsmn   <- (1 / sigWN^2) 
        ned2l.dvarSD <- 
          if ( .var.arg ) 1 / (2 * sigWN^4) else 2 / sigWN^2
        
        gammas.y <- apply(y, 2, function(x) {
          cross.gammas(x = x, lags = 0)
        })
        
        ned2l.dthe <- 
           matrix(gammas.y, n, NOS * ( .ARord ), byrow = TRUE) /
                              matrix(sigWN^2, n, NOS * ( .ARord ))
        
        dThe.deta <- matrix(dThe.deta, n, NOS * ( .ARord ))
        
        ned2l.dphi <- matrix(1, n, NOS * ( .MAord ))
        dPhi.deta  <- matrix(dPhi.deta, n , NOS * ( .MAord ))
        
        wz <- c(w) * cbind( if ( .nodrift ) NULL else 
             ned2l.dsmn * ddrif.deta^2, ned2l.dvarSD * dsewn.deta^2,
             if (!( .ARord2 )) NULL else ned2l.dthe * dThe.deta^2 ,
             if (!( .MAord2 )) NULL else ned2l.dphi * dPhi.deta^2)
        
        wz <- wz[, interleave.VGAM( NOS * M1 , M1 = M1)]
    }
      
      wz[which(wz[, 1:(NOS * M1)] < .Machine$double.eps)] <- 
                                           .Machine$double.eps^(0.75)
      
      wz
      
    }), list( .var.arg = var.arg , .nOrder   = nOrder, 
              .ARord   = ARord   , .MAord   = MAord ,
              .realM1  = realM1  , .nodrift = nodrift ,
              .idrift = idrift ,
              .iARcoeff = iARcoeff , .iMAcoeff = iMAcoeff ,
              .ivar = ivar, .isd = isd , .type.EIM = type.EIM ,
              .ARord2 = ARord2 , .MAord2 = MAord2 )))
    ) # End of alo
  
   alo
  
  } 
  
