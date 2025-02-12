#####################################################################
# These functions are
# Copyright (C) 2014-2024 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0
# Supports the Wald, score, and lrt tests (20180209)


ARMA.EIM.G2 <- function(y, arCoe, maCoe,  estRes, sdError,
                         var.arg = TRUE, order = c(1, 1),
                         nodrift = FALSE, addRidge = 0.0001) {
  
  if (length(order) !=2 )
    stop("Currently, this function only handles ARMA (p, q) models.")
  
  Order <- order; rm(order)
  arOrd <- Order[1]; maOrd <- Order[2]
  y <- cbind(y); nn <- nrow(y)
  arCoe <- matrix(arCoe, nrow = nn, ncol = arOrd)
  
  # R matrix
  RMat  <- diag(c(sdError^2))
  dRMat <- cbind( if(nodrift) NULL else rep(0, nn),
                  if (var.arg)  1 else 2 * sdError,
                  matrix(0, nn, arOrd + maOrd))
  
  # dm(t) matrix
  arLags <- WN.lags(y = cbind(y), lags = arOrd,
                    to.complete = 1/(1 - arCoe[1, ] + 0.01)^2)
  maLags <- WN.lags(y = cbind(estRes), lags = maOrd,
                    to.complete = 1/(1 - maCoe[1, ] + 0.01)^2)
  dmt <- cbind( if (nodrift) NULL else rep(1, nn),
                rep(0, nn),
                WN.lags(y = cbind(y), lags = arOrd,
                        to.complete = 1/(1 - arCoe[1, ] + 0.01)^2),
                WN.lags(y = cbind(estRes), lags = maOrd,
                        to.complete = 1/(1 - maCoe[1, ] + 0.01)^2))
  
  
  M <- 2 + arOrd + maOrd # ARMA model
  try.comb <- iam(NA, NA, M = M - nodrift, both = TRUE)
  try.comb <- cbind(try.comb$row.index, try.comb$col.index)
  
  ### This is equation A.6, Porat and Friedlander (1986)
  finMat <- apply(try.comb, 1, function(x){
    kl <- x
    invRMat <- diag(solve(RMat))
    dRthel  <- c(dRMat[, kl[1]])
    dRthek  <- c(dRMat[, kl[2]])
    term.1  <- (0.5) * invRMat * dRthel * invRMat * dRthek
    term.2  <- c(dmt[, kl[1]]) * invRMat * c(dmt[, kl[2]])
    
    term.1 + term.2
  })
  
  finMat[, 1:M] <- (1 + addRidge) * finMat[, 1:M]
  
  finMat
  
}




# ARMApq  density
dARMA <- function(x, mean = 0, sd = 1, log = FALSE) 
  dnorm(x = x, mean = mean, sd = sd, log = log)





ARMAXff <- 
  function(order     = c(1, 1),
           zero      = c(if (nodrift) NULL else "drift",
                         "ARcoeff", "MAcoeff"),
           xLag      = 0,
           type.EIM  = c("exact", "approximate")[1],
           var.arg   = TRUE,
           nodrift   = FALSE,
           noChecks  = FALSE,
           ldrift    = "identitylink",
           lsd       = "loglink",
           lvar      = "loglink",
           lARcoeff  = "identitylink",
           lMAcoeff  = "identitylink",
           idrift    = NULL,      # Must be a vector of length NOS
           isd       = NULL,      # Must be a vector of length NOS
           ivar      = NULL,      # Must be a vector of length NOS
           iARcoeff  = NULL,      # Must be a vectors
           iMAcoeff  = NULL) {
    
    
    if ( !Is.Numeric(order, length.arg = 2, isInteger = TRUE )  )
      stop("Invalid 'order'. Should be a vector (p, q) of integers.")
    
    if (xLag < 0)
      stop("Wrong input for argument 'xLag'")
    
    ARord  <- order[1]
    MAord  <- order[2]; rm(order)
    nOrder <- ARord + MAord
    realM1 <- 2 + ARord + MAord - nodrift
    
    if ( ARord == 0 || MAord == 0 )
      stop("Only ARMAX(p, q) models, p, q > 0, handled." ,
           " Please, refer to ARXff() or MAXff().")
    
    if (!is.logical(noChecks) && !is.logical (var.arg) &&
          !is.logical(nodrift))
      stop("'noChecks', 'nodrift', and 'var.arg' must be logical")
    
    type.likelihood <- "exact"
    type.EIM <- match.arg(type.EIM, c("exact", "approximate"))[1]
    
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
      new( "vgltsmff",
          blurb = c("VGLTMs: ", nOrder + 1 + !nodrift,
                    "-parameter ARMAX model of order",
                    "(", ARord, ",", MAord, ").",
                    "\n", "\n",
                    "LINKS: ",
                    blurb.vec2, "\n",
                    "Model: Y[t] = ", 
                    if (nodrift) NULL else "drift + ",
                    " B^T * X_t + ",
                    if (ARord == 1) preAR else
               c("(ARcoeff1) * Y[t - 1] +... + (ARcoeff", ARord,")",
                        " * Y[t -", ARord,"] + "),
                    if (MAord == 1) preMA else
               c(" (MAcoeff1) * e[t - 1] + ... + (MAcoeff", MAord,")",
                        " * e[t - ", MAord,"]") , "+ e[t].", "\n"),
                   #"E[Yt] = drift / (1 - ARcoeff1 - ... - ARcoeff",
                   #ARord, ")"),
          
          
          
          
          
   constraints = eval(substitute(expression({
      
      M1      <- .realM1
      constraints <- cm.zero.VGAM(constraints, x = x, 
                                  zero = .zero , M = M,
                                  # Old : predictors.names...
                                  predictors.names = parameters.names,
                                  M1 = M1)
    }), 
    list( .zero = zero , .realM1 = realM1 ) )),
   
   
   
   first = eval(substitute(expression({
     
     if ((NCOL(x) == 1) && 
              (colnames(x) == "(Intercept)") && ( .xLag > 0))
       stop("Set  xLag = 0. No covariates entered.")
     
     if ( .xLag ) {
       x.matrix <- x
       for (jj in 1:(NCOL(x) - 1)) {
         temp1 <-  WN.lags(cbind(x[, 1 + jj]), .xLag )
         temp2 <- paste(colnames(x)[1 + jj], "Lag", sep = "")
         colnames(temp1) <- paste(temp2, 1:( .xLag), sep = "")
         x.matrix<- cbind(x.matrix, temp1)
       }
       rm(temp1, temp2)
       list.names <- vector("list", NCOL(x.matrix))
       names(list.names) <- colnames(x.matrix)
       for (ii in 1:(NCOL(x.matrix))) 
         list.names[[ii]] <- ii
       attr(x.matrix, "assign") <- list.names
       
       x <- x.matrix
     }
     
   }), list( .xLag = xLag ))),
   
   
   
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
               .ARord = ARord , .MAord = MAord , 
               .nodrift = nodrift , .type.EIM = type.EIM ,
               .zero = zero ))),
   
   
   
   
   
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
      M   <- NOS * M1; n <- nrow(y)
      nOrder <- .nOrder
      
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
      
      parameters.names <- c(if ( .nodrift ) NULL else drift.n,
                            varsd.n, ar.n, ma.n)
   
      predictors.names <-  c( if( .nodrift ) NULL else 
      namesof(drift.n, link = .ldrift , earg = .edrift, tag = FALSE),
      namesof(varsd.n, link = if ( .var.arg ) .lvar  else .lsd ,
                           earg = if ( .var.arg ) .evar  else .esd ),
     namesof(ar.n, link = .lARcoeff , earg = .eARcoeff , tag = FALSE),
     namesof(ma.n, link = .lMAcoeff , earg = .eMAcoeff , tag = FALSE))
      
     predictors.names <- 
               predictors.names[interleave.VGAM(M1*NOS, M1 = M1)]
     
      if (!length(etastart)) {
        init.dr  <- matrix(0.0, nrow = n, ncol = NOS)
        init.sig <- matrix(0.0, nrow = n, ncol = NOS)
        init.AR  <- array( 0.0, dim = c( n, NOS , .ARord ))
        init.MA  <- array( 0.0, dim = c( n, NOS , .MAord ))
        
        y.sc <- y# scale(y, center = TRUE, scale = FALSE)
        #extra$y.sc <- scale(y, center = TRUE, scale = FALSE)
        extra$res <- matrix(NA_real_, nrow = n, ncol = NOS)
        res2.mat  <- matrix(NA_real_, nrow = n, ncol = NOS)
        
        for (rsp in 1:NOS) {
          to.fit <- cbind(y.sc[, rsp, drop = FALSE],
                         WN.lags(y = cbind(y.sc[, rsp, drop = FALSE]),
                         lags = .ARord + .MAord + 4))
          
          to.fit <- lsfit(x = to.fit[, -1, drop = FALSE],
                          y = to.fit[,  1, drop = FALSE],
                          intercept = FALSE)
          
          extra$res <- res2.mat[, rsp] <- residuals(to.fit)
          # Updating residuals not allowed.
          #res2.mat[, rsp] <- residuals(to.fit)
          
          # Initial values.
          if (TRUE) {
          to.fit <- cbind(y.sc[, rsp, drop = FALSE],
                        WN.lags(y = cbind(y.sc[, rsp, drop = FALSE]),
                          lags = .ARord ),
                        WN.lags(y = res2.mat[, rsp, drop = FALSE],
                                  lags = .MAord ))
          to.fit <- lsfit(x = to.fit[, -1, drop = FALSE],
                          y = to.fit[,  1, drop = FALSE],
                          intercept = FALSE)
          }
          
          initCoe <- coef(to.fit)
          
      extra$res <- res2.mat[, rsp] <- residuals(to.fit)
      just.hlp <-matrix(if (length( .iARcoeff )) iniar[ rsp, ] else
                     initCoe[1:( .ARord )], n, .ARord , byrow = TRUE)
      init.AR[, rsp, ] <- just.hlp
          
      just.hlp <-matrix(if (length( .iMAcoeff )) inima[ rsp, ] else
                initCoe[-c(1:( .ARord ))], n, .MAord , byrow = TRUE)
      init.MA[, rsp, ] <- just.hlp
          
      init.dr[, rsp] <- if (length( .idrift )) myIdrift[rsp] else
                   mean(y[, rsp]) * (1 - sum(initCoe[1:( .ARord )]))
          
      init.sig[, rsp] <- if( .var.arg ) var(res2.mat[, rsp]) else
                                              sd(res2.mat[, rsp])
        }
        
  #res2.mat[c(1: max( .ARord, .MAord )), 1] <- median(res2.mat[, 1])
        extra$res <- res2.mat
        
        etastart  <- cbind(if ( .nodrift ) NULL else init.dr ,
                           init.sig, 
                           matrix(init.AR, n, NOS * ( .ARord )), 
                           matrix(init.MA, n, NOS * ( .MAord )))
        
        etastart <- cbind(if ( .nodrift ) NULL else
                theta2eta(init.dr, .ldrift , .edrift ),
                theta2eta(init.sig, if ( .var.arg ) .lvar else .lsd ,
                            if ( .var.arg ) .evar else .esd ),
                    theta2eta(matrix(init.AR, n, NOS * ( .ARord )),
                             .lARcoeff , .eARcoeff ), 
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
              .nodrift = nodrift ))),
   
   
   
   
   
    linkinv = eval(substitute(function(eta, extra = NULL){
      
      M1  <- .realM1; n <- nrow(eta)
      NOS <- ncol(eta)/M1; NNro <- 1
      
      drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
             eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                          link =  .ldrift , earg = .edrift )
      
      sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder, drop = FALSE],
                         link = if ( .var.arg ) .lvar else .lsd ,
                         earg = if ( .var.arg ) .evar else .esd )
      
      ars <- array(NA_real_, dim = c(n, NOS, .ARord ))
      mas <- array(NA_real_, dim = c(n, NOS, .MAord ))
      
      for (ai in 1:( .ARord )) 
        ars[, , ai] <- eta2theta(eta[, M1 * (1:NOS) -
                            ( .ARord + .MAord ) + ai, drop = FALSE], 
                                 link = .lARcoeff , earg = .eARcoeff )
        
      for ( mi in 1:( .MAord )) 
        mas[, , mi] <- eta2theta(eta[, M1 * (1:NOS) -
                           ( .MAord ) + mi , drop = FALSE],
                                 link = .lMAcoeff , earg = .eMAcoeff )
      
      y.est <- matrix(NA_real_, n , NOS)
      
      for (ii in 1:NOS) {
        y.est[, ii] <- drifts[, ii]  +  rowSums(cbind(ars[, ii, ] * 
                  WN.lags(y = cbind(extra$y[, ii]), lags = .ARord ),
                    mas[, ii, ] * WN.lags(y = cbind(extra$res[, ii]),
                            lags = .MAord ))) # + extra$res[, ii]
        y.est[c(1:NNro), ii] <- (NNro + 0.10) *
                               y.est[c(2:(NNro + 1)), ii]
      }  
      
      y.est
      
    }, list( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
             .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
             .edrift = edrift , .evar = evar , .esd = esd ,
             .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff , 
             .idrift = idrift , .ivar = ivar , .isd = isd ,
             .iARcoeff = iARcoeff , .iMAcoeff = iMAcoeff ,
             .var.arg = var.arg , .realM1 = realM1 ,
             .nodrift = nodrift , .nOrder = nOrder ,
             .ARord = ARord , .MAord = MAord ))),
   
   
   
   
   
    last = eval(substitute(expression({
      
      M1 <- .realM1
      res.mat <- extra$res
      drifts <- if ( !(.nodrift) ) 
         eta2theta(eta[, M1*(1:NOS) - .nOrder - 1 , drop = FALSE], 
                  link = .ldrift , earg = .edrift ) else 
                    if (.nodrift ) matrix(0.0, nrow = n, ncol = NOS)
      
      sigWN <- if ( .var.arg )
        eta2theta(eta[, M1*(1:NOS) - .nOrder  , drop = FALSE], 
                  link = .lvar , earg = .evar ) else 
          eta2theta(eta[ , M1*(1:NOS) - .nOrder , drop = FALSE], 
                    link = .lsd , earg = .esd )
      
      ars <- array(NA, dim = c(n, NOS, .ARord ))
      mas <- array(NA, dim = c(n, NOS, .MAord ))
      for (ai in 1:.ARord) 
        ars[, , ai ] <- 
        eta2theta(eta[, M1*(1:NOS) - (.MAord + .ARord) + ai, 
                      drop = FALSE], 
                  link = .lARcoeff , earg = .eARcoeff )
      for ( mi in 1:.MAord )
        mas[, , mi ] <-
        eta2theta(eta[, M1*(1:NOS) - ( .MAord ) + mi , drop = FALSE],
                  link = .lMAcoeff , earg = .eMAcoeff )
      
      prov.names <- c(if ( .nodrift ) NULL else drift.n , varsd.n)
      for (kk in 1:.ARord )
        prov.names <- c(prov.names, ar.n[[kk]])
      for (kk in 1:.MAord )
        prov.names <- c(prov.names, ma.n[[kk]])
      
      misc$link <- rep( .ldrift , times = M1 * NOS )
      names(misc$link) <- 
        prov.names[interleave.VGAM(M1 * NOS , M1 = M1)]
      misc$earg <- vector("list", length = M1 * NOS)
      names(misc$earg) <- 
        prov.names[interleave.VGAM(M1 * NOS , M1 = M1)]
      
      for (jj in 1:NOS) {
        
  if (!(.nodrift ))
    misc$link[ M1 * jj - .nOrder - 1 ] <- .ldrift
    misc$link[M1 * jj - .nOrder] <- if ( .var.arg ) .lvar else .lsd 
        
  if (!(.nodrift ))
    misc$earg[[ M1 * jj - .nOrder - 1 ]] <- .edrift 
    misc$earg[[M1 * jj - .nOrder]] <- if ( .var.arg ) .evar else .esd
        
        for (kk in 1:.ARord ) {
          misc$link[  M1 * jj - .nOrder + kk ]  <- .lARcoeff 
          misc$earg[[ M1 * jj - .nOrder + kk ]] <- .eARcoeff 
        }
        for( ll in 1:.MAord ) {
          misc$link[  M1 * jj - .MAord + ll ]  <- .lMAcoeff 
          misc$earg[[ M1 * jj - .MAord + ll ]] <- .eMAcoeff
        }  
      }
      
      extra$res[c(1: max( .ARord, .MAord )), 1] <- median(res.mat)
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
      misc$flagAR    <- FALSE
      misc$residuals <- extra$res
      misc$theta.names <- parameters.names
      misc$multipleResponses <- TRUE
      misc$control$wzepsilon <- control$wzepsilon
      fit$prior.weights <- w  # 20180315
      
      # Checks removed on March/2/2016. Now, refer to 'summary'.
      if (!( .noChecks )) {
        
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
    flag1 <- ((all(checkMA > 0) && any(checkMA < 1 + 5e-3)) || flag1)
    
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
              .nOrder = nOrder,
              .type.likelihood = type.likelihood ))),
   
   
   
   
   
    loglikelihood = 
      eval(substitute(function(mu, y, w, residuals = FALSE,
                               eta, extra = NULL, summation = TRUE) {
      
      #NOS <- extra$NOS
      y   <- extra$y
      M1  <- .realM1; n <- nrow(eta)
      NOS <- ncol(eta)/M1
      eim.ap <- ( .type.EIM == "approximate")
      
      drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
              eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                             link =  .ldrift , earg = .edrift )
     
      if ( .var.arg ) {
        sigWN <- sqrt(eta2theta(eta[, M1*(1:NOS) - .nOrder,
                     drop = FALSE], .lvar , .evar ))
      } else {
        sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder,
                               drop = FALSE], .lsd , .esd )
      }
     
      ars <- array(NA_real_, dim = c(n, NOS, .ARord ))
      mas <- array(NA_real_, dim = c(n, NOS, .MAord ))
     
      for (ai in 1:( .ARord ))
        ars[, , ai] <- eta2theta(eta[, M1 * (1:NOS) -
                           ( .ARord + .MAord ) + ai, drop = FALSE], 
                            link = .lARcoeff , earg = .eARcoeff )
     
      for ( mi in 1:( .MAord ))
        mas[, , mi] <- eta2theta(eta[, M1 * (1:NOS) -
                                 ( .MAord ) + mi , drop = FALSE],
                               link = .lMAcoeff , earg = .eMAcoeff )
     
      y.est <- matrix(NA_real_, n , NOS)
      print(head(cbind(ars)))
      print(head(cbind(mas)))
      
      
      
      for (ii in 1:NOS) 
        y.est[, ii] <- drifts[, ii] + rowSums(cbind(ars[, ii, ] *
                    WN.lags(y = cbind(extra$y[, ii]), lags = .ARord ),
                  mas[, ii, ] * WN.lags(y = cbind(extra$res[, ii]),
                              lags = .MAord ))) #+  extra$res
      
      
       if (residuals) {
         stop("Loglikelihood not implemented yet to",
              " handle residuals.") 
       } else {
         loglik.terms <- c(w) * dARMA(x = y, mean = y.est,
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
             .var.arg = var.arg , .realM1 = realM1 ))),
    
   
   
   
   
    validparams = eval(substitute(function(eta, y, extra = NULL) {
     
     M1  <- .realM1; n <- nrow(eta)
     NOS <- ncol(eta)/M1
     
     drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
       eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                 link =  .ldrift , earg = .edrift )
     
     sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder, drop = FALSE],
                        link = if ( .var.arg ) .lvar else .lsd ,
                        earg = if ( .var.arg ) .evar else .esd )
     
     ars <- array(NA_real_, dim = c(n, NOS, .ARord ))
     mas <- array(NA_real_, dim = c(n, NOS, .MAord ))
     
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
   
   
   
    #vfamily = c("vgtsff", "vgltsff-class", "ARMAff"),
    vfamily = c("ARMAXff", "ARMAvgltsmff"),
   
   
    
    deriv = eval(substitute(expression({
      
      M1  <- .realM1; n <- nrow(eta)
      NOS <- ncol(eta)/M1
      M   <- M1 * NOS
      
      drifts <- if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
        eta2theta(eta[, M1*(1:NOS) - .nOrder - 1, drop = FALSE],
                  link =  .ldrift , earg = .edrift )
      
      if ( .var.arg ) {
        sigWN <- sqrt(eta2theta(eta[, M1*(1:NOS) - .nOrder,
                                    drop = FALSE], .lvar , .evar ))
      } else {
        sigWN <- eta2theta(eta[, M1*(1:NOS) - .nOrder,
                               drop = FALSE], .lsd , .esd )
      }
      
      ars <- array(NA_real_, dim = c(n, NOS, .ARord ))
      mas <- array(NA_real_, dim = c(n, NOS, .MAord ))
      
      for (ai in 1:( .ARord ))
        ars[, , ai] <- eta2theta(eta[, M1 * (1:NOS) -
                        ( .ARord + .MAord ) + ai, drop = FALSE], 
                                 link = .lARcoeff , earg = .eARcoeff )
      
      for ( mi in 1:( .MAord ))
        mas[, , mi] <- eta2theta(eta[, M1 * (1:NOS) -
                        ( .MAord ) + mi , drop = FALSE],
                                link = .lMAcoeff , earg = .eMAcoeff )
      
      y.est <- matrix(NA_real_, n , NOS)
      for (ii in 1:NOS) 
        y.est[, ii] <- drifts[, ii] + rowSums(cbind(ars[, ii, ] *
                  WN.lags(y = cbind(extra$y[, ii]), lags = .ARord ),
                              mas[, ii, ] *
                WN.lags(y = cbind(extra$res[, ii]), lags = .MAord )))
      
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
      dsewn.deta <- 
        if (.var.arg) dtheta.deta(sigWN^2, .lvar , .evar) else
                                      dtheta.deta(sigWN, .lsd , .esd)
      dThe.deta  <- dtheta.deta(ars , .lARcoeff , .eARcoeff)
      dPhi.deta  <- dtheta.deta(mas , .lMAcoeff , .eMAcoeff)
      
      myderiv <- 
       c(w) * cbind(if ( .nodrift ) NULL else dl.drfmean * ddrif.deta,
        if ( .var.arg ) dl.dvar * dsewn.deta else dl.dsd * dsewn.deta,
                matrix( dl.dThe * dThe.deta, n, NOS * ( .ARord )),
                  matrix( dl.dPhi * dPhi.deta, n, NOS * ( .MAord )))
      
      myderiv <- 
         myderiv[, interleave.VGAM(( .realM1 )*NOS , M1 = .realM1)]
      
      myderiv
      
    }), list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
              .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .esd = esd , .evar = evar ,
              .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
              .ARord = ARord , .MAord = MAord , .nOrder = nOrder ,
              .var.arg = var.arg , .realM1 = realM1,
              .nodrift = nodrift ))),
   
   
   
   
   
    weight = eval(substitute(expression({
      
      if (( .type.EIM ) == "exact" ) {
        
        final.eim <- numeric(0)
        for (ii in 1:NOS) {
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
          
      final.eim <- cbind(final.eim, c(w[, ii]) * exact.eim * dthdeta)
    }
    
        
    final.eim <-
      final.eim[, interleave.VGAM( NOS * M1 * (M1 + 1) / 2, M1 = NOS)]
               final.eim <- array(c(final.eim), 
                                    dim = c(n, NOS, M1 * (M1 + 1)/ 2))
        
        wz <- arwz2wz(final.eim, M = M1 * NOS, M1 = M1)
        
      } else {
        
        ned2l.dsmn   <- 1 / sigWN^2
        ned2l.dvarSD <- if ( .var.arg ) 1 / (2 * sigWN^4) else 
                                                        2 / sigWN^2
        
        gammas.y <- apply(y, 2, function(x) {
          cross.gammas(x = x, lags = 0)
        })
        
        ned2l.dthe <- matrix(gammas.y, n, NOS * ( .ARord ), 
                  byrow = TRUE) / matrix(sigWN^2, n, NOS * ( .ARord ))
        
        dThe.deta <- matrix(dThe.deta, n, NOS * ( .ARord ))
        
        ned2l.dphi <- matrix(1.0, n, NOS * ( .MAord ))
        dPhi.deta  <- matrix(dPhi.deta, n , NOS * ( .MAord ))
        
        wz <- c(w) * cbind( if ( .nodrift ) NULL else 
             ned2l.dsmn * ddrif.deta^2, ned2l.dvarSD * dsewn.deta^2,
             ned2l.dthe * dThe.deta^2 , ned2l.dphi * dPhi.deta^2)
        
        wz <- wz[, interleave.VGAM( NOS * M1 , M1 = M1)]
        
      }
      
      
       wz[which(wz[, 1:(NOS * M1)] < .Machine$double.eps)] <-
                                          .Machine$double.base^(0.75)
       
       wz
      
    }), list( .var.arg = var.arg , .nOrder   = nOrder, 
              .ARord   = ARord   , .MAord   = MAord ,
              .realM1  = realM1  , .nodrift = nodrift ,
              .idrift = idrift ,
              .iARcoeff = iARcoeff , .iMAcoeff = iMAcoeff ,
              .ivar = ivar, .isd = isd ,
              .type.EIM = type.EIM )))
    ) # End of alo
  
  alo
  
} # End of ARMAXff
  





## ARMAXff control function.
ARMAXff.control <- function(save.weights = TRUE,
                            summary.HDEtest = FALSE,...) { 
  #criterion <- "loglikelihood"
  list(save.weights = save.weights,
       summary.HDEtest = summary.HDEtest,...)
}





quick.check.coeffs <- function(arc, mac, NOS, arOrd, maOrd) {
  
  if (length(arc) && (length(arc) != arOrd * NOS))
    stop("Conflicting number of initial values for the AR component.")
  
  if ( length(mac) && ( length(mac) != maOrd * NOS) )
    stop("Conflicting number of initial values for the MA component.")
  
  if (length(arc)) {
    
    arc <- matrix(arc, NOS, arOrd, byrow = TRUE)
    initChecks <- apply(arc, 1, function(x) {
         checkTS.VGAMextra(thetaEst = x,
                           tsclass  = "AR",
                           NofS     = 1,
                           chOrder  = arOrd,
                           pRoots   = FALSE, 
                           retmod   = TRUE)
      })
    if(any( initChecks < 1 + 1e-8))
      warning("Initial values for the AR coefficients for response ",
           " do not belong to a stationary model.")
  }
  
  
  
  if (length(mac)) {
    mac <- matrix(mac, NOS, maOrd, byrow = TRUE)
    initChecks <- apply(mac, 1, function(x) {
      checkTS.VGAMextra(thetaEst = x, tsclass  = "MA",
                        NofS     = 1, chOrder  = maOrd,
                        pRoots   = FALSE, retmod   = TRUE)
    })
    
    if(any( initChecks < 1 + 1e-8)) 
      warning("Initial values for the MA coefficients for response ",
              " do not belong to an invertible model.")
  }
  
  invisible(NULL)
}


