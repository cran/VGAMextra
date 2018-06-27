##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.
# Supports the Wald, score, and lrt tests (20180209)

MAqEIM.G2 <- function(y, mean, sdError, MAqcoeff,
                      var.arg = TRUE, order = 1,
                      nomean = FALSE, addRidge = 0.001) {
  
  Order <- order; rm(order)
  y <- cbind(y); nn <- nrow(y)
  MAqcoeff <- matrix(MAqcoeff, nrow = nn, ncol = Order)
  
  RMat  <- diag(c(sdError^2))
  dRMat <- cbind( if(!nomean) matrix(0, nrow = nn, ncol = 1) else NULL, 
                  if (var.arg) matrix(1, nrow = nn, ncol = 1) else
                    2 * sdError,  
                  matrix(0, nrow = nn, ncol = Order))
  
  y.lags <- WN.lags(y = y, lags = Order,
                    to.complete =  0.25 / (1 - MAqcoeff[1, ])^2)
  
  dmt0 <- matrix(1, nrow = nn, ncol = 1)
  dmt  <- cbind(if (!nomean) dmt0 else NULL, rep_len(0, nn), y.lags)
  
  M <- 2 + Order # For MAq 
  try.comb <- iam(NA, NA, M = M - nomean, both = TRUE)
  try.comb <- cbind(try.comb$row.index, try.comb$col.index)
  
  ### This is equation A.6, Porat and Friedlander (1986)
  finMat <- apply(try.comb, 1, function(x) {
    kl <- x
    invRMat <- diag(solve(RMat))
    dRdthel <- c(dRMat[, kl[1]])
    dRdthek <- c(dRMat[, kl[2]])
    term.1 <- (0.5) * invRMat * dRdthel * invRMat * dRdthek
    term.2 <- c(dmt[, kl[1]]) * invRMat * c(dmt[, kl[2]])
    
    term.1 + term.2
  })
  
  finMat[ , 1:M] <- (1 + addRidge) * finMat[ , 1:M]
  finMat
}





# MAq density
dMAq <- function(x, mean = 0, sd = 1, log = FALSE)
    dnorm(x = x, mean = mean, sd = sd, log = log)





### Family function MAXff
MAXff <-
  function(order     = 1,
           zero      = c(if (nomean) NULL else "Mean", "MAcoeff"),
           xLag      = 0,
           type.EIM  = c("exact", "approximate")[1],
           var.arg   = TRUE,
           nomean    = FALSE,
           noChecks  = FALSE,
           lmean     = "identitylink", 
           lsd       = "loge",
           lvar      = "loge",
           lMAcoeff  = "identitylink",
           imean     = NULL,
           isd       = NULL,
           ivar      = NULL,
           # Must be a VECTOR
           iMAcoeff  = NULL) {
    
    if (xLag < 0)
      stop("Wrong input for argument 'xLag'")
    
    if (any(order == 0))
      stop("Refer to uninormalff() to estimate the 2--parameter",
           "\n", " univariate normal distribution.")
    
    if ( !Is.Numeric(x = order, Nnegative = TRUE,
                     isInteger = TRUE))
      stop("Wrong input for argument 'order'.")
    Order <- order; rm(order)
    
    if ( length(imean) && !Is.Numeric(imean) )
      stop("Bad entry for argument 'imean'.") 
    if ( length(isd) && !Is.Numeric(isd, Nnegative = TRUE) )
      stop("Wrong input for argument 'isd'.") 
    if ( length(ivar) && !Is.Numeric(ivar, Nnegative = TRUE) )
      stop("Wrong input for argument 'ivar'.") 
    
    if ( !is.logical(var.arg) ||  !is.logical(nomean) )
      stop("Arguments 'var.arg' and 'nodrift' must \n",
           "be logical.")
    
    if (!is.logical(noChecks))
      stop("Wrong input for argument 'noChecks'")
    
    
    type.likelihood <- "exact"
    type.EIM  <- match.arg(type.EIM, 
                           c("exact", "approximate"))[1]
    
    eimMAff <- (type.EIM == "exact") 
    nodrift <- nomean
    idrift  <- imean
    ldrift  <- lmean
    
    
    ldrift  <- match.arg(ldrift, "identitylink")
    ldrift  <- as.list(substitute(ldrift))
    edrift  <- link2list(ldrift)
    ldrift  <- attr(edrift, "function.name")
    
    lsd <- as.list(substitute(lsd))
    esd <- link2list(lsd)
    lsd <- attr(esd, "function.name")
    
    lvar <- as.list(substitute(lvar))
    evar <- link2list(lvar)
    lvar <- attr(evar, "function.name")
    
    lMAcoeff <- as.list(substitute(lMAcoeff))
    eMAcoeff <- link2list(lMAcoeff)
    lMAcoeff <- attr(eMAcoeff,"function.name")
    
    
    BlurbMa <-blAr <- max(Order)[1]
    pre.blurb <- paste("MAcoeff", 1:blAr, sep ="")
    blurb.vec <- 
      c(namesof("Mean", ldrift, earg = edrift, tag = FALSE),
        if (var.arg)
          namesof("noiseVar", lvar, earg = evar, tag = FALSE)
        else
          namesof("noiseSD", lsd, earg = esd, tag = FALSE),
        namesof(pre.blurb, lMAcoeff , earg = eMAcoeff , tag = FALSE))
    
    blurb.vec2 <- NULL
    for (ii in 3:(blAr + 2)) 
      if (ii == (blAr + 2))
        blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii])) else
          blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii], ", "))
    
    alo <-
      new("vgltsmff", 
          blurb = c( "VGLTMs: ", ifelse( nodrift, blAr + 1, blAr + 2),
                      "-parameter MAX model of order-",blAr, "\nLinks: ",
                    #c(ifelse( nodrift, blAr + 1, blAr+2), 
                    #  "-parameter MAX model",
                    #  " of order-",blAr,".\n\n",
                    #  "Links:    "), 
                    blurb.vec[1], ", ",
                    blurb.vec[2], ", ",
                    blurb.vec2, ". \n",
                    c(ifelse(nomean, 
                    ifelse( BlurbMa == 1, 
                    paste("Model:    Y_t =   B^T * X_t + 
                          (MAcoeff1)*e_{t - 1}", "\n"),
                        paste("Model:  Y_t =   B^T * X_t +
                              (MAcoeff1)*e_{t - 1} + ",
                             "... +(MAcoeff",BlurbMa,")*e_{t - ",
                             BlurbMa,"} + e_{t} ,", "\n", sep = "")),
                      ifelse( BlurbMa == 1, 
                        paste("Model:    Y_t = Mean +  B^T * X_t +", 
                                "(MAcoeff1)*e_{t - 1} ", "\n", sep = ""),
                        paste("Model:    Y_t = Mean +  B^T * X_t +" , 
                              " (MAcoeff1)*e_{t - 1} + ",
                              "... + (MAcoeff",BlurbMa,")*e_{t - ",
                              BlurbMa,"} + e_{t} ,", "\n", sep = "")))),
                    "where e_{t} iid ~ N(0, noiseSD^2), and ",
                    "E[Y_t] = Mean", ". \n"),
                    #"Var[Y_t] = noiseSD^2 * ",
                    #c(ifelse(BlurbMa == 1,"(1 + MAcoeff1^2)" , 
                    #   paste("(1 + MAcoeff1^2 +...+ MAcoeff",
                     #         BlurbMa,"^2).", sep = "")))),
  #"\n", "Here, F_t denotes the information available at",
  #"time 't - 1'."),
    
    
    
    
    constraints = eval(substitute(expression({
      
      M1vec <-  2 + nOrder - .nomean
      constraints <- 
        cm.zero.VGAM(constraints, x = x, zero = .zero , M = M,
                     predictors.names = parameters.names,
                     #predictors.names = predictors.names,
                     M1 = M1vec)
      
    }), list( .zero = zero, .Order = Order, .nomean = nomean ))),
    
    
    
    
    
    infos = eval(substitute(function(...) {
      
      ## M1 was numeric(NA). Changed on 2018/01/19 to handle the
      ## Hauck - Donner effects. VGAM 1.0-5.
      
      list(M1       =  ( sum( .Order + 2) - .nodrift ) /length( .Order ),
           eimMAff  = .eimMAff ,
           Order    = .Order,
           Q1       = 1,
           expected = TRUE,
           multipleResponse = TRUE,
           ldrift   = .ldrift ,
           lsd      = .lsd ,
           lvar     = .lvar  ,
           edrift   = .edrift ,
           esd      = .esd ,
           evar     = .evar ,
           nomean   = .nomean ,
           lMAcoeff = .lMAcoeff ,
           eMAcoeff = .eMAcoeff ,
           type.likelihood = .type.likelihood ,
           type.EIM = .type.EIM ,
           zero  = .zero )
    }, list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
             .lMAcoeff = lMAcoeff ,
             .edrift = edrift , .esd = esd , .evar = evar ,
             .eMAcoeff = eMAcoeff , .type.EIM =type.EIM ,
             .type.likelihood = type.likelihood ,
             .Order = Order, .nodrift = nodrift ,
             .nomean = nomean ,
             .zero = zero , .eimMAff = eimMAff ))),
    
    
    first = eval(substitute(expression({
      # Not needed, at the moment... 2018-Jan
      M1 <- sum( .Order + 2) - .nodrift
      
      if ((NCOL(x) == 1) && (colnames(x) == "(Intercept)") && ( .xLag > 0))
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
      
    }), list( .Order = Order , .nodrift = nodrift , .xLag = xLag ))),
  
    
    initialize = eval(substitute(expression({
      
      check <- w.y.check(w = w, y = y,
                         Is.positive.y = FALSE,
                         ncol.w.max = Inf,
                         ncol.y.max = Inf,
                         out.wy     = TRUE,
                         colsyperw  = 1,
                         maximize   = TRUE)
      
      w <- check$w
      y <- extra$y <- check$y
      NOS <- extra$NOS <- ncoly <- ncol(y)
      n <- nrow(y)
      
      # Conformability condition. Changed on 2018/01/19 to handle the
      ## Hauck - Donner effects. VGAM 1.0-5.
      if ( length( .Order ) != NOS )
      stop("Non-conformable orders. ",
           "Enter the desired order per response.")
             
      nOrder <- rep( .Order , length = NOS)[1:NOS]
      extra$nOrder <- nOrder
      M1 <- (2 - .nodrift) + nOrder 
      
      for (jj in 1:NOS)
        w[1:nOrder[jj], jj] <- 1e0
      
      m.max  <- max(nOrder)
      m1.max <- max(M1)
      
      init.the <- vector("list", length = NOS)
      
      if ( length( .iMAcoeff ) ) {
        
        if( !Is.Numeric( .iMAcoeff ) )
          stop("'iMAcoeff' must be a vector with numeric entries.") 
        
        if ( length( .iMAcoeff ) != sum( nOrder ) )  
          stop("Invalid number of entries in 'iMAcoeff'.")
        
        auuxx <- 0
        for ( ii in 1:NOS ) {
          niAux <- .iMAcoeff[ auuxx + 1:nOrder[ii] ] 
          initChecks <- checkTS.VGAMextra(thetaEst = c( niAux ),
                                          tsclass  = "MA",
                                          NofS     = 1,
                                          chOrder  = nOrder[ii],
                                          pRoots   = FALSE, 
                                          retmod   = TRUE)
          if( any( initChecks < 1 + 1e-8 ) ) 
          warning("Initial values of the MA coefficients for response ",
                 ii, " do not belong to an invertible model. ", "\n",
                 "Get the summary() for more details.")
          
          init.the[[ii]] <- matrix(niAux, nrow = n , ncol = nOrder[ii],
                                   byrow = TRUE)
          auuxx <- auuxx + nOrder[ii]
        }
      }     
      
      names.1 <- if ( .nodrift ) NULL else
        paste("Mean", 1:NOS, sep = "")
      names.2 <- if ( .var.arg ) paste("noiseVar", 1:NOS, sep="") else
        paste("noiseSD", 1:NOS, sep = "")
      names.3 <- character(0)
      for (jj in 1:m.max) 
        names.3 <- c(names.3, paste(paste("MAcoeff", jj, sep = ""),
                                    1:NOS, sep = ""))
      
      pars.names <- c(names.1, names.2, names.3)
      pars.names <- pars.names[interleave.VGAM(m1.max * NOS, M1 = m1.max)]
      pars.names <- matrix(pars.names, NOS, ncol = m1.max, byrow = TRUE)
      
      pres.names <- c(if ( .nodrift ) NULL else
        namesof(names.1, .ldrift, .edrift, tag = FALSE),
        namesof(names.2, if ( .var.arg ) .lvar else .lsd ,
                if ( .var.arg ) .evar else .esd , tag = FALSE),
        namesof(names.3, .lMAcoeff , .eMAcoeff , tag = FALSE))
      
      pres.names <- pres.names[interleave.VGAM(m1.max * NOS, M1 = m1.max)]
      pres.names <- matrix(pres.names, NOS, ncol = m1.max, byrow = TRUE)
      parameters.names <- predictors.names <- character(0)
      
      for (ii in 1:NOS) {
        parameters.names <- c(parameters.names, pars.names[ii, 1:M1[ii]])
        predictors.names <- c(predictors.names, pres.names[ii, 1:M1[ii]])
      }
      
      nidrift <- rep( .idrift , NOS)[1:NOS]
      nisd    <- rep( .isd , NOS)[1:NOS]
      nivar   <- rep( .ivar , NOS)[1:NOS]
      
      init.sd <- if ( length ( nisd ) )
        matrix( nisd, nrow = n , ncol = NOS, byrow = TRUE ) else
          matrix( 1.0 , nrow = n , ncol = NOS) 
      
      init.var <- if ( length ( nivar ) )
        matrix( nivar , nrow = n , ncol = NOS, byrow = TRUE ) else
          matrix( 1.0 , nrow = n , ncol = NOS)
      
      if (!length(etastart)) {
        
        init.mean <- if ( length ( nidrift ) )
          matrix( nidrift, nrow = n , ncol = NOS, byrow = TRUE ) else
            matrix( colMeans(y) , nrow = n , ncol = NOS, byrow = TRUE)
        
        y.sc <- scale(y, center = TRUE, scale = FALSE)
        extra$res <- matrix(NA_real_, nrow = n, ncol = NOS)
        res2.mat  <-  matrix(NA_real_, nrow = n, ncol = NOS)
        
        for (rsp in 1:NOS) {
         to.fit  <- cbind(y.sc[, rsp],
                          WN.lags(y = y.sc[, rsp, drop = FALSE],
                                  lags = nOrder[rsp] + 1e1)) #[-(1:6), , drop = FALSE]
          ini.fit <- lsfit(x = to.fit[, -1, drop = FALSE],
                           y = to.fit[,  1, drop = FALSE],
                           intercept = FALSE)
          
          e.rsp <- cbind(y.sc[, rsp], #[-(1:6), rsp],
                         WN.lags(y = cbind(residuals(ini.fit)),
                                 lags = nOrder[rsp] ))
          
          #extra$res[, rsp] <- cbind(c(rep(NULL, 6), residuals(ini.fit)))
          #e.rsp   <- e.rsp[-(1:nOrder[rsp]), , drop = FALSE]
          
          ini.fit <- lsfit(x = e.rsp[, -1, drop = FALSE],
                           y = e.rsp[,  1, drop = FALSE],
                           intercept = FALSE)
          res2.mat[, rsp] <- residuals(ini.fit)
          
          if (!length( .iMAcoeff )) {
            ma.cfs <- c(coef(ini.fit),
                        if (length(coef(ini.fit)) < m.max)
                        rep(0.1, m.max - length(coef(ini.fit))) else NULL)
            init.the[[rsp]] <- matrix(ma.cfs, nrow = n,
                                      ncol = m.max, byrow = TRUE)
          }
          
          if (!length( nisd )) 
            init.sd[, rsp] <- sqrt(var(y[, rsp]) /
                                      (1 + sum(coef(ini.fit)^2)))  
          if (!length( nivar )) 
            init.var[, rsp] <- var(y[, rsp]) / (1 + sum(coef(ini.fit)^2))
        }
        
        extra$res <- res2.mat
        etastart <- cbind(if ( .nodrift ) NULL else 
          theta2eta(init.mean, .ldrift , earg = .edrift),
          if ( .var.arg )
            theta2eta(init.var, .lvar , earg = .evar ) else 
              theta2eta(init.sd, .lsd , earg = .esd ))
        
        pre.eta <- numeric(0)
        for (kk in 1:NOS) 
          pre.eta <- cbind(pre.eta, theta2eta(cbind(init.the[[kk]]),
                                          .lMAcoeff , earg = .eMAcoeff ))
        pre.eta  <- pre.eta[, interleave.VGAM(m.max * NOS, M1 = m.max,
                                              inverse = TRUE)]
      } else {
        pre.eta <- NULL
      }
     
      etastart <- cbind(etastart, pre.eta)
      etastart <- etastart[, interleave.VGAM(m1.max * NOS, M1 = m1.max)]
      pre.eta  <- array(etastart, dim = c(n, m1.max, NOS))
      etastart <- numeric(0)
      for (ii in 1:NOS) 
        etastart <- cbind(etastart, pre.eta[, , ii][, 1:M1[ii]])
      
      etastart
      
    }), list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
              .lMAcoeff = lMAcoeff,
              .edrift = edrift , .esd = esd , .evar = evar ,
              .eMAcoeff = eMAcoeff , .idrift = idrift ,
              .imean = imean, .ivar = ivar , .isd = isd ,
              .iMAcoeff = iMAcoeff , 
              .nodrift  = nodrift , .var.arg = var.arg ,
              .type.likelihood = type.likelihood ,
              .Order = Order ))),
    
    
    
    
    linkinv = eval(substitute(function(eta, extra = NULL ) {
      
      n   <- nrow(eta)
      NOS <- extra$NOS
      nOrder <- extra$nOrder
      ordMax <- max(nOrder)[1]
      M1 <- 2 + nOrder - .nodrift
      
      
      lkInv <- break.VGAMextra(eta      = eta, 
                               M1       = M1, 
                               noInter  = .nodrift ,
                               bOrder   = nOrder,
                               NOS      = NOS,
                               lInter   = .ldrift ,
                               lvar     = .lvar ,
                               lsd      = .lsd ,
                               lcoeff1  = .lMAcoeff ,
                               lcoeff2  = .lMAcoeff ,
                               typeTS   = "MA",
                               Complete = TRUE,
                               varArg   = .var.arg )
      maq.drMean <-  if (!( .nodrift )) lkInv[[1]] else
                              matrix(0.0, nrow = n, ncol = NOS)
      maq.Phi <- lkInv[[4]]
      update  <- FALSE
      
      if (update) {
        upd.rs <- numeric(0)
        for (jj in 1:NOS) 
          upd.rs <- cbind(upd.rs,  WN.lags(y = cbind(extra$res[, jj]),
                                           lags = ordMax))
        upd.rs    <- array(upd.rs * maq.Phi, dim = c(n, ordMax, NOS))
        extra$res <- extra$y - (maq.drMean +
                              apply(upd.rs, c(1, 3), function(x) sum(x)))
      }
      
      e.lags  <- array(0, dim = c(n, max(nOrder), NOS))
      for (ii in 1:NOS)
        e.lags[,  1:nOrder[ii], ii] <- 
              WN.lags(y = cbind(extra$res[, ii]), lags = nOrder[ii]) 
      
      e.lags <- matrix(e.lags, nrow = n, ncol = max(nOrder) * NOS)
      fitted.ts <- array(maq.Phi * e.lags, dim = c(n, max(nOrder), NOS))
      fitted.ts <- maq.drMean +
                apply(fitted.ts, c(1, 3), function(x) sum(x)) + extra$res
      
      fitted.ts
      
    }, list ( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
              .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .evar = evar , .esd = esd ,
              .eMAcoeff = eMAcoeff ,
              .var.arg = var.arg , .Order = Order ,
              .nodrift = nodrift ))),
          
          
          
          
    last = eval(substitute(expression({
      
      n   <- nrow(eta)
      NOS <- extra$NOS
      nOrder <- extra$nOrder
      M1  <- 2 - .nodrift + nOrder
      M   <- sum(M1)
      
      sllast <- 
        break.VGAMextra(eta     = eta, 
                        M1      = M1, 
                        noInter = .nodrift ,
                        NOS     = NOS,
                        bOrder  = nOrder,
                        lInter  = .ldrift ,
                        lvar    = .lvar ,
                        lsd     = .lsd ,
                        lcoeff1 = .lMAcoeff ,
                        lcoeff2 = .lMAcoeff ,
                        typeTS  = "MA",
                        varArg  = .var.arg )
      maq.drMean <- if (!( .nodrift )) sllast[[1]] else
                               matrix(0.0, nrow = n, ncol = NOS)
      maq.sd  <- sllast[[2]]
      maq.var <- sllast[[3]]
      
      misc$link <- c( rep( .lsd , M ) )
      names(misc$link) <- parameters.names
      misc$earg <- vector("list", length = M )
      names(misc$earg) <- parameters.names
      
      auxlast <- 0
      for ( ll in 1:NOS ) {
        if ( .nodrift ) {
          misc$link[ auxlast + 1 ]   <- if ( .var.arg ) .lvar else .lsd
          misc$earg[[auxlast + 1 ]]  <- if ( .var.arg ) .evar else .esd
        } else {
          misc$link[ auxlast + 1 ]   <- .ldrift
          misc$earg[[ auxlast + 1 ]] <- .edrift
          
          misc$link[ auxlast + 2 ]   <- if ( .var.arg ) .lvar else .lsd
          misc$earg[[ auxlast + 2 ]] <- if ( .var.arg ) .evar else .esd
        }
        
        auxlast <- auxlast + 2 - .nodrift 
        
        for ( ii in 1:nOrder[ll] ) {
          misc$link[ auxlast  + ii]    <- .lMAcoeff
          misc$earg[[ auxlast  + ii]]  <- .eMAcoeff
        }
        
        auxlast <- auxlast + nOrder[ll]
      }
      
      misc$var.arg   <- .var.arg
      misc$expected  <- TRUE
      misc$process   <- "MA"
      misc$Order     <- nOrder
      misc$zero      <- .zero
      misc$NOS       <- NOS
      misc$M1        <- M1
      misc$nomean    <- .nomean
      misc$flagAR   <- FALSE
      misc$eimMAff   <- .eimMAff
      misc$residuals <- extra$res
      misc$theta.names       <- parameters.names
      misc$control$wzepsilon <- control$wzepsilon
      misc$type.likelihood   <- .type.likelihood
      misc$multipleResponses <- TRUE
      fit$prior.weights <- w  # 20180315
      
      if ( !(.noChecks ) ) {
        
        saveNames <- names( fit$coefficients )
        Flag <- grep("(Intercept)", names(fit$coefficients) )
        
        if (length(Flag) == M) {
          fin.coe <- numeric(0)
          # At this point, No constraints set up
          myCoeffs <- coef(fit)[Flag]
          
          for (jj in 1:NOS) {
            pre.coe  <- myCoeffs[1:M1[jj]]
            myCoeffs <- myCoeffs[-(1:M1[jj])]
            if ( .nodrift ) {
              pre.coe <- eta2theta(pre.coe[-1], .lMAcoeff , .eMAcoeff )
            } else {
              pre.coe <- eta2theta(pre.coe[-(1:2)], .lMAcoeff , .eMAcoeff )
            }
            
            if (length(pre.coe) < max(nOrder))
              pre.coe <- c(pre.coe, rep(0, max(nOrder)- length(pre.coe)))
            
            fin.coe <- c(fin.coe, pre.coe)
          }
          
          MAroots <- checkTS.ffs(thetaEst = fin.coe,  # a vector.
                                 tsclass  = "MA",
                                 NofS     = NOS,
                                 pRoots   = FALSE,
                                 chOrder  = max(nOrder))
          
          flag <- (all( MAroots > 0 ) &&  any( MAroots <= 1 + 5e-3 ))
          # 2016/March/01. Checks on stationarity/ invertibility are 
          # now carried out by 'summary()'.
          if (!flag)
            cat("\nChecks on stationarity / invertibility successfully",
                "performed. \nNo roots lying inside the unit circle.",
                "\nFurther details within the 'summary' output.\n") else
                  
              cat("\nChecks on stationarity / invertibility successfully",
                  "performed. \nNOTE: Some roots lie inside the unit",
                  "circle. \nFurther details within the 'summary'", 
                  "output.\n")
          rm(flag)
        } else {
          cat("\nApparently, some constraints using 'cm.ARMA' have been",
              "set. \nDetails on stationarity / invertibility checks",
              "whitin the \n'summary' output.")
        } # End of (length(Flag) == M)
      }   # End ( !(.noChecks ) )
      
    }), list( .ldrift = ldrift, .lsd = lsd , .lvar = lvar ,
              .lMAcoeff = lMAcoeff , .nomean = nodrift ,
              .edrift = edrift , .esd = esd , .evar = evar , 
              .eMAcoeff = eMAcoeff , 
              .var.arg = var.arg , .nodrift = nodrift ,
              .type.likelihood = type.likelihood , .eimMAff = eimMAff ,
              .Order = Order , .zero = zero , .noChecks = noChecks ))),
    
    
    
    
    loglikelihood = eval(substitute(
            function(mu, y, w, residuals = FALSE, eta, 
                     extra = NULL, summation = TRUE) {
              
      NOS <- extra$NOS
      y   <- extra$y
      nOrder <- extra$nOrder
      M1  <- if (length(extra$M1)) extra$M1 else 2 - .nodrift + nOrder
      n   <- nrow(eta)
      ordMax <- max(nOrder)[1]
      
      slloglik <- break.VGAMextra(eta      = eta, 
                                  M1       = M1, 
                                  noInter  = .nodrift ,
                                  NOS      = NOS,
                                  bOrder   = nOrder,
                                  lInter   = .ldrift ,
                                  lvar     = .lvar ,
                                  lsd      = .lsd ,
                                  lcoeff1  = .lMAcoeff ,
                                  lcoeff2  = .lMAcoeff ,
                                  typeTS   = "MA",
                                  Complete = TRUE,
                                  varArg   = .var.arg )
      maq.drMean <- if (!( .nodrift ))  slloglik[[1]] else 
        matrix(0.0, nrow = n, ncol = NOS) 
      maq.sd  <- slloglik[[2]]
      maq.var <- slloglik[[3]]
      maq.Phi <- slloglik[[4]]
      
      e.lags  <- array(0, dim = c(n, max(nOrder), NOS))
      
      for (ii in 1:NOS)
        e.lags[,  1:nOrder[ii], ii] <- 
        WN.lags(y = cbind(extra$res[, ii]), lags = nOrder[ii]) 
      
      e.lags <- matrix(e.lags, nrow = n, ncol = max(nOrder) * NOS)
      
      mean.ts <- array(maq.Phi * e.lags, dim = c(n, max(nOrder), NOS))
      mean.ts <- maq.drMean + apply(mean.ts, c(1, 3), sum) #+ extra$res 
      names(mean.ts) <- NULL
      
      if (residuals) {
        stop("Loglikelihood not implemented yet", 
             " to handle residuals.")
      } else {
        loglik.terms <-
          c(w) * dMAq(x = y, mean = mean.ts,
                       sd = maq.sd, log = TRUE)
        sum( loglik.terms )
      }
      
    }, list (  .ldrift = ldrift , .lsd = lsd , .lvar = lvar ,
               .lMAcoeff = lMAcoeff , 
               .edrift = edrift , .esd = esd , .evar = evar ,
               .MAcoeff = eMAcoeff ,
               .var.arg = var.arg , .Order = Order ,
               .type.likelihood = type.likelihood ,
               .nodrift = nodrift ))),
    
    
    
    
    validparams = eval(substitute(function(eta, y, extra = NULL) {
      
      NOS <- ncol(y)
      n   <- nrow(eta)
      nOrder <- extra$nOrder
      ordMax <- max(nOrder)[1]
      M1  <- 2 - .nodrift + nOrder
      
      slloglik <- 
        break.VGAMextra(eta      = eta, 
                        M1       = M1, 
                        noInter  = .nomean ,
                        NOS      = NOS,
                        bOrder   = nOrder,
                        lInter   = .ldrift ,
                        lvar     = .lvar ,
                        lsd      = .lsd ,
                        lcoeff1  = .lMAcoeff ,
                        lcoeff2  = .lMAcoeff ,
                        typeTS   = "MA",
                        Complete = TRUE,
                        varArg   = .var.arg )
      maq.drMean <- if (!(.nodrift)) slloglik[[1]] else 
        matrix(0.0, nrow = n, ncol = NOS) 
      maq.sd  <- slloglik[[2]]
      maq.var <- slloglik[[3]]
      maq.Phi <- slloglik[[4]]
      
      okay1 <- (all(is.finite(maq.sd)) && all(0 < maq.sd) &&
                  all(is.finite(maq.var)) && all(0 < maq.var) &&
                  all(is.finite(maq.drMean)) && all(is.finite(maq.Phi)))
      okay1
      
    }, list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
             .lMAcoeff = lMAcoeff , 
             .edrift = edrift , .esd = esd , .evar = evar ,
             .eMAcoeff = eMAcoeff ,
             .var.arg = var.arg , .Order = Order ,
             .type.likelihood = type.likelihood ,
             .nodrift = nodrift , .nomean = nomean ))),
    
    
    
    #vfamily = c("vgtsff", "vgltsff-class", "MAff"),
     vfamily = c("MAXff", "ARMAvgltsmff"),
    
    
    
    deriv = eval(substitute(expression({
      
      NOS <- ncol(y)
      n   <- nrow(eta)
      nOrder <- extra$nOrder
      M1  <- 2 - .nodrift + nOrder
      M   <- sum(M1)
      ordMax <- max(nOrder)[1]
      
      slderiv <- 
        break.VGAMextra(eta      = eta, 
                        M1       = M1, 
                        noInter  = .nomean ,
                        NOS      = NOS,
                        bOrder   = nOrder,
                        lInter   = .lmean ,
                        lvar     = .lvar ,
                        lsd      = .lsd ,
                        lcoeff1  = .lMAcoeff ,
                        lcoeff2  = .lMAcoeff ,
                        typeTS   = "MA",
                        Complete = TRUE, 
                        varArg   = .var.arg )
      
      maq.drMean <-  if (!(.nodrift )) slderiv[[1]] else 
                       matrix(0.0, nrow = n, ncol = NOS)
      maq.sd  <- slderiv[[2]]
      maq.var <- slderiv[[3]]
      maq.Phi <- pre.Phi <- slderiv[[4]]
      pre.Phi <- array(pre.Phi, dim = c(n, ordMax, NOS))
      
      e.lags  <- array(0, dim = c(n, max(nOrder), NOS))
      for (ii in 1:NOS) 
        e.lags[,  1:nOrder[ii], ii] <- 
                 WN.lags(y = cbind(extra$res[, ii]), lags = nOrder[ii]) 
      
      e.lags <- matrix(e.lags, nrow = n, ncol = max(nOrder) * NOS)
      #y.lags <- matrix(y.lags, nrow = n, ncol = max(nOrder) * NOS)
      y.means <- array(maq.Phi * e.lags, dim = c(n, max(nOrder), NOS))
      y.means <- y - (maq.drMean + apply(y.means, c(1, 3), sum))
      
      dl.drfmean <- y.means / maq.var
      if ( .var.arg ) {
        dl.dvar <- y.means^2 / (2 * maq.var^2) - 1 / (2 * maq.var)
      } else {
        dl.dsd <- y.means^2 / maq.sd^3 - 1 / maq.sd
      }
      
      e.lags  <- array(e.lags, dim = c(n, max(nOrder), NOS))
      #y.lags <- array(y.lags, dim = c(n, max(nOrder), NOS))
      dl.dPhi <- array(NA_real_, dim = c(n, max(nOrder), NOS))
      for (ii in 1:NOS)
        dl.dPhi[, , ii] <- (array( y.means, dim = c(n, 1, NOS))[, , ii] *
                              e.lags[, , ii]) / maq.var[, ii]
      
      dl.dPhi <- interleaveArray.VGAMextra(dl.dPhi)
      
      # Partial derivatives of theta wrt eta #
      drfm.deta <- dtheta.deta(maq.drMean, .lmean , earg = .emean )
      if ( .var.arg ) {
        dvar.deta <- dtheta.deta(maq.var , .lvar , earg = .evar )
        dVarSD    <- dvar.deta
      } else {
        dsd.deta  <- dtheta.deta(maq.sd  , .lsd  ,  earg = .esd )
        dVarSD    <- dsd.deta
      }
      
      dPhidEta <- interleaveArray.VGAMextra(array(maq.Phi,
                                    dim = c(n, ordMax, NOS)))
      dPhidEta <- dtheta.deta(dPhidEta, .lMAcoeff , earg = .eMAcoeff)
      
      # dlog / deta
      dl.deta1 <- c(w) * dl.drfmean * drfm.deta
      if ( .var.arg ){
        dl.deta2 <- c(w) * dl.dvar * dvar.deta
      } else {
        dl.deta2 <- c(w) * dl.dsd * dsd.deta
      }
      
      a.help  <- interleaveArray.VGAMextra(c(w) * dl.dPhi * dPhidEta,
                                           inverse = TRUE)
      fDeriv <- numeric(0)
      for (ii in 1:NOS) {
        fDeriv <- cbind(fDeriv, if ( .nodrift ) NULL else  dl.deta1[, ii],
                        dl.deta2[, ii],
                        matrix(c(a.help[, 1:nOrder[ii] , ii]), nrow = n,
                               ncol = nOrder[ii]))
      }
      
      fDeriv
      
    }), list( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
              .lMAcoeff  = lMAcoeff , .lmean = ldrift , .emean = edrift ,
              .edrift = edrift , .evar = evar , .esd = esd ,
              .eMAcoeff = eMAcoeff , 
              .var.arg = var.arg , .Order = Order , 
              .nodrift = nodrift , .nomean = nomean ))),
    
    
    
    
    weight = eval(substitute(expression({
      
      dPhidEta <- interleaveArray.VGAMextra(dPhidEta, inverse = TRUE)
      
      ### Exact EIMs.
      if ( .eimMAff ) {
        exact.eim <- vector("list", NOS); col.count <- 0
        for (jj in 1:NOS) {
          exact.eim[[jj]] <-  MAqEIM.G2(y  = extra$res[, jj],
                              mean  = maq.drMean[, jj, drop = FALSE],
                              sdError  = maq.sd[, jj, drop = FALSE],
                              MAqcoeff = pre.Phi[, 1:nOrder[jj], jj],
                              var.arg  = .var.arg ,
                              order    = nOrder[jj],
                              nomean   = .nomean )

          comb.wz <- combVGAMextra(1:M1[jj], nodrift = FALSE)
          a.bind  <- cbind(if ( .nodrift ) NULL else drfm.deta[, jj],
                           dVarSD[, jj],
                           dPhidEta[, 1:nOrder[jj], jj])
          
          dthdeta <- apply(comb.wz, 1, function(x) {
            a.bind[, x[1]] * a.bind[, x[2]]
          })
          
          exact.eim[[jj]] <- c(w[, jj]) * exact.eim[[jj]] * dthdeta
          col.count <- col.count + ncol(exact.eim[[jj]])
        }
        
        if (NOS == 1) {
          wz <- exact.eim[[1]]
        } else {
          if (max(nOrder) == min(nOrder)) {
            wz <- matrix(unlist(exact.eim), nrow = n, ncol = col.count)
            wz <- wz[, interleave.VGAM(col.count, M1 = col.count / NOS,
                                       inverse = TRUE)]
            wz <- array(wz, dim = c(n, NOS ,col.count / NOS))
            wz <- arwz2wz(wz, M = M1[1] * NOS, M1 = M1[1])
          } else {
            wz <- numeric(0)
            for (jj in 1:NOS) 
              wz <- cbind(wz, exact.eim[[jj]][, 1:M1[jj]])
          } 
        }
        
      } else {
        ned2l.dsmn   <- (1 / maq.var) 
        ned2l.dvarSD <- if ( .var.arg ) 1/(2 * maq.var^2) else 2 / maq.var
        wz <- cbind( if ( .nodrift ) NULL else 
                  ned2l.dsmn * drfm.deta^2, ned2l.dvarSD * dVarSD^2)
        
        comp.ders <- vector("list", NOS)
        for (ii in 1:NOS) {
          comp.ders[[ii]] <- WN.lags(y = cbind(extra$res[, ii]), 
                                     lags = ordMax, #nOrder[ii],
                                     to.complete = rep(1, ordMax))
        }
        comp.ders <- matrix(unlist(comp.ders)^2, n, NOS * ordMax) / 
                          matrix(maq.var, nrow = n, ncol = NOS * ordMax)
        
        # Approximate EIMs built up from E(e^2)/sigma^2
        comp.ders <- matrix(1, n, NOS * ordMax) 
          #matrix(maq.var, nrow = n, ncol = NOS * ordMax)
        
        dPhi.bis <- matrix(dPhidEta, nrow = n, ncol = ordMax * NOS)
        
        if (FALSE) { # keep for ARMA
          gammas.y <- apply(y, 2, function(x) {
            cross.gammas(x = x, lags = 0)
          })
          
          comp.ders <- matrix(c(gammas.y), nrow = n, ncol = NOS * ordMax,
                              byrow = TRUE) / matrix(maq.var, nrow = n,
                                                     ncol = NOS * ordMax)
        }
        
        comp.ders <- comp.ders * matrix(dPhi.bis^2, nrow = n,
                                        ncol = ordMax * NOS)
        
        comp.ders <- comp.ders[, interleave.VGAM(ncol(comp.ders),
                                            M1 = ordMax, inverse = FALSE)]
        
        wz <- c(w) * cbind(wz, comp.ders)
        pre.wz <- wz[, interleave.VGAM(ncol(wz), M1 = max(M1))]
        pre.wz <- array(pre.wz, dim = c(n, max(M1), NOS))
        wz <- numeric(0)
        for (jj in 1:NOS) 
          wz <- cbind(wz, pre.wz[, 1:M1[jj], jj] )
        wz[wz < 1e-10] <- 1e-8
      }
      
      wz
      
    }), list( .var.arg  = var.arg , .Order   = Order , 
              .nodrift  = nodrift , .eimMAff = eimMAff ,
              .idrift = idrift , .nomean = nodrift ,
              .ivar = ivar , .isd = isd , 
              .iMAcoeff = iMAcoeff ))) ) # End of new("vgltsmff") 
    
    alo
} # End of MAXff






## MAXff control function.
MAXff.control <- function(save.weights = TRUE,
                          summary.HDEtest = FALSE,
                          ...) { 
  list(save.weights = save.weights,
       summary.HDEtest = summary.HDEtest,...)
}


# MAqEIM.G2 back up using matrix algebra.
if (FALSE)
  MAqEIM.G2extra <- function(y, mean, sdError, MAqcoeff, var.arg = TRUE,
                             order = 1, nomean = FALSE) {
    
    Order <- order; rm(order)
    y <- cbind(y); nn <- nrow(y)
    MAqcoeff <- matrix(MAqcoeff, nrow = nn, ncol = Order)
    
    RMat  <- diag(c(sdError^2))
    dRMat <- cbind( if(!nomean) matrix(0, nrow = nn, ncol = 1) else NULL, 
                    if (var.arg) matrix(1, nrow = nn, ncol = 1) else
                      2 * sdError,  
                    matrix(0, nrow = nn, ncol = Order))
    
    y.lags <- WN.lags(y = y, lags = Order,
                      to.complete =  1 / (1 - MAqcoeff[1, ])^2)
    
    dmt0 <- matrix(1, nrow = nn, ncol = 1)
    dmt  <- cbind(if (!nomean) dmt0 else NULL, rep_len(0, nn), y.lags)
    
    M <- 2 + Order # For MAq 
    try.comb <- iam(NA, NA, M = M - nomean, both = TRUE)
    try.comb <- cbind(try.comb$row.index, try.comb$col.index)
    
    ### This is equation A.6, Porat and Friedlander (1986)
    finMat <- apply(try.comb, 1, function(x) {
      kl <- x
      invRMat <- solve(RMat)
      dRdthel <- diag(dRMat[, kl[1]])
      dRdthek <- diag(dRMat[, kl[2]])
      term.1 <- (0.5) * invRMat %*% dRdthel %*% invRMat %*% dRdthek
      term.2 <- diag(dmt[, kl[1]]) %*%
        invRMat %*% diag(dmt[, kl[2]])
      diag(term.1 + term.2)
    })
    
    ## Return the EIMs
    finMat[ , -(1:M)] <- (0.99) * finMat[ , -(1:M)]
    finMat
  }
