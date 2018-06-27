##########################################################################
# These functions are 
# Copyright (C) 2014-2018 V. Miranda & T. W. Yee, University of Auckland.
# All rights reserved.
# Supports the Wald, score, and lrt tests (20180209)

ARpEIM.G2 <- function(y, drift, sdError, ARpcoeff = NULL,
                      var.arg = TRUE, order = 1,
                      nodrift = FALSE, addRidge = 0.001) {
  
  Order <- order; rm(order)
  y <- cbind(y); nn <- nrow(y)
  ARpcoeff <-
    if (!length(ARpcoeff)) matrix(0.5, nrow = nn, ncol = Order) else
      matrix(ARpcoeff, nrow = nn, ncol = Order)
  
  RMat  <- diag(c(sdError^2))
  dRMat <- cbind( if(!nodrift) matrix(0, nrow = nn, ncol = 1) else NULL, 
                  if (var.arg) matrix(1, nrow = nn, ncol = 1) else
                    2 * sdError,  
                  matrix(0, nrow = nn, ncol = Order))
  
  y.lags <- WN.lags(y = y, lags = Order,
                    to.complete = 1 / (1 - ARpcoeff[1, ])^2)
  
  dmt0 <- matrix(1, nrow = nn, ncol = 1)
  dmt  <- cbind(if (!nodrift) dmt0 else NULL, rep_len(0, nn), y.lags)
  
  M <- 2 + Order  # For ARp
  try.comb <- iam(NA, NA, M = M - nodrift, both = TRUE)
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





# Computes the positions for every entry in the EIMs conforming
# with VGAM.
combVGAMextra <- function(x, nodrift = FALSE) {
  
  x.lth <- length(x)
  if (nodrift) {
    x <- x[-x.lth]
    x.lth <- length(x)
  }
  
  if (length(x) < 2)
    return(x)
  
  my.list <- apply(cbind(1:x.lth), 1, function (ii) {
    x.1 <- x[1:(x.lth - ii + 1)]
    x.2 <- x[ii:x.lth]
    cbind(x.1, x.2)
  })
  
  ret.mat <- my.list[[1]]
  for (ll in 2:x.lth)
    ret.mat <- rbind(ret.mat, my.list[[ll]])
  ret.mat
}





# Computes a matrix with the lagged values for every y[, jj].
WN.lags <- function(y, lags, to.complete = NULL) {
  
  #if (!(is.data.frame(y) || is.matrix(y)) )
  #  stop("'y' must be a data.frame or matrix.")
  if (NCOL(y) > 1)
    stop("Univariates data frames handled only.")
  if (lags <= 0 || lags%%1 != 0)
    stop("Invalid input for argument 'lags'.")
  
  nn <- NROW(y)
  y <- cbind(c(y))
  ret.matrix <- matrix(0, nrow = nn, ncol = lags)
  fill.t <- to.complete
  for (ll in 1:lags) {
    ret.matrix[-(1:ll), ll] <- y[1:(nn - ll), 1]
    # 2017/03/31 Old: ret.matrix[1:ll, ll] <- mean(y[1:(nn - ll), 1])
    ret.matrix[1:ll, ll] <- if (length(fill.t)) fill.t[1:ll] else
      0 * mean(y[1:(nn - ll), 1])
  }
  rm(fill.t)
  
  ret.matrix
}





# Useful for AR1extra() only, to arrange the EIMs accordingly.
# This is done manually in ARpff.
arwzTS <- function(wz, w, M1, dTHE.dETA, print.EIM = FALSE) {
  
  if (length(dim(wz)) != 3)
    stop("'wz' must be an array")
  
  nn    <- dim(wz)[1]
  dim.y <- dim(wz)[2]
  nresp <- dim(wz)[3]
  
  if (nresp == 1) {
    return(dTHE.dETA * c(w) *
             matrix(wz, nrow = nn, ncol = dim.y, byrow = FALSE))
  }
  
  fin.wz <- matrix(c(wz), nrow = nn, ncol = nresp * dim.y, byrow = FALSE)
  fin.wz <- fin.wz[, interleave.VGAM(nresp * dim.y, M1 = dim.y,
                                     inverse = TRUE)]
  
  fin.wz <- dTHE.dETA * matrix(c(w), nrow = nn, ncol = nresp * dim.y,
                               byrow = FALSE) * fin.wz
  if ( print.EIM ) {
    fin.wz <- fin.wz[, interleave.VGAM(nresp * dim.y, M1 = dim.y,
                                       inverse = FALSE)]
    fin.wz <- array(fin.wz, dim = c(nn, dim.y, nresp))
    return(fin.wz)
  }
  
  fin.wz <- array(fin.wz, dim = c(nn, nresp, dim.y))
  
  ## Here, M1 = number of parameters, .e.g M1 = 3 - nodrift in the AR1
  fin.wz <- arwz2wz(fin.wz, M = M1 * nresp, M1 = M1)
  fin.wz
}





# Same as interleave.VGAM() but applied to arrays.
interleaveArray.VGAMextra <- function(x, inverse = FALSE) {
  if (!is.array(x))
    stop("'x' must be an array")
  
  if (!inverse) {
    nn    <- dim(x)[1]
    dim.x <- dim(x)[2]
    nresp <- dim(x)[3]
    
    fin.x <- matrix(x, nrow = nn, ncol = nresp * dim.x, byrow = FALSE)
    fin.x <- fin.x[, interleave.VGAM(nresp * dim.x, M1 = dim.x,
                                     inverse = TRUE)]
    
    return(array(fin.x, dim = c(nn, nresp, dim.x)))
  } else {
    nn    <- dim(x)[1]
    nresp <- dim(x)[2]
    dim.x <- dim(x)[3]
    
    fin.x <- matrix(x, nrow = nn, ncol = nresp * dim.x, byrow = FALSE)
    fin.x <- fin.x[, interleave.VGAM(nresp * dim.x, M1 = dim.x,
                                     inverse = FALSE)]
    
    return(array(fin.x, dim = c(nn, dim.x, nresp)))
    
  }
}





# ARp density
dARp <- function(x, mean = 0, sd = 1, log = FALSE)
  dnorm(x = x, mean = mean, sd = sd, log = log)





### Family function ARXff
ARXff <- 
  function(order     = 1,
           zero      = c( if (nodrift) NULL else "ARdrift", "ARcoeff"),
           xLag      = 0,
           type.EIM  = c("exact", "approximate")[1],
           var.arg   = TRUE,
           nodrift   = FALSE,
           noChecks  = FALSE,
           ldrift    = "identitylink", 
           lsd       = "loge",
           lvar      = "loge",
           lARcoeff  = "identitylink",
           idrift    = NULL,
           isd       = NULL,
           ivar      = NULL,  
           # Must be a VECTOR
           iARcoeff  = NULL) {
    
    
    if (xLag < 0)
      stop("Wrong input for argument 'xLag'")
    
    if (any(order == 0))
      stop("Refer to uninormalff() to estimate the 2--parameter",
           "\n", " univariate normal distribution.")
    
    if ( !Is.Numeric(x = order,
                     Nnegative = TRUE,
                     isInteger = TRUE))
      stop("Wrong input for argument 'order'.")
    Order <- order; rm(order)
    
    if ( length(idrift) && !Is.Numeric(idrift) )
      stop("Bad entry for argument 'idrift'.") 
    if ( length(isd) && !Is.Numeric(isd, Nnegative = TRUE) )
      stop("Wrong input for argument 'isd'.") 
    if ( length(ivar) && !Is.Numeric(ivar, Nnegative = TRUE) )
      stop("Wrong input for argument 'ivar'.") 
    
    if ( !is.logical(var.arg) ||  !is.logical(nodrift) )
      stop("Arguments 'var.arg' and 'nodrift' must \n",
           "be logical.")
    
    if (!is.logical(noChecks))
      stop("Wrong input for argument 'noChecks'")
    
    type.likelihood <- "exact"
    type.EIM  <- match.arg(type.EIM, 
                           c("exact", "approximate"))[1]
    
    eimARff  <- (type.EIM == "exact") 
    
    ldrift <- match.arg(ldrift, "identitylink")
    ldrift <- as.list(substitute(ldrift))
    edrift <- link2list(ldrift)
    ldrift <- attr(edrift, "function.name")
    
    lsd <- as.list(substitute(lsd))
    esd <- link2list(lsd)
    lsd <- attr(esd, "function.name")
    
    lvar <- as.list(substitute(lvar))
    evar <- link2list(lvar)
    lvar <- attr(evar, "function.name")
    
    lARcoeff <- as.list(substitute(lARcoeff))
    eARcoeff <- link2list(lARcoeff)
    lARcoeff <- attr(eARcoeff,"function.name")
    
    blAr <- max(Order)
    pre.blurb <- paste("ARcoeff", 1:blAr, sep ="")
    blurb.vec <- 
      c(namesof("ARdrift", ldrift, earg = edrift, tag = FALSE),
        if (var.arg)
          namesof("noiseVar", lvar, earg = evar, tag = FALSE)
        else
          namesof("noiseSD", lsd, earg = esd, tag = FALSE),
        namesof(pre.blurb, lARcoeff , earg = eARcoeff , tag = FALSE))
    
    blurb.vec2 <- NULL
    for (ii in 3:(blAr + 2)) 
      if (ii == (blAr + 2))
        blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii])) else
          blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii], ", "))
    
    alo <- 
      new("vgltsmff", 
          blurb = c("VGLTMs: ", ifelse( nodrift, blAr + 1, blAr + 2),
                    "-parameter ARX model of order-",blAr,".\n\n",
                      "Links:    ", 
                    blurb.vec[1], ", ",
                    blurb.vec[2], ", ",
                    blurb.vec2, ". \n",
                    c("Model:    Y_t = ",
                      if (nodrift) NULL else "ARdrift + B^T * X_t +",
                      " (ARcoeff1)*Y_{t - 1} + 
                      ... +(ARcoeff",blAr,")*Y_{t - ",
                      blAr,"} + e_{t} ,", "\n"),
                      "where e_{t} ~ N(0, noiseSD^2), and ", "\n"),
                      #"E[Y_t]   = ARdrift/ (1 - ARcoeff1 -...- ARcoeff",
                      #blAr, ")"
                      #),
                    #"\n",
                    #"Var[Y_t ] = ARCoeff1*gamma1 +...+ ARCoeff", 
                    #blAr,"*gamma",blAr, 
                    #"+ noiseSD^2,", "\n",
                    #"where gamma{k} = Cov(Y_{t}, Y_{t+k}) =", 
                    #"Cov(Y_{t}, Y_{t-k})."
                    #),
    
    
    
    
    
    first = eval(substitute(expression({
      
      y <- as.matrix(y)
      mat.res <- matrix(NA_real_, NROW(y), NCOL(y))
      for (ii in 1:(NCOL(y))) 
        mat.res[, ii] <- residuals(lsfit(x = WN.lags(y = cbind(y[, ii]),
                                                lags = max( .Order , 10)),
                                         y = cbind(y[, ii, drop = FALSE]),
                                         intercept = TRUE))
      extra$res <- mat.res
      
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
      
      
    }),list( .Order = Order , .xLag = xLag ))),
    
    
    
    
    constraints = eval(substitute(expression({
      
      M1vec <-  2 + nOrder - .nodrift
      constraints <- 
        cm.zero.VGAM(constraints, x = x, zero = .zero , M = M,
                     predictors.names = parameters.names,
                     #predictors.names = predictors.names,
                     M1 = M1vec)
    }), list( .zero = zero, .Order = Order, .nodrift = nodrift ))),
          
          
          
          
          
    infos = eval(substitute(function(...) {
      
      ## M1 was numeric(NA). Changed on 2018/01/19 to handle the
      ## Hauck - Donner effects. VGAM 1.0-5.
      list(M1       = ( sum( .Order + 2) - .nodrift ) /length( .Order ),
           eimARff  = .eimARff ,
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
           nodrift  = .nodrift ,
           lARcoeff = .lARcoeff ,
           eARcoeff = .eARcoeff ,
           type.likelihood = .type.likelihood ,
           type.EIM = .type.EIM ,
           zero  = .zero )
    }, list(.ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
            .lARcoeff = lARcoeff ,
            .edrift = edrift , .esd = esd , .evar = evar , 
            .eARcoeff = eARcoeff , .type.EIM = type.EIM ,
            .type.likelihood = type.likelihood ,
            .Order = Order, .nodrift = nodrift ,
            .zero = zero , .eimARff = eimARff ))),
          
          
          
          
          
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
      
      # Conformability condition Changed on 2018/01/19 to handle the
      ## Hauck - Donner effects. VGAM 1.0-5.
      
      if ( length( .Order ) != NOS )    #if ( length( .Order ) > NOS )
        stop("Non-conformable orders. ",
             "Enter the desired order per response.")
      
      nOrder <- rep( .Order , length = NOS)[1:NOS]
      extra$nOrder <- nOrder
      M1 <- (2 - .nodrift) + nOrder 
      
      for (jj in 1:NOS)
        w[1:nOrder[jj], jj] <- 1e0
      m.max  <- max(nOrder)
      m1.max <- max(M1)
      
      if ( length( .iARcoeff ) ) {
        
        if( !Is.Numeric( .iARcoeff ) )
          stop("'iARcoeff' must be a vector with numeric entries.") 
        
        if ( length( .iARcoeff ) != sum( nOrder ) )  
          stop("Invalid number of entries in 'iARcoeff'.")
        
        init.the <- vector("list", length = NOS)
        auuxx <- 0
        for ( ii in 1:NOS ) {
          niAux <- .iARcoeff[ auuxx + 1:nOrder[ii] ] 
          initChecks <- checkTS.ffs(thetaEst = c( niAux ),
                                    tsclass  = "AR",
                                    NofS     = 1,
                                    chOrder  = nOrder[ii],
                                    pRoots   = FALSE, 
                                    retmod   = TRUE)
          if( any( initChecks < 1 + 1e-2 ) ) 
            stop("Initial values of the AR coefficients for response ",
                 ii, " do not belong to a stationary model.")
          
          init.the[[ii]] <- matrix(niAux, nrow = n , ncol = nOrder[ii],
                                   byrow = TRUE)
          auuxx <- auuxx + nOrder[ii]
        }
      }     
      
      names.1 <- if ( .nodrift ) NULL else
        paste("ARdrift", 1:NOS, sep = "")
      names.2 <- if ( .var.arg ) paste("noiseVar", 1:NOS, sep="") else
        paste("noiseSD", 1:NOS, sep = "")
      names.3 <- character(0)
      for (jj in 1:m.max) 
        names.3 <- c(names.3, paste(paste("ARcoeff", jj, sep = ""),
                                    1:NOS, sep = ""))
      
      pars.names <- c(names.1, names.2, names.3)
      pars.names <- pars.names[interleave.VGAM(m1.max * NOS, M1 = m1.max)]
      pars.names <- matrix(pars.names, NOS, ncol = m1.max, byrow = TRUE)
      
      pres.names <- c(if ( .nodrift ) NULL else
        namesof(names.1, .ldrift, .edrift, tag = FALSE),
        namesof(names.2, if ( .var.arg ) .lvar else .lsd ,
                if ( .var.arg ) .evar else .esd , tag = FALSE),
        namesof(names.3, .lARcoeff , .eARcoeff , tag = FALSE))
      
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
      
      if ( !length( .iARcoeff ))
        init.the <- vector("list", length = NOS)
      
      if (!length(etastart)) {
        
        init.mean <- if ( length ( nidrift ) )
          matrix( nidrift, nrow = n , ncol = NOS, byrow = TRUE ) else
            matrix( 0.01 , nrow = n , ncol = NOS) 
        
        init.sd <- if ( length ( nisd ) )
          matrix( nisd, nrow = n , ncol = NOS, byrow = TRUE ) else
            matrix( 1.0 , nrow = n , ncol = NOS) 
        
        init.var <- if ( length ( nivar ) )
          matrix( nivar , nrow = n , ncol = NOS, byrow = TRUE ) else
            matrix( 1.0 , nrow = n , ncol = NOS) 
        
        for (spp. in 1:NOS) {
          GammaS <- cross.gammas(x = y[, spp. ], lags = nOrder[spp. ])
          myToeplitz <- toeplitz(GammaS[1:nOrder[ spp. ] ] )
          myresp     <- c(GammaS[2:(nOrder[ spp. ] + 1)])
          initTheta  <- justHelp <- solve(myToeplitz, myresp)
          if (length(initTheta) < m.max)
            initTheta <- c(initTheta, rep(0, m.max - length(initTheta)))
          
          if ( !length( .iARcoeff ) ) 
            init.the[[spp.]] <- matrix(initTheta, nrow = n,
                                ncol = length(initTheta), byrow = TRUE)
          
          if (  !length( .idrift )  ) 
            init.mean[, spp. ] <- mean(y[, spp. ]) * (1 - sum(initTheta))
          
          justHelp <- GammaS[2:(nOrder[spp.] + 1)] * justHelp
          if (  !length( .ivar )  ) 
            init.var[, spp. ] <- max(0.01, GammaS[1] - (sum(justHelp)) )
          
          if (max(0.01, GammaS[1] - (sum(justHelp)) ) == 0.01 )
            warning("The noise variance is too close to zero.")
          
          if (  !length( .isd )  ) {
            init.sd[, spp. ] <- 
              max( sqrt(0.01), sqrt(init.var[, spp. ]))
            if (max( sqrt(0.01), sqrt(init.var[, spp. ]))  ==  sqrt(0.01))
              warning("The noise standard deviation",
                      " is too close to zero.")
          } 
        }   
        
        etastart <- cbind(if (.nodrift ) NULL else 
          theta2eta(init.mean, .ldrift , earg = .edrift),
          if ( .var.arg )
            theta2eta(init.var, .lvar , earg = .evar ) else 
              theta2eta(init.sd, .lsd , earg = .esd ))
        
        pre.eta <- numeric(0)
        for (kk in 1:NOS) 
          pre.eta <- cbind(pre.eta, theta2eta(cbind(init.the[[kk]]),
                                           .lARcoeff , earg = .eARcoeff ))
        
        pre.eta <- pre.eta[, interleave.VGAM(m.max * NOS, M1 = m.max,
                                             inverse = TRUE)]
      } else {
        pre.eta <- NULL
      }
      
      etastart <- cbind(etastart, pre.eta)
      etastart <- etastart[, interleave.VGAM(m1.max * NOS, M1 = m1.max)]
      pre.eta <- array(etastart, dim = c(n, m1.max, NOS))
      etastart <- numeric(0)
      
      for (ii in 1:NOS) 
        etastart <- cbind(etastart, pre.eta[, , ii][, 1:M1[ii]])
      
      etastart
    }), list(.ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
             .lARcoeff = lARcoeff,
             .edrift = edrift , .esd = esd , .evar = evar ,
             .eARcoeff = eARcoeff ,
             .idrift = idrift, .ivar = ivar , .isd = isd ,
             .iARcoeff = iARcoeff ,
             .nodrift  = nodrift , .var.arg = var.arg ,
             .type.likelihood = type.likelihood ,
             .Order = Order ))),
          
          
          
          
    linkinv = eval(substitute(function(eta, extra = NULL ) {
      n   <- nrow(eta)
      NOS <- extra$NOS
      nOrder <- extra$nOrder
      M1 <- 2 + nOrder - .nodrift
      
      lkInv <- break.VGAMextra(eta      = eta, 
                               M1       = M1, 
                               noInter  = .nodrift ,
                               bOrder   = nOrder,
                               NOS      = NOS,
                               lInter   = .ldrift ,
                               lvar     = .lvar ,
                               lsd      = .lsd ,
                               lcoeff1  = .lARcoeff ,
                               lcoeff2  = .lARcoeff ,
                               typeTS   = "AR",
                               Complete = TRUE,
                               varArg   = .var.arg )
      arp.drMean <- 
        if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else lkInv[[1]]
      arp.sd  <- lkInv[[2]]
      arp.var <- lkInv[[3]]
      arp.The <- lkInv[[4]]
      y.lags  <- array(0, dim = c(n, max(nOrder), NOS))
      
      for (ii in 1:NOS) 
        y.lags[,  1:nOrder[ii], ii] <- 
        WN.lags(y = cbind(extra$y[, ii]), lags = nOrder[ii]) 
      
      y.lags <- matrix(y.lags, nrow = n, ncol = max(nOrder) * NOS)
      fitted.ts <- array(arp.The * y.lags, dim = c(n, max(nOrder), NOS))
      fitted.ts <- arp.drMean +
        apply(fitted.ts, c(1, 3), function(x) sum(x))
      
      fitted.ts + extra$res
      
      
    }, list ( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
              .lARcoeff = lARcoeff ,
              .edrift = edrift , .evar = evar , .esd = esd ,
              .eARcoeff = eARcoeff ,
              .var.arg = var.arg , .Order = Order ,
              .nodrift = nodrift ))),
          
          
          
          
    last = eval(substitute(expression({
      n   <- nrow(eta)
      NOS <- extra$NOS
      nOrder <- extra$nOrder
      M1 <- 2 - .nodrift + nOrder
      M  <- sum(M1)
      
      sllast <- 
        break.VGAMextra(eta     = eta, 
                        M1      = M1, 
                        noInter = .nodrift ,
                        NOS     = NOS,
                        bOrder  = nOrder,
                        lInter  = .ldrift ,
                        lvar    = .lvar ,
                        lsd     = .lsd ,
                        lcoeff1 = .lARcoeff ,
                        lcoeff2 = .lARcoeff ,
                        typeTS  = "AR",
                        varArg  = .var.arg )
      arp.drMean <- if ( !( .nodrift ) ) sllast[[1]] else
        matrix(0.0, nrow = n, ncol = NOS)
      
      arp.sd  <- sllast[[2]]
      arp.var <- sllast[[3]]
      arp.The <- sllast[[4]]
      
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
          misc$link[ auxlast  + ii]    <- .lARcoeff
          misc$earg[[ auxlast  + ii]]  <- .eARcoeff
        }
        
        auxlast <- auxlast + nOrder[ll]
      }
      
      misc$var.arg  <- .var.arg
      misc$expected <- TRUE
      misc$process  <- "AR"
      misc$Order    <- nOrder
      misc$zero     <- .zero
      misc$NOS      <- NOS
      misc$M1       <- M1
      misc$nomean   <- .nodrift
      misc$flagAR   <- FALSE
      misc$eimARff  <- .eimARff
      misc$theta.names       <- parameters.names
      misc$control$wzepsilon <- control$wzepsilon
      misc$type.likelihood   <- .type.likelihood
      misc$multipleResponses <- TRUE
      fit$prior.weights <- w #20180315
      
      if ( !(.noChecks ) ) {
        
        saveNames <- names( fit$coefficients )
        Flag <- grep("(Intercept)", names(fit$coefficients) )
        
        if (length(Flag) == M) {
          fin.coe <- numeric(0)
          # At this point, No constraints have been established.
          myCoeffs <- coef(fit)[Flag]
          
          for (jj in 1:NOS) {
            pre.coe  <- myCoeffs[1:M1[jj]]
            myCoeffs <- myCoeffs[-(1:M1[jj])]
            if ( .nodrift ) {
              pre.coe <- eta2theta(pre.coe[-1], .lARcoeff , .eARcoeff )
            } else {
              pre.coe <- eta2theta(pre.coe[-(1:2)], .lARcoeff , .eARcoeff )
            }
            
            if (length(pre.coe) < max(nOrder))
              pre.coe <- c(pre.coe, rep(0, max(nOrder)- length(pre.coe)))
            
            fin.coe <- c(fin.coe, pre.coe)
          }
          
          ARroots <- checkTS.ffs(thetaEst = fin.coe,  # a vector.
                                 tsclass  = "AR",
                                 NofS     = NOS,
                                 pRoots   = FALSE,
                                 chOrder  = max(nOrder))
          
          flag <- (all( ARroots > 0 ) && any( ARroots <= 1 + 5e-3 ))
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
              .lARcoeff = lARcoeff ,
              .edrift = edrift , .esd = esd , .evar = evar , 
              .eARcoeff = eARcoeff ,
              .var.arg = var.arg , .nodrift = nodrift ,
              .type.likelihood = type.likelihood , .eimARff = eimARff ,
              .Order = Order , .zero = zero , .noChecks = noChecks ))),
          
          
          
          
    loglikelihood = eval(substitute( function(mu, y, w, residuals = FALSE,
                    eta, extra = NULL, summation = TRUE) {
      NOS <- extra$NOS
      y   <- extra$y
      nOrder <- extra$nOrder
      M1  <- if (length(extra$M1)) extra$M1 else 2 - .nodrift + nOrder
      n   <- nrow(eta)
      ordMax <- max(nOrder)[1]
      
      slloglik <- 
        break.VGAMextra(eta      = eta, 
                        M1       = M1, 
                        noInter  = .nodrift ,
                        NOS      = NOS,
                        bOrder   = nOrder,
                        lInter   = .ldrift ,
                        lvar     = .lvar ,
                        lsd      = .lsd ,
                        lcoeff1  = .lARcoeff ,
                        lcoeff2  = .lARcoeff ,
                        typeTS   = "AR",
                        Complete = TRUE,
                        varArg   = .var.arg )
      
      arp.drMean <- if ( .nodrift ) matrix(0.0, nrow = n,
                                           ncol = NOS) else slloglik[[1]]
      arp.sd  <- slloglik[[2]]
      arp.var <- slloglik[[3]]
      arp.The <- slloglik[[4]]
      
      y.lags  <- array(0, dim = c(n, max(nOrder), NOS))
      
      for (ii in 1:NOS) 
        y.lags[,  1:nOrder[ii], ii] <- 
        WN.lags(y = cbind(extra$y[, ii]), lags = nOrder[ii])
      
      y.lags  <- matrix(y.lags, nrow = n, ncol = max(nOrder) * NOS)
      mean.ts <- array(arp.The * y.lags, dim = c(n, max(nOrder), NOS))
      mean.ts <- arp.drMean + apply(mean.ts, c(1, 3), sum)
      names(mean.ts) <- NULL
      
      if (residuals) {
        stop("Loglikelihood not implemented yet", 
             " to handle residuals.")
      } else {
        loglik.terms <-
          c(w) * dARp(x = y, mean = mean.ts,
                      sd = arp.sd, log = TRUE)
        sum( loglik.terms )
      }
      
    }, list (  .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
               .lARcoeff = lARcoeff , 
               .edrift = edrift , .esd = esd , .evar = evar ,
               .eARcoeff = eARcoeff ,
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
                        noInter  = .nodrift ,
                        NOS      = NOS,
                        bOrder   = nOrder,
                        lInter   = .ldrift ,
                        lvar     = .lvar ,
                        lsd      = .lsd ,
                        lcoeff1  = .lARcoeff ,
                        lcoeff2  = .lARcoeff ,
                        typeTS   = "AR",
                        Complete = TRUE,
                        varArg   = .var.arg )
      arp.drMean <- if ( .nodrift ) matrix(0.0, nrow = n,
                                           ncol = NOS) else slloglik[[1]]
      arp.sd  <- slloglik[[2]]
      arp.var <- slloglik[[3]]
      arp.The <- slloglik[[4]]
      
      okay1 <- (all(is.finite(arp.sd)) && all(0 < arp.sd) &&
                  all(is.finite(arp.drMean)) ) 
      okay1
      
    }, list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
             .lARcoeff = lARcoeff , 
             .edrift = edrift , .esd = esd , .evar = evar ,
             .eARcoeff = eARcoeff ,
             .var.arg = var.arg , .Order = Order ,
             .type.likelihood = type.likelihood ,
             .nodrift = nodrift   ))),
          
          
    
    
          
   #vfamily = c("vgtsff", "vgltsff-class", "ARff"),
   vfamily  = c("ARXff", "ARMAvgltsmff"),
          
  
   
   
   
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
                       noInter  = .nodrift ,
                       NOS      = NOS,
                       bOrder   = nOrder,
                       lInter   = .ldrift ,
                       lvar     = .lvar ,
                       lsd      = .lsd ,
                       lcoeff1  = .lARcoeff ,
                       lcoeff2  = .lARcoeff ,
                       typeTS   = "AR",
                       Complete = TRUE, 
                       varArg   = .var.arg )
     arp.drMean <- 
       if ( .nodrift ) matrix(0.0, nrow = n, ncol = NOS) else 
         slderiv[[1]] 
     arp.sd  <- slderiv[[2]]
     arp.var <- slderiv[[3]]
     arp.The <- pre.The <- slderiv[[4]]
     pre.The <- array(pre.The, dim = c(n, ordMax, NOS))
     
     y.lags  <- array(0, dim = c(n, max(nOrder), NOS))
     for (ii in 1:NOS) 
       y.lags[,  1:nOrder[ii], ii] <- 
       WN.lags(y = cbind(extra$y[, ii]), lags = nOrder[ii]) 
     
     y.lags <- matrix(y.lags, nrow = n, ncol = max(nOrder) * NOS)
     y.means <- array(arp.The * y.lags, dim = c(n, max(nOrder), NOS))
     y.means <- y - (arp.drMean + 
                       apply(y.means, c(1, 3), function(x) sum(x)))
     
     dl.drfmean <- y.means / arp.var
     if ( .var.arg ) {
       dl.dvar <- y.means^2 / (2 * arp.var^2) - 1 / (2 * arp.var)
     } else {
       dl.dsd <- y.means^2 / arp.sd^3 - 1 / arp.sd
     }
     
     y.lags <- array(y.lags, dim = c(n, max(nOrder), NOS))
     dl.dThe <- array(NA_real_, dim = c(n, max(nOrder), NOS))
     for (ii in 1:NOS)
       dl.dThe[, , ii] <- (array( y.means, dim = c(n, 1, NOS))[, , ii] *
                             y.lags[, , ii]) / arp.var[, ii]
     
     dl.dThe <- interleaveArray.VGAMextra(dl.dThe)
     
     # Partial derivatives of theta wrt eta #
     drfm.deta <- dtheta.deta(arp.drMean, .ldrift , earg = .edrift )
     if ( .var.arg ) {
       dvar.deta <- dtheta.deta(arp.var , .lvar , earg = .evar )
       dVarSD    <- dvar.deta
     } else {
       dsd.deta  <- dtheta.deta(arp.sd  , .lsd  ,  earg = .esd )
       dVarSD    <- dsd.deta
     }
     
     dThedEta <- interleaveArray.VGAMextra(array(arp.The,
                                                 dim = c(n, ordMax, NOS)))
     dThedEta <- dtheta.deta(dThedEta, .lARcoeff , earg = .eARcoeff)
     
     # The chain rule to obtain the total derivative of 'l' wrt
     # eta: (dl/deta_i) = (dl / dtheta_i) * (dtheta_i / deta_i)
     dl.deta1 <- c(w) * dl.drfmean * drfm.deta
     if ( .var.arg ){
       dl.deta2 <- c(w) * dl.dvar * dvar.deta
     } else {
       dl.deta2 <- c(w) * dl.dsd * dsd.deta
     }
     
     a.help  <- interleaveArray.VGAMextra(c(w) * dl.dThe * dThedEta,
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
             .lARcoeff  = lARcoeff , 
             .edrift = edrift , .evar = evar , .esd = esd ,
             .eARcoeff = eARcoeff , 
             .var.arg = var.arg , .Order = Order , 
             .nodrift = nodrift ))),
          
          
          
          
    weight = eval(substitute(expression({
      dThedEta <- interleaveArray.VGAMextra(dThedEta, inverse = TRUE)
      
      ### Exact EIMs.
      if ( .eimARff ) {
        exact.eim <- vector("list", NOS); col.count <- 0
        for (jj in 1:NOS) {
          exact.eim[[jj]] <- 
                  ARpEIM.G2(y  = extra$y[, jj],
                            drift   = arp.drMean[, jj, drop = FALSE],
                            sdError = arp.sd[, jj, drop = FALSE],
                            ARpcoeff = pre.The[, 1:nOrder[jj], jj],
                            var.arg = .var.arg ,
                            order   = nOrder[jj],
                            nodrift = .nodrift )
          
          comb.wz <- combVGAMextra(1:M1[jj], nodrift = FALSE)
          a.bind  <- cbind(if ( .nodrift ) NULL else drfm.deta[, jj],
                           dVarSD[, jj],
                           dThedEta[, 1:nOrder[jj], jj])
          
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
        ned2l.dsmn   <- (1 / arp.var) 
        ned2l.dvarSD <- if ( .var.arg ) 1/(2 * arp.var^2) else 2 / arp.var
        wz <- cbind( if ( .nodrift ) NULL else 
          ned2l.dsmn * drfm.deta^2, ned2l.dvarSD * dVarSD^2)
        
        dThe.bis <- matrix(dThedEta, nrow = n, ncol = ordMax * NOS)
        
        gammas.y <- apply(y, 2, function(x) {
          cross.gammas(x = x, lags = 0)
        })
        
        comp.ders <- matrix(c(gammas.y), nrow = n, ncol = NOS * ordMax,
                            byrow = TRUE) / matrix(arp.var, nrow = n,
                                                   ncol = NOS * ordMax)
        comp.ders <- comp.ders *
          matrix(dThe.bis^2, nrow = n, ncol = ordMax * NOS)
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
              .nodrift  = nodrift , .eimARff = eimARff ,
              .idrift = idrift ,
              .ivar = ivar , .isd = isd , 
              .iARcoeff = iARcoeff ))) ) # End of new("vgltsmff")
    alo
  }  # End of ARXff












### Family function ARMAX.GARCHff
ARMAX.GARCHff <- 
  function(ARMAorder   = c(1, 0),
           GARCHorder  = c(1, 0),
           type.TS     = c("ARCH", "GARCH", "IGARCH",
                           "Taylor-Schwert", "A-GARCH",
                           "Log-GARCH", "M-GARCH")[1],
           type.param = c("residuals", "observed")[1],
           Cov.on.Var = FALSE,
           noChecks   = FALSE,
           G1.transform = NULL,
           G2.transform = NULL,
           lvar       = NULL,
           lsd        = NULL,
           ldrift     = "identitylink", 
           lARcoeff   = "identitylink",
           lMAcoeff   = "identitylink",
           idrift     = NULL,
           iARcoeff   = NULL,
           iMAcoeff   = NULL) {
    
    Include.cov <- c("ARMAmodel", "GARCHmodel")[1]
    if (length(ARMAorder) != 2)
      stop("'ARMAorder is a length-2 vector.'")
    ARorder <- ARMAorder[1]
    MAorder <- ARMAorder[2]; rm(ARMAorder)
    
    if (length(G1.transform) && !is.function(G1.transform))
      stop("Wrong input for argument 'G1.transform'.
           It must be a function.",
           " Enter NULL for the identity function.")
    
    if (length(G2.transform) && !is.function(G2.transform))
      stop("Wrong input for argument 'G1.transform'.
           It must be a function.",
           " Enter NULL for the identity function.")
    
    if ( !Is.Numeric(x = ARorder, isInteger = TRUE) || ARorder < 0)
      stop("Wrong input for argument 'ARorder'.")
    
    if ( !Is.Numeric(x = MAorder, isInteger = TRUE) || MAorder < 0)
      stop("Bad input for argument 'MAorder'.")
    
    GarchOrd <- GARCHorder; rm(GARCHorder)
    if ( !Is.Numeric(x = GarchOrd, length.arg = 2))
      stop("Wrong input for argument 'GarchOrd'.")
    
    ###  FURTHER ARGUMENTS from ARXff()
    nodrift  <- FALSE
    flagAR   <- (ARorder == 0)
    flagMA   <- (MAorder == 0)
    type.EIM <- "exact"
    noChecks <- noChecks
    var.arg  <- NULL
    
    g1f  <- G1.transform; rm(G1.transform)
    g2f  <- G2.transform; rm(G2.transform)
    ord1 <- GarchOrd[1]
    ord2 <- GarchOrd[2]
    type.TS  <- match.arg(type.TS, c("ARCH", "GARCH", "IGARCH",
                                     "Taylor-Schwert", "A-GARCH",
                                     "Log-GARCH", "M-GARCH",
                                     "EGARCH"))[1]
    
    if ((type.TS == "ARCH") && ord2) {
      ord2 <- 0
      warning("Arguments 'GARCHord' and 'type.TS' do not match. ",
              " Set 'type.TS = GARCH' to fit GARCH models." )
    }
    
    if ((type.TS == "ARCH") && !ord1) {
      stop("ARCH order required" )
    }

    fin.blurb <- (1 + ord1 + ord2) + (1 + ARorder + MAorder)
    Order   <- c(ifelse(flagAR, 1 , ARorder))
    OrderMA <- c(ifelse(flagMA, 1 , MAorder))
    type.param <- match.arg(type.param, c("residuals", "observed"))[1]
    
    if (!length(var.arg))
      var.arg <- switch(type.TS ,
                        "ARCH"    = TRUE,
                        "GARCH"   = TRUE,
                        "IGARCH"  = TRUE,
                        "Taylor-Schwert" = FALSE,
                        "A-GARCH" = TRUE,
                        "Log-GARCH" = TRUE,
                        "M-GARCH" = TRUE,
                        "EGARCH"  = TRUE)
    
    if (!length(lvar)) 
      lvar <- switch(type.TS ,
                     "ARCH"   = "identitylink",
                     "GARCH"  = "identitylink",
                     "IGARCH" = "identitylink",
                     "Taylor-Schwert" = "identitylink",
                     "A-GARCH"   = "identitylink",
                     "Log-GARCH" = "loge",
                     "M-GARCH" = "loge",
                     "EGARCH"  = "loge")
      #lvar <-  if (type.TS == "GARCH") "identitylink" else "loge"
    
    
    if (!length(lsd))
      lsd  <- switch(type.TS ,
                     "ARCH"   = "identitylink",
                     "GARCH"  = "identitylink",
                     "IGARCH" = "identitylink",
                     "Taylor-Schwert" = "identitylink",
                     "A-GARCH"   = "identitylink",
                     "Log-GARCH" = "loge",
                     "M-GARCH" = "loge",
                     "EGARCH"  = "loge")
    
    ivar <- NULL
    isd  <- NULL
    
    temp.drift <- if (nodrift) NULL else
      if (flagAR) "Mean" else
        if (flagMA) "drift" else "drift.mean"
    realM1  <- 2 - nodrift + ARorder + MAorder
    
    if (FALSE) {
      I.cov <- match.arg(Include.cov, c("ARMAmodel", "GARCHmodel"))[1]
      temp.zero <- if (I.cov == "GARCHmodel") temp.drift else
        if (I.cov == "ARMAmodel") "Var" else NULL  
    }
    
    
    zero <- c(if (var.arg) "Var" else "sd",
              if ( flagAR ) NULL else "ARcoeff",
              if ( flagMA ) NULL else "MAcoeff")
    
    # Not used yet.
    real.zero <- c(temp.drift,
                  if ( flagAR ) NULL else "ARcoeff",
                  if ( flagMA ) NULL else "MAcoeff")
    
    
    if (flagAR && length(iARcoeff)) {
      iARcoeff <- 2e-1
      warning("Initial values for the AR mean-model component are ignored")
    }
    
    if (flagMA && length(iMAcoeff)) {
      iMAcoeff <- 2e-1
      warning("Initial values for the MA mean-model component are ignored")
    }
    
    if ( length(idrift) && !Is.Numeric(idrift) )
      stop("Bad entry for argument 'idrift'.") 
    
    if ( length(isd) && !Is.Numeric(isd, Nnegative = TRUE) )
      stop("Wrong input for argument 'isd'.") 
    
    if ( length(ivar) && !Is.Numeric(ivar, Nnegative = TRUE) )
      stop("Wrong input for argument 'ivar'.") 
    
    if ( !is.logical(var.arg) ||  !is.logical(nodrift) )
      stop("Arguments 'var.arg' and 'nodrift' must \n",
           "be logical.")
    
    if (!is.logical(noChecks))
      stop("Wrong input for argument 'noChecks'")
    
    type.likelihood <- "exact"
    type.EIM  <- match.arg(type.EIM, 
                           c("exact", "approximate"))[1]
    
    eimARff  <- (type.EIM == "exact")
    
    ldrift <- match.arg(ldrift, "identitylink")
    ldrift <- as.list(substitute(ldrift))
    edrift <- link2list(ldrift)
    ldrift <- attr(edrift, "function.name")
    
    lsd <- as.list(substitute(lsd))
    esd <- link2list(lsd)
    lsd <- attr(esd, "function.name")
    
    lvar <- as.list(substitute(lvar))
    evar <- link2list(lvar)
    lvar <- attr(evar, "function.name")
    
    lARcoeff <- as.list(substitute(lARcoeff))
    eARcoeff <- link2list(lARcoeff)
    lARcoeff <- attr(eARcoeff,"function.name")
    
    lMAcoeff <- as.list(substitute(lMAcoeff))
    eMAcoeff <- link2list(lMAcoeff)
    lMAcoeff <- attr(eMAcoeff,"function.name")
    
    
    blAr <- max(Order)
    blMa <- max(OrderMA)
    pre.blurb  <- paste("ARcoeff", 1:blAr, sep ="")
    pre.blurb3 <- paste("MAcoeff", 1:blMa, sep ="")
    blurb.vec <- 
      c(namesof("ARdrift", ldrift, earg = edrift, tag = FALSE),
        if (var.arg)
          namesof("noiseVar", lvar, earg = evar, tag = FALSE)
        else
          namesof("noiseSD", lsd, earg = esd, tag = FALSE),
        namesof(pre.blurb, lARcoeff , earg = eARcoeff , tag = FALSE),
        namesof(pre.blurb3, lMAcoeff , earg = eMAcoeff , tag = FALSE))
    
    blurb.vec2 <- NULL
    for (ii in 3:(blAr + 2)) 
      #if (ii == (blAr + 2))
      #  blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii])) else
      blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii], ", "))
    
    for (ii in (3 + blAr):(blAr + + blMa + 2)) 
      if (ii == (blAr + blMa + 2))
        blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii])) else
          blurb.vec2 <- c(blurb.vec2, c(blurb.vec[ii], ", "))
    
    alo <- 
      new("vgltsmff", 
          blurb = c(c("VGLTSMs: ", fin.blurb, 
                      "-parameter ARMAX model",
                  " of order-(", ARorder, ",", MAorder  , ") with GARCH(",
                      ord1, ", ", ord2, ") errors",
                      "\n\n",
                      "Links:    "), 
                    blurb.vec[1], ", ",
                    blurb.vec[2], ", ",
                    blurb.vec2, ". \n",
              c("Model:    Y_t = ",
                if (nodrift) NULL else "drift +",
                " (ARcoeff1)*Y_{t - 1} + ... + (ARcoeff",blAr,")*Y_{t - ",
                      blAr,"} + ", "(MAcoeff1)*e_{t - 1} + ... 
                      + (MAcoeff",blMa,")*e_{t - ",
                      blMa, "} + e_{t} ,", "\n",
                      "where e_{t} ~ N(0, noiseSD^2_t), and ", "\n"),
                      #"E[Y_t]   = ARdrift/ (1 - Sum(ARcoeff))"),
                    "\n",
                    "  g(noiseSD^2_t) = w + A1 * e^2_{t - 1} + ...",
                    " + Ap * e^2_{t - r} + ", "\n", 
                    "                 B1 * noiseSD^2_{t - 1} + ...",
                    " + Bq * noiseSD^2_{t - s}"),
          
          
          
          
    constraints = eval(substitute(expression({
      M1vec <-  2 + nOrder - .nodrift
      
      zero2 <- if (extra$NOxx) .real.zero else .zero
      
      constraints <- cm.zero.VGAM(constraints, x = x, zero = .zero ,
                                  M = M,
                                  predictors.names = parameters.names,
                                  #predictors.names = predictors.names,
                                  M1 = M1vec)
      
      col.x <- c(grep("ARCH", colnames(x)),
                 grep("G1", colnames(x)),
                 grep("G2", colnames(x)))
      constra.temp <- colnames(x)[col.x]
      for (ii in constra.temp)
        constraints[[ii]] <- cbind(c(0, 1, rep(0, .ARord + .MAord )))
      
      constra.temp <- colnames(x)[-c(1, col.x)]
      for (ii in constra.temp)
        constraints[[ii]] <- cbind(c(1, ( .Cov.on.Var ),
                                     rep(0, .ARord + .MAord ) ))
      
      
    }), list( .zero = zero, .Order = Order, .nodrift = nodrift ,
              .temp.drift = temp.drift , .real.zero = real.zero ,
              .ARord = ARorder , .MAord = MAorder ,
              .Cov.on.Var = Cov.on.Var ))),
          
          
          
          
          
    first = eval(substitute(expression({
      intercept.only <- (ncol(x) == 1 && colnames(x) == "(Intercept)")
      
      if (!intercept.only && FALSE)
        stop("Currently, this family function only handles ",
             "intercept-only models.")
      
      
      g1f <- .g1f
      g2f <- .g2f
      extra$NOxx <- intercept.only
      
      iniOrd <- max(10, .ord2 + 1)
      nn <- NROW(y)
      
      if (NCOL(y) > 1)
        stop("Currently, only one-response handled.")
      
      if ( !( .flagMA ) ) {
        lm.help <- lsfit(x = WN.lags(y = cbind(y),
                                     lags = max(8, .ARorder + .MAorder )),
                         y = y, intercept = TRUE)
      } else {
        if ( .flagAR ) {
          lm.help   <- lm(y ~ 1, data = data.frame(y = y))
        } else {
          lm.help <- lsfit(x = WN.lags(y = cbind(y), lags = .ord1 ),
                           y = y, intercept = TRUE)
        }
      }
      
      extra$res <- saved.res <- cbind(residuals(lm.help))
      
      x.matrix  <- x
      x1.mat    <- x2.mat <- NULL
      
      
      
      if ( .type.TS == "ARCH") {
        
        mytype <- if ( .type.param == "observed" ) y^2 else saved.res^2
        x1.mat <- WN.lags(y = if (!length(g1f)) cbind(mytype) else
                                             g1f(cbind(sqrt(mytype))),
                          lags = .ord1 ); rm(mytype)
        
        colnames(x1.mat) <- 
          if (!length(g1f)) paste("ARCH(", 1:( .ord1 ), ")", sep = "") else
            if ( .type.param == "observed" )
             paste("G1(ys)-(", 1:( .ord1 ), ")", sep = "") else 
               paste("G1(errors)-(", 1:( .ord1 ), ")", sep = "")
        
        x.matrix   <- cbind(x.matrix, x1.mat, x2.mat)
        list.names <- vector("list", NCOL(x) + ( .ord1 ) + ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        
        for (ii in 1:(NCOL(x) + ( .ord1 ) + ( .ord2 )  ))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
      }
      
      
      if ( .type.TS == "GARCH" ) {
        
        if ( .ord1 ) {
          
          mytype <- if ( .type.param == "observed" ) y^2 else saved.res^2
          x1.mat <- WN.lags(y = if (!length(g1f)) cbind(mytype) else
                     g1f(cbind(sqrt(mytype))), lags = .ord1 )
          rm(mytype)
          
          colnames(x1.mat) <- 
          if (!length(g1f)) paste("ARCH(", 1:( .ord1 ), ")", sep = "") else
              if ( .type.param == "observed" )
                paste("G1(ys)-(", 1:( .ord1 ), ")", sep = "") else 
                  paste("G1(errors)-(", 1:( .ord1 ), ")", sep = "")
          
          #mytype <- if ( .type.param == "observed" ) y^2 else saved.res^2
          #x1.mat <- WN.lags(y = cbind(mytype),
          #                  lags = .ord1 ); rm(mytype)
          #colnames(x1.mat) <- paste("ARCH(", 1:( .ord1 ), ")", sep = "")
        }
        
        ## Estimate sigma_t, if required. Here, e_t ~ N(0, sigma_t^2)
        ## Then normalsdff() is utilized, link = "identitylink"
        if ( .ord2 ) {
          temp.data <- data.frame(y = saved.res,
                                  WN.lags(data.frame(saved.res^2),
                                          lags = iniOrd))
          names(temp.data) <- c("y", paste("x", 2:(iniOrd + 1), sep = ""))
          
          myform <- character(0)
          for (ii in 2:(iniOrd + 1)) {
            pre <- paste(paste("x", ii, sep = ""),
                         if (ii == iniOrd + 1) "" else " + ", sep ="")
            myform <- paste(myform, pre, sep = "")
          }
          
          myform <- paste("y ~ 1 +", myform, sep = " ")
          ersfit <- vglm(as.formula(myform),
                         family = normal1sdff(fixed.mean = 0, zero = NULL,
                                   link = "identitylink", var.arg = TRUE),
                         trace = FALSE, smart = FALSE,
                         data = temp.data, eps = 5e-5)
          coess <- c(coef(ersfit, matrix= TRUE))
          #Sigmas <- cbind(rep(1, nn),
          #                WN.lags(y = cbind(saved.res^2),
          #                         lags = iniOrd)) %*% coess
          
          Sigmas <- fitted.values(ersfit)
          
          x2.mat <- WN.lags(y = if (!length(g2f)) cbind(Sigmas) else
                               g2f(cbind(sqrt(Sigmas))), lags = .ord2)
          
          colnames(x2.mat) <-  if (!length(g2f)) paste("GARCH(", 
                                       1:( .ord2 ), ")", sep = "") else
                  paste("G2(Sigma)-(", 1:( .ord2 ), ")", sep = "")
          
          #x2.mat <- WN.lags(y = cbind(Sigmas), lags = .ord2 )
          #colnames(x2.mat) <- paste("GARCH(", 
          #                          1:( .ord2 ), ")", sep = "")
        }
        
        
        x.matrix   <- cbind(x.matrix, x1.mat, x2.mat)
        list.names <- vector("list", NCOL(x) + ( .ord1 ) + ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        for (ii in 1:(NCOL(x) + ( .ord1 ) + ( .ord2 )  ))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
        
      }
      
      
      
      
      if ( .type.TS == "IGARCH") {
        
        if ( .ord1 ) {
          x1.mat <- WN.lags(y = cbind(saved.res^2), lags = .ord1 )
          colnames(x1.mat) <- 
            paste("IGARCH-e^2(", 1:( .ord1 ), ")", sep = "")
        } 
        
        ## Estimate sigma_t^2 sum(a) + sum(b) = 1, e_t ~ N(0, sigma_t^2)
        ## Then normalsdff(), link = "identitylink"
        if ( .ord2 ) {
          
          iniOrd   <- .ord1
          res2.mat <- WN.lags(y = cbind(saved.res^2), lags = iniOrd + 1)
          x2.mat   <- matrix(NA_real_, nn, iniOrd)
          for (ii in 1:iniOrd) 
            x2.mat[, ii] <- res2.mat[, ii + 1] - res2.mat[, 1]
          
          temp.data <- data.frame(y = saved.res, offs = res2.mat[, 1],
                                  cbind(x2.mat))
          names(temp.data) <- c("y", "offset",
                                paste("x", 2:(iniOrd + 1), sep = ""))
          
          myform <- character(0)
          for (ii in 2:(iniOrd + 1)) {
            pre <- paste(paste("x", ii, sep = ""),
                         if (ii == iniOrd + 1) "" else " + ", sep ="")
            myform <- paste(myform, pre, sep = "")
          }
          
          myform <- paste("y ~ 1 + offset(offset) + ", myform)
          ersfit <- vglm(as.formula(myform),
                         family = normal1sdff(fixed.mean = 0, zero = NULL,
                                   link = "identitylink", var.arg = TRUE),
                         trace = FALSE, smart = FALSE,
                         data = temp.data, eps = 5e-5)
          #  coess <- c(coef(ersfit, matrix= TRUE))
          #Sigmas <- cbind(rep(1, nn),
          #                WN.lags(y = cbind(saved.res^2),
          #                        lags = iniOrd)) %*% coess + res2.mat[, 1]
          Sigmas <- fitted.values(ersfit)
          
          x2.mat <- WN.lags(y = cbind(Sigmas), lags = .ord2 ) 
          colnames(x2.mat) <- paste("IGARCH-s^2(", 
                                    1:( .ord2 ), ")", sep = "")
        }
        
        x.matrix <- cbind(x.matrix, x1.mat, x2.mat)
        list.names <- vector("list", NCOL(x) + ( .ord1 ) + ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        
        for (ii in 1:(NCOL(x) + ( .ord1 ) + ( .ord2 )))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
      }
      
      
      
      
      if (.type.TS == "Taylor-Schwert") {
        
        if ( .ord1 ) {
          x1.mat <- WN.lags(y = cbind(abs(saved.res)), lags = .ord1 )
          colnames(x1.mat) <- paste("TS-ARCH(", 1:( .ord1 ), ")",sep = "")
        }
        
        ## Estimate sigma_t, if required. Here, e_t ~ N(0, sigma_t^2)
        ## Then normalsdff() is utilized, link = "identitylink"
        if ( .ord2 ) {
          temp.data <- data.frame(y = saved.res,
                            WN.lags(cbind(abs(saved.res)), lags = iniOrd))
          names(temp.data) <- c("y", paste("x", 2:(iniOrd + 1), sep = ""))
          myform <- character(0)
          for (ii in 2:(iniOrd + 1)) {
            pre <- paste(paste("x", ii, sep = ""),
                         if (ii == iniOrd + 1) "" else " + ", sep ="")
            myform <- paste(myform, pre, sep = "")
          }
          
          myform <- paste("y ~ 1 +", myform, sep = " ")
          ersfit <- vglm(as.formula(myform),
                         family = normal1sdff(fixed.mean = 0, zero = NULL,
                                  link = "identitylink", var.arg = FALSE),
                         trace = FALSE, smart = FALSE,
                         data = temp.data, eps = 5e-5)
          rm(temp.data)
          Sigmas <- fitted.values(ersfit)
          #coess <- c(coef(ersfit, matrix= TRUE))
          #Sigmas <- cbind(rep(1, nn),
          #                WN.lags(y = cbind(abs(saved.res)),
          #                        lags = iniOrd)) %*% coess
          
          x2.mat <- WN.lags(y = cbind(Sigmas), lags = .ord2 )
          colnames(x2.mat) <- paste("TS-GARCH(", 
                                    1:( .ord2 ), ")", sep = "")
        }
        
        x.matrix   <- cbind(x.matrix, x1.mat, x2.mat)
        list.names <- vector("list", NCOL(x) + ( .ord1 ) + ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        
        for (ii in 1:(NCOL(x) + ( .ord1 ) + ( .ord2 )  ))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
      }
      
      
      
      
      if ( .type.TS == "A-GARCH") {
        
        if ( .ord1 ) {
          x1.mat <- cbind(WN.lags(y = cbind(saved.res^2), lags = .ord1 ),
                          WN.lags(y = cbind(saved.res),   lags = .ord1 ))
          colnames(x1.mat) <- c(paste("A-GARCH-e^2", 1:( .ord1 ), sep = ""),
                                paste("Lev-AGARCH", 1:( .ord1 ), sep = ""))
        }
        
        
        if ( .ord2 ) {
          iniOrd <- 2
          temp.data <- data.frame(y = saved.res,
                     x1 = WN.lags(y = cbind(saved.res^2), lags = iniOrd ),
                     y1 = WN.lags(y = cbind(saved.res), lags = iniOrd ))
          names(temp.data) <- c("y", paste("x2", 2:(iniOrd + 1), sep = ""),
                                paste("e", 2:(iniOrd + 1), sep = ""))
          
          myform <- character(0)
          for (ii in 2:length(names(temp.data))) {
            pre <- paste(paste(names(temp.data)[ii], "", sep = ""),
              if (ii ==length(names(temp.data))) "" else " + ", sep ="")
            myform <- paste(myform, pre, sep = "")
          }
          
          myform <- paste("y ~ 1 +", myform, sep = " ")
          ersfit <- vglm(as.formula(myform),
                         family = normal1sdff(fixed.mean = 0, zero = NULL,
                                   link = "identitylink", var.arg = TRUE),
                         trace = FALSE, smart = FALSE,
                         data = temp.data, eps = 5e-5)
          
          Sigmas <- fitted.values(ersfit)
          #coess <- c(coef(ersfit, matrix= TRUE))
          #Sigmas <- cbind(rep(1, nn),
          #                as.matrix(temp.data[, -1])) %*% coess
          
          x2.mat <- WN.lags( y = cbind(Sigmas), lags = .ord2 )
          colnames(x2.mat) <- paste("A-GARCH-s^2", 1:( .ord2 ), sep = "")
        }
        
        x.matrix   <- cbind(x.matrix, x1.mat, x2.mat)
        list.names <- vector("list", NCOL(x) + 2 * ( .ord1 ) + ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        
        for (ii in 1:(NCOL(x) + 2 * ( .ord1 ) + ( .ord2 )  ))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
      }
      
      
      
      
      if ( .type.TS == "Log-GARCH" ) {
        
        if ( .ord1 ) {
          
          mytype <- if ( .type.param == "observed" ) y else saved.res
          x1.mat <- WN.lags(y = cbind(abs(mytype)), 
                            lags = .ord1 ); rm(mytype)
          #x1.mat <- WN.lags(y = cbind(abs(saved.res)), lags = .ord1 )
          colnames(x1.mat) <- paste("Log-ARCH(", 1:( .ord1 ), ")",sep = "")
        }
        
        
        ## Estimate sigma_t, if required. Here, e_t ~ N(0, sigma_t^2)
        ## Then normalsdff() is utilized, link = "loge", since
        ## log(sigma_t) is modeled.
        if ( .ord2 ) {
          temp.data <- data.frame(y = saved.res,
                                  WN.lags(data.frame(abs(saved.res)),
                                          lags = iniOrd))
          names(temp.data) <- c("y", paste("x", 2:(iniOrd + 1), sep = ""))
          
          myform <- character(0)
          for (ii in 2:(iniOrd + 1)) {
            pre <- paste(paste("x", ii, sep = ""),
                         if (ii == iniOrd + 1) "" else " + ", sep ="")
            myform <- paste(myform, pre, sep = "")
          }
          
          myform <- paste("y ~ 1 +", myform, sep = " ")
          ersfit <- vglm(as.formula(myform),
                         family = normal1sdff(fixed.mean = 0, zero = NULL,
                                          link = "loge", var.arg = FALSE),
                         trace = FALSE, smart = FALSE,
                         data = temp.data, eps = 5e-5)
          rm(temp.data)
          Sigmas <- fitted.values(ersfit)
          #coess <- c(coef(ersfit, matrix= TRUE))
          #Sigmas <- exp(cbind(rep(1, nn),
          #                    WN.lags(y = cbind(abs(saved.res)),
          #                            lags = iniOrd)) %*% coess)
          
          x2.mat <- WN.lags(y = cbind(log(Sigmas)), lags = .ord2 )
          colnames(x2.mat) <- paste("Log-GARCH(", 
                                    1:( .ord2 ), ")", sep = "")
        }
        
        x.matrix   <- cbind(x.matrix, x1.mat, x2.mat)
        list.names <- vector("list", NCOL(x) + ( .ord1 ) + ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        
        for (ii in 1:(NCOL(x) + ( .ord1 ) + ( .ord2 )  ))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
      }
      
      
      
      
      
      if ( .type.TS == "M-GARCH" ) {
        
        if ( .ord1 ) {
          mytype <- if ( .type.param == "observed" ) y^2 else saved.res^2
          x1.mat <- WN.lags(y = cbind(log(mytype)),
                            lags = .ord1 ); rm(mytype)
          #x1.mat <- WN.lags(y = cbind(log(saved.res^2)), lags = .ord1 )
          colnames(x1.mat) <- paste("M-ARCH(", 1:( .ord1 ), ")",
                                    sep = "")
        }
        
        
        ## Estimate sigma_t, if required. Here, e_t ~ N(0, sigma_t^2)
        ## Then normalsdff() is utilized, link = "loge", since
        ## log(sigma_t^2) is modeled.
        if ( .ord2 ) {
          temp.data <- data.frame(y = saved.res,
                                  WN.lags(data.frame(log(saved.res^2)),
                                          lags = iniOrd))
          names(temp.data) <- c("y", paste("x", 2:(iniOrd + 1), sep = ""))
          
          myform <- character(0)
          for (ii in 2:(iniOrd + 1)) {
            pre <- paste(paste("x", ii, sep = ""),
                         if (ii == iniOrd + 1) "" else " + ", sep ="")
            myform <- paste(myform, pre, sep = "")
          }
          
          myform <- paste("y ~ 1 +", myform, sep = " ")
          ersfit <- vglm(as.formula(myform),
                         family = normal1sdff(fixed.mean = 0, zero = NULL,
                                           link = "loge", var.arg = TRUE),
                         trace = FALSE, smart = FALSE,
                         data = temp.data, eps = 5e-5)
          rm(temp.data)
          Sigmas <- fitted.values(ersfit)
          #coess <- c(coef(ersfit, matrix= TRUE))
          #Sigmas <- exp(cbind(rep(1, nn),
          #                    WN.lags(y = cbind(log(saved.res^2)),
          #                            lags = iniOrd)) %*% coess)
          
          x2.mat <- WN.lags(y = cbind(log(Sigmas)), lags = .ord2 )
          colnames(x2.mat) <- paste("M-GARCH(", 
                                    1:( .ord2 ), ")", sep = "")
        }
        
        x.matrix   <- cbind(x.matrix, x1.mat, x2.mat)
        list.names <- vector("list", NCOL(x) + ( .ord1 ) + ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        
        for (ii in 1:(NCOL(x) + ( .ord1 ) + ( .ord2 )  ))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
      }
      
      
      
      if ( .type.TS ==  "EGARCH" ) {
        
        gams <- numeric(nn)
        gams[1] <- var(saved.res)
        for (ii in 2:(nn - 1)) {
          gams[ii] <- var(saved.res[1:(nn - ii + 1)])
        }
        gams[nn] <- var(saved.res)
        
        if (FALSE) {
          if ( .ord1 ) 
            x1.mat <- cbind(WN.lags(y = cbind(log(gams)),
                                    lags = .ord1 ))
          if ( .ord2 ) {
            x2.mat <- cbind(WN.lags(y = abs(saved.res)/sqrt(gams) -
                                      sqrt(2/pi), lags = .ord2 ),
                         WN.lags(y = saved.res/sqrt(gams), lags = .ord2 ))
          }
          x.mat <- cbind(x1.mat, x2.mat)
          colnames(x.mat) <- c(if ( .ord1 ) paste("Log EGARCH(", 
                                    1:( .ord1 ), ")", sep = "") else NULL,
                             if ( .ord2 ) paste("abs EGARCH(", 1:( .ord2 ),
                                                 ")",  sep = "") else NULL,
                               if ( .ord2 ) paste("leverage(", 1:( .ord2 ),
                                                  ")", sep = "") else NULL)
          coess <- lsfit(x = x.mat[-1, ], y = log(gams[-1]),
                         intercept = TRUE)$coefficients
          Sigmas <- exp(cbind(rep(1, nn), x.mat) %*% coess)
        }
        
        Sigmas <- extra$Sigmas <- gams
        x1.mat <- x2.mat <- NULL
        
        if ( .ord1 ) {
          x1.mat <- WN.lags(y = cbind(log(Sigmas)), lags = .ord1 ,
                         to.complete = rep(mean(log(Sigmas)) * 0, .ord1 ))
        }
        
        if (  ( .ord2 ) && ( .ord1 )  ) {
          x2.mat <- cbind(WN.lags(y = abs(saved.res)/sqrt(Sigmas) -
                                    sqrt(2/pi), lags = .ord2 ),
                      WN.lags(y = saved.res/sqrt(Sigmas), lags = .ord2 ))
        }
        
        x.matrix <- cbind(x.matrix, x1.mat, x2.mat)
        colnames(x.matrix) <- c(colnames(x), 
                        c(if ( .ord1 ) paste("Log EGARCH(", 
                                   1:( .ord1 ), ")", sep = "") else NULL,
                          if ( .ord2 ) paste("abs EGARCH(", 1:( .ord2 ),
                                            ")",  sep = "") else NULL,
                          if ( .ord2 ) paste("leverage(", 1:( .ord2 ),
                                             ")", sep = "") else NULL))
        
        list.names <- vector("list", NCOL(x) + ( .ord1 ) + 2 * ( .ord2 ))
        names(list.names) <- colnames(x.matrix)
        
        for (ii in 1:(NCOL(x) + ( .ord1 ) + 2 * ( .ord2 )  ))
          list.names[[ii]] <- ii
        attr(x.matrix, "assign") <- list.names
        
      }
      
      x <- x.matrix
      
    }), list( .ord1 = ord1 , .ord2 = ord2, 
              .type.param = type.param , .type.TS = type.TS ,
              .flagAR = flagAR , .flagMA = flagMA ,
              .ARorder = ARorder, .MAorder = MAorder ,
              .g1f = g1f , .g2f = g2f ))),
          
          
          
          
          
    infos = eval(substitute(function(...) {
      
      ## M1 was numeric(NA). Changed on 2018/01/19 to handle the
      ## Hauck - Donner effects. VGAM 1.0-5.
      list(M1       = .realM1 ,
           eimARff  = .eimARff ,
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
           nodrift  = .nodrift ,
           lARcoeff = .lARcoeff ,
           lMAcoeff = .lMAcoeff ,
           eARcoeff = .eARcoeff ,
           eMAcoeff = .eMAcoeff ,
           type.likelihood = .type.likelihood ,
           type.EIM = .type.EIM ,
           zero  = .zero )
      
    }, list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
             .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
             .edrift = edrift , .esd = esd , .evar = evar , 
             .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
             .type.EIM = type.EIM , .type.likelihood = type.likelihood ,
             .Order = Order, .nodrift = nodrift ,
             .zero = zero , .eimARff = eimARff ,
             .realM1 = realM1 ))),
          
          
          
          
          
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
      maOrder <- .OrderMA
      # Conformability condition #
      if ( length( .Order ) > NOS )
        stop("Non-conformable orders")
      nOrder <- rep( .Order , length = NOS)[1:NOS]
      extra$nOrder <- nOrder
      M1 <- (2 - .nodrift) + nOrder + maOrder
      
      for (jj in 1:NOS)
        w[1:nOrder[jj], jj] <- 1e0
      m.max  <- max(nOrder)
      m1.max <- max(M1)
      
      if ( length( .iARcoeff ) ) {
        
        if( !Is.Numeric( .iARcoeff ) )
          stop("'iARcoeff' must be a vector with numeric entries.") 
        
        if ( length( .iARcoeff ) != sum( nOrder ) )  
          stop("Invalid number of entries in 'iARcoeff'.")
        
        init.the <- init.theMA <-  vector("list", length = NOS)
        auuxx <- 0
        for ( ii in 1:NOS ) {
          niAux <- .iARcoeff[ auuxx + 1:nOrder[ii] ] 
          initChecks <- checkTS.ffs(thetaEst = c( niAux ),
                                    tsclass  = "AR",
                                    NofS     = 1,
                                    chOrder  = nOrder[ii],
                                    pRoots   = FALSE,
                                    retmod   = TRUE)
          if( any( initChecks < 1 + 1e-2 ) ) 
            stop("Initial values of the AR coefficients for response ",
                 ii, " do not belong to a stationary model.")
          
          init.the[[ii]] <- matrix(niAux, nrow = n , ncol = nOrder[ii],
                                   byrow = TRUE)
          auuxx <- auuxx + nOrder[ii]
        }
      }     
      
      temp.name <- if ( ( .flagAR ) ) "Mean" else
         if ( ( .flagMA ) ) "drift" else "drift.mean"
      
      names.1 <- if ( .nodrift ) NULL else
        paste(temp.name, 1:NOS, sep = "")
      
      rm(temp.name)
      names.2 <- if ( .var.arg ) paste("noiseVar", 1:NOS, sep="") else
        paste("noiseSD", 1:NOS, sep = "")
      names.3 <- names.4 <- character(0)
      for (jj in 1:m.max) 
        names.3 <- c(names.3, paste(paste("ARcoeff", jj, sep = ""),
                                    1:NOS, sep = ""))
      
      for (jj in 1:maOrder) 
        names.4 <- c(names.4, paste(paste("MAcoeff", jj, sep = ""),
                                    1:NOS, sep = ""))
      
      pars.names <- c(names.1, names.2, names.3, names.4)
      pars.names <- pars.names[interleave.VGAM(m1.max * NOS, M1 = m1.max )]
      pars.names <- matrix(pars.names, NOS, ncol = m1.max, byrow = TRUE)
      
      pres.names <- c(if ( .nodrift ) NULL else
        namesof(names.1, .ldrift, .edrift, tag = FALSE),
        namesof(names.2, if ( .var.arg ) .lvar else .lsd ,
                if ( .var.arg ) .evar else .esd , tag = FALSE),
        namesof(names.3, .lARcoeff , .eARcoeff , tag = FALSE),
        namesof(names.4, .lMAcoeff , .eMAcoeff , tag = FALSE))
      
      pres.names <- pres.names[interleave.VGAM(m1.max * NOS, M1 = m1.max)]
      pres.names <- matrix(pres.names, NOS, ncol = m1.max, byrow = TRUE)
      parameters.names <- predictors.names <- character(0)
      
      for (ii in 1:NOS) {
        parameters.names <- c(parameters.names, pars.names[ii, 1:M1[ii]])
        predictors.names <- c(predictors.names, pres.names[ii, 1:M1[ii]])
      }
      
      ts.w <- ts2.w <- NULL
      
      if ( ( .flagAR ) )
        ts.w <- grep("ARcoeff", parameters.names)
      
      if ( ( .flagMA ) )
        ts2.w <- grep("MAcoeff", parameters.names)
      
      if (length(ts.w) || length( ts2.w)) {
        parameters.names <- parameters.names[-unique(c(ts.w, ts2.w))]
        predictors.names <- predictors.names[-unique(c(ts.w, ts2.w))]  
      }
      
      nidrift <- rep( .idrift , NOS)[1:NOS]
      nisd    <- rep( .isd , NOS)[1:NOS]
      nivar   <- rep( .ivar , NOS)[1:NOS]
      
      if ( !length( .iARcoeff ))
        init.the <- vector("list", length = NOS)
      
      if ( !length( .iMAcoeff ))
        init.theMA <- vector("list", length = NOS)
      
      if (!length(etastart)) {
        
        init.mean <- if ( length ( nidrift ) )
          matrix( nidrift, nrow = n , ncol = NOS, byrow = TRUE ) else
            matrix( 0.01 , nrow = n , ncol = NOS) 
        
        init.sd <- if ( length ( nisd ) )
          matrix( nisd, nrow = n , ncol = NOS, byrow = TRUE ) else
            matrix( 1.0 , nrow = n , ncol = NOS) 
        
        init.var <- if ( length ( nivar ) )
          matrix( nivar , nrow = n , ncol = NOS, byrow = TRUE ) else
            matrix( 1.0 , nrow = n , ncol = NOS) 
        
        for (spp. in 1:NOS) {
          GammaS <- cross.gammas(x = y[, spp. ], lags = nOrder[spp. ])
          myToeplitz <- toeplitz(GammaS[1:nOrder[ spp. ] ] )
          myresp     <- c(GammaS[2:(nOrder[ spp. ] + 1)])
          initTheta  <- justHelp <- solve(myToeplitz, myresp)
          
          y.res <- cbind(y , WN.lags(y = y, lags = .Order ), 
                         WN.lags(y = cbind(extra$res), lags = .OrderMA))
          
          ini.fit <- lsfit(x = y.res[, -1, drop = FALSE],
                           y = y.res[,  1, drop = FALSE],
                           intercept = TRUE)  # nodrift = FALSE...
          
          init.MA <- coef(ini.fit)[-(1:(1 + .Order))]
          #extra$res <- residuals(ini.fit) ## Commented out on 20171221
          
          if (length(initTheta) < m.max)
            initTheta <- c(initTheta, rep(0, m.max - length(initTheta)))
          
          if ( !length( .iARcoeff ) ) 
            init.the[[spp.]] <- matrix(initTheta, nrow = n,
                                 ncol = length(initTheta), byrow = TRUE)
          
          if ( !length( .iMAcoeff ) ) 
            init.theMA[[spp.]] <- matrix(init.MA, nrow = n,
                                   ncol = length(init.MA), byrow = TRUE)
          
          if (  !length( .idrift )  ) 
            init.mean[, spp. ] <- mean(y[, spp. ]) * (1 - sum(initTheta))
          
          justHelp <- GammaS[2:(nOrder[spp.] + 1)] * justHelp
          if (  !length( .ivar )  ) 
            init.var[, spp. ] <- max(0.01, GammaS[1] - (sum(justHelp)) )
          
          if (max(0.01, GammaS[1] - (sum(justHelp)) ) == 0.01 )
            warning("The noise variance is too close to zero.")
          
          if (  !length( .isd )  ) {
            init.sd[, spp. ] <- 
              max( sqrt(0.01), sqrt(init.var[, spp. ]))
            if (max( sqrt(0.01), sqrt(init.var[, spp. ]))  ==  sqrt(0.01))
              warning("The noise standard deviation",
                      " is too close to zero.")
          } 
        }   
        
        etastart <- cbind(if (.nodrift ) NULL else 
          theta2eta(init.mean, .ldrift , earg = .edrift),
          if ( .var.arg )
            theta2eta(init.var, .lvar , earg = .evar ) else 
              theta2eta(init.sd, .lsd , earg = .esd ))
        
        pre.eta <- numeric(0)
        for (kk in 1:NOS) 
          pre.eta <- cbind(pre.eta, theta2eta(cbind(init.the[[kk]]),
                                           .lARcoeff , earg = .eARcoeff ))
        for (kk in 1:NOS) 
          pre.eta <- cbind(pre.eta, theta2eta(cbind(init.theMA[[kk]]),
                                           .lMAcoeff , earg = .eMAcoeff ))
        pre.eta <- pre.eta[, interleave.VGAM((m.max + .OrderMA) * NOS,
                                             M1 = m.max + .OrderMA,
                                             inverse = TRUE), drop = FALSE]
      } else {
        pre.eta <- NULL
      }
      
      etastart <- cbind(etastart, pre.eta)
      etastart <- etastart[, interleave.VGAM(m1.max * NOS, M1 = m1.max)]
      if ( ( .flagAR ) || ( .flagMA ) ) 
        etastart <- etastart[, -unique(c(ts.w, ts2.w)), drop = FALSE]
      
      
      etastart
      
    }), list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
              .lARcoeff = lARcoeff, .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .esd = esd , .evar = evar ,
              .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
              .idrift = idrift, .ivar = ivar , .isd = isd ,
              .iARcoeff = iARcoeff , .iMAcoeff = iMAcoeff ,
              .nodrift  = nodrift , .var.arg = var.arg ,
              .type.likelihood = type.likelihood ,
              .Order = Order , .OrderMA = OrderMA ,
              .flagMA = flagMA , .flagAR = flagAR ,
              .ARord = ARorder , .MAord = MAorder ))),
          
          
          
          
          
    linkinv = eval(substitute(function(eta, extra = NULL ) {
      
      eta <- cbind(eta)
      n   <- nrow(eta)
      NOS <- extra$NOS
      nOrder <- extra$nOrder
      M1  <- 2 - .nodrift + .ARord + .MAord
      M   <- 2 + .ARord + .MAord #sum(M1) - nOrder - .OrderMA
      NNro <- 1
      
      if ( ( .flagAR ) && ( .flagMA ) ) {
        
        fitted.ts <- eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
        if (.var.arg ) {
          arp.var <- eta2theta(eta[, 2 - ( .nodrift ), drop = FALSE],
                               .lvar , .evar )
          arp.sd  <- sqrt(arp.var)
        } else {
          arp.sd <- eta2theta(eta[, 2 - ( .nodrift ), drop = FALSE],
                              .lsd , .esd )
          arp.var  <- arp.sd^2
        }
        
      } else {
        
        arp.drMean <- eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
        
        if ( .var.arg ) {
          arp.var <- eta2theta(eta[, 2, drop = FALSE], .lvar, .evar )
          arp.sd  <- sqrt(arp.var)
        } else {
          arp.sd  <- eta2theta(eta[, 2, drop = FALSE], .lsd , .esd )
          arp.var <- arp.sd^2
        }
        
        arp.The <- matrix(0.0, n, .Order )
        maq.The <- matrix(0.0, n, .OrderMA )
        
        kk <- 0
        if ( !( .flagAR ) ) {
          arp.The <- eta2theta(eta[, 3:(2 + .Order ), drop = FALSE],
                               .lARcoeff , .eARcoeff )
          kk <- .Order
        }
        
        
        if ( !( .flagMA ) )
          maq.The <- eta2theta(eta[, -(1:( 2 + kk )), drop = FALSE], 
                               .lMAcoeff , .eMAcoeff )
        
        y.lags  <- WN.lags(cbind(extra$y), lags = .Order )
        ma.lags <- WN.lags(cbind(extra$res), lags = .OrderMA )
        
        fitted.ts <- rowSums(arp.The * y.lags) +
          rowSums(maq.The * ma.lags) + arp.drMean
        #fitted.ts <- arp.drMean +
        #  apply(fitted.ts, c(1, 3), function(x) sum(x))
      }
      
      fitted.ts <- fitted.ts + extra$res
      fitted.ts[c(1:NNro), 1] <- (NNro + 0.10) * extra$y[1]
      
      fitted.ts
      
    }, list ( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
              .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .evar = evar , .esd = esd ,
              .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
              .OrderMA = OrderMA , .Order = Order ,
              .var.arg = var.arg , .nodrift = nodrift , 
              .flagAR = flagAR , .flagMA = flagMA ,
              .ARord = ARorder , .MAord = MAorder ))),
          
          
          
          
          
    last = eval(substitute(expression({
      eta <- cbind(eta)
      n   <- nrow(eta)
      NOS <- extra$NOS
      nOrder <- extra$nOrder
      M1  <- 2 - .nodrift + .ARord + .MAord
      M   <- 2 + .ARord + .MAord #sum(M1) - nOrder - .OrderMA
      
      arp.drMean <- eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
      
      if ( .var.arg ) {
        arp.var <- eta2theta(eta[, 2, drop = FALSE], .lvar , .evar )
        arp.sd  <- sqrt(arp.var)
      } else {
        arp.sd  <- eta2theta(eta[, 2, drop = FALSE], .lsd , .esd )
        arp.var <- arp.sd^2
      }
      
      misc$link <- character(0)
      misc$earg <- vector("list", length = 0 )
      
      kk <- 0
      misc$link[1]   <- .ldrift
      misc$earg[[1]] <- .edrift
      misc$link[2]   <- if ( .var.arg ) .lvar else .lsd
      misc$earg[[2]] <- if ( .var.arg ) .evar else .esd
      
      
      if ( !( .flagAR ) ) {
        misc$link <- c(misc$link, rep( .lARcoeff , .Order ))
        for (ii in 3:( 2 + .Order )) 
          misc$earg[[ii]] <- .eARcoeff  
        kk <- .Order
      }
      
      
      if ( !( .flagMA ) ) {
        misc$link <- c(misc$link, rep( .lMAcoeff , .OrderMA ))
        for (jj in (3 + kk):(2 + kk + .OrderMA))
          misc$earg[[jj]] <- .eMAcoeff
      }
      
      
      names(misc$link) <- parameters.names
      names(misc$earg) <- parameters.names
      
      my.process <- "pureGARCH"
      if ( ( .flagMA ) && !( .flagAR ) ) {
        my.process <- "AR"
      } else {
        if ( !( .flagMA ) && ( .flagAR) ) {
          my.process <- "MA"
        } else {
          if ( !( .flagMA ) && !( .flagAR)  ) {
            my.process <- "ARMA"
          } else {
            
          }
        }
      }
      
      pre.Order <- c( .ARord , .MAord )
      
      if ( ( .flagMA ) && !( .flagAR ) ) {
        # nOrder is the AR order
        pre.Order <- nOrder
      } else {
        if ( !( .flagMA ) && ( .flagAR ) ) {
          pre.Order <- .MAord
        }  else {
        }
      }
      
      misc$var.arg  <- .var.arg
      misc$expected <- TRUE
      misc$process  <- my.process
      misc$OrderMA  <- .OrderMA
      misc$Order    <- pre.Order
      rm(pre.Order, my.process)
      misc$varT     <- arp.var # extra$Sigmas
      misc$zero     <- .zero
      misc$NOS      <- NOS
      misc$M1       <- M1
      misc$M        <- M
      misc$ARord    <- .ARord
      misc$MAord    <- .MAord
      misc$flagAR   <- .flagAR
      misc$flagMA   <- .flagMA
      misc$nomean   <- .nodrift
      misc$eimARff  <- .eimARff
      misc$residuals <- extra$res
      misc$theta.names       <- parameters.names
      misc$control$wzepsilon <- control$wzepsilon
      misc$type.likelihood   <- .type.likelihood
      misc$multipleResponses <- FALSE
      fit$prior.weights <- w  #20180315
      
      noChecks <- .noChecks
      if (( .flagAR ) && (.flagMA ) && !( .noChecks )) {
        noChecks <- TRUE
      cat("\n No checks on stationarity/invertibility performed, since\n",
              "no ARMA component involved. \n")
      }
      
      
      if ( !( .noChecks ) ) {
        
        saveNames <- names( fit$coefficients )
        Flag <- grep("(Intercept)", names(fit$coefficients) )
        
        if (length(Flag) > 2) {
          
          fin.coe <- numeric(0)
          # At this point, No constraints have been established.
          myCoeffs <- coef(fit)[Flag][-(1:2)]
          flag1   <- TRUE
          
          if ( !( .flagAR )) {
            Ar.c <- eta2theta(myCoeffs[1:nOrder], .lARcoeff, .eARcoeff)
            ARroots <- checkTS.ffs(thetaEst = Ar.c,  # a vector.
                                   tsclass  = "AR",
                                   NofS     = NOS,
                                   pRoots   = FALSE,
                                   chOrder  = nOrder)
            
            myCoeffs <- myCoeffs[-(1:nOrder)]
            flag1 <- (all(ARroots > 1 - 1e-4) && flag1)
            
          }
          
          if ( !( .flagMA )) {
            
            Ma.c <- eta2theta(myCoeffs, .lMAcoeff, .eMAcoeff)
            MAroots <- checkTS.ffs(thetaEst = Ma.c,  # a vector.
                                   tsclass  = "MA",
                                   NofS     = NOS,
                                   pRoots   = FALSE,
                                   chOrder  = .OrderMA )
            flag1 <- (all(MAroots > 1 - 1e-4) && flag1)
          } 
          
          
          
          
          if (flag1)
            cat("\nChecks on stationarity / invertibility successfully",
                "performed. \nNo roots liying inside the unit circle.",
                "\nFurther details within the 'summary' output.\n") else
                  
                cat("\nChecks on stationarity / invertibility successfully",
                    "performed. \nNOTE: Some roots lie inside the unit ",
                    "circle. \nFurther details within the 'summary' ", 
                    "output.\n")
        } else {
          if (FALSE)
            cat("\nApparently, some constraints using 'cm.ARMA' have been",
                "set. \nDetails on stationarity / invertibility checks",
                "whitin the \n'summary' output.")
        } # End of (length(Flag) > 2)
      }   # End ( !(.noChecks ) )
    }), list( .ldrift = ldrift, .lsd = lsd , .lvar = lvar , 
              .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
              .edrift = edrift , .esd = esd , .evar = evar , 
              .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
              .flagAR = flagAR , .flagMA = flagMA ,
              .var.arg = var.arg , .nodrift = nodrift ,
              .type.likelihood = type.likelihood ,
              .eimARff = eimARff ,
              .Order = Order , .OrderMA = OrderMA ,
              .zero = zero , .noChecks = noChecks ,
              .ARord = ARorder , .MAord = MAorder ))),
          
          
          
          
          
          
    loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE,
                eta, extra = NULL, summation = TRUE) {
      
      eta <- cbind(eta)
      NOS <- extra$NOS
      y   <- extra$y
      nOrder <- extra$nOrder
      M1  <- 2 - .nodrift + .ARord + .MAord
      M   <- 2 + .ARord + .MAord #sum(M1) - nOrder - .OrderMA
      n   <- nrow(eta)
      
      if ( ( .flagAR ) && ( .flagMA ) ) {
        mean.ts <- if ( .nodrift ) matrix(0.0, n, NOS) else 
          eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
        
        if (.var.arg ) {
          arp.var <- eta2theta(eta[, 2 - ( .nodrift ), drop = FALSE],
                               .lvar , .evar )
          arp.sd  <- sqrt(arp.var)
        } else {
          arp.sd <- eta2theta(eta[, 2 - ( .nodrift ), drop = FALSE],
                              .lsd , .esd )
          arp.var  <- arp.sd^2
        }
        
      } else {
        
        arp.drMean <- eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
        
        if ( .var.arg ) {
          arp.var <- eta2theta(eta[, 2, drop = FALSE], .lvar, .evar )
          arp.sd  <- sqrt(arp.var)
        } else {
          arp.sd  <- eta2theta(eta[, 2, drop = FALSE], .lsd , .esd )
          arp.var <- arp.sd^2
        }
        
        arp.The <- matrix(0.0, n, .Order )
        maq.The <- matrix(0.0, n, .OrderMA )
        kk <- 0
        
        if ( !( .flagAR ) ) {
          arp.The <- eta2theta(eta[, 3:(2 + .Order ), drop = FALSE],
                               .lARcoeff , .eARcoeff )
          kk <- .Order
        }
        
        
        if ( !( .flagMA ) )
          maq.The <- eta2theta(eta[, -(1:( 2 + kk )), drop = FALSE], 
                               .lMAcoeff , .eMAcoeff )
        
        y.lags  <- WN.lags(cbind(extra$y), lags = .Order )
        ma.lags <- WN.lags(cbind(extra$res), lags = .OrderMA )
        mean.ts <- rowSums(arp.The * y.lags) +
          rowSums(maq.The * ma.lags) + arp.drMean
        #mean.ts <- array(arp.The * y.lags, dim = c(n, max(nOrder), NOS))
        #mean.ts <- arp.drMean + apply(mean.ts, c(1, 3), sum)
        names(mean.ts) <- NULL
      }
      
      
      if (residuals) {
        stop("Loglikelihood not implemented yet", 
             " to handle residuals.")
      } else {
        loglik.terms <-
          c(w) * dARMA(x = y, mean = mean.ts,
                       sd = arp.sd, log = TRUE)
        sum( loglik.terms )
      }
      
    }, list (  .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
               .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
               .edrift = edrift , .esd = esd , .evar = evar ,
               .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
               .var.arg = var.arg , .nodrift = nodrift ,
               .Order = Order , .OrderMA = OrderMA ,
               .type.likelihood = type.likelihood ,
               .flagAR = flagAR , .flagMA = flagMA ,
               .ARord = ARorder , .MAord = MAorder ))),
          
          
          
          
          
    validparams = eval(substitute(function(eta, y, extra = NULL) {
      
      NOS <- ncol(y)
      eta <- cbind(eta)
      n   <- nrow(eta)
      nOrder <- extra$nOrder
      M1  <- 2 - .nodrift + .ARord + .MAord
      M   <- 2 + .ARord + .MAord #sum(M1) - nOrder - .OrderMA
      
      if ( ( .flagAR ) && ( .flagMA ) ) {
        
        mean.ts <- if ( .nodrift ) matrix(0.0, n, NOS) else 
          eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
        
        if (.var.arg ) {
          arp.var <- eta2theta(eta[, 2, drop = FALSE], .lvar , .evar )
          arp.sd  <- sqrt(arp.var)
        } else {
          arp.sd <- eta2theta(eta[, 2, drop = FALSE], .lsd , .esd )
          arp.var  <- arp.sd^2
        }
        
      } else {
        
        arp.drMean <- eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
        
        if ( .var.arg ) {
          arp.var <- eta2theta(eta[, 2, drop = FALSE], .lvar , .evar )
          arp.sd  <- sqrt(arp.var)
        } else {
          arp.sd  <- eta2theta(eta[, 2, drop = FALSE], .lsd , .esd )
          arp.var <- arp.sd^2
        }
        
        arp.The <- matrix(0.0, n, .Order )
        maq.The <- matrix(0.0, n, .OrderMA )
        kk <- 0
        
        if ( !( .flagAR ) ) {
          
          arp.The <- eta2theta(eta[, 3:(2 + .Order ), drop = FALSE],
                               .lARcoeff , .eARcoeff )
          kk <- .Order
        }
        
        if ( !( .flagMA ) )
          maq.The <- eta2theta(eta[, -(1:( 2 + kk )), drop = FALSE], 
                               .lMAcoeff , .eMAcoeff )
        
        y.lags  <- WN.lags(cbind(extra$y), lags = .Order )
        ma.lags <- WN.lags(cbind(extra$res), lags = .OrderMA )
        mean.ts <- rowSums(arp.The * y.lags) + 
          rowSums(maq.The * ma.lags) + arp.drMean
        #mean.ts <- array(arp.The * y.lags, dim = c(n, max(nOrder), NOS))
        #mean.ts <- arp.drMean + apply(mean.ts, c(1, 3), sum)
        names(mean.ts) <- NULL
      }
      
      okay1   <- (all(is.finite(arp.var)) && all(arp.var > 0) &&
                    all(is.finite(mean.ts)))
      
      okay1
      
    }, list( .ldrift = ldrift , .lsd = lsd , .lvar = lvar , 
             .lARcoeff = lARcoeff , .lMAcoeff = lMAcoeff ,
             .edrift = edrift , .esd = esd , .evar = evar ,
             .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
             .var.arg = var.arg , .nodrift = nodrift , 
             .Order = Order , .OrderMA = OrderMA ,
             .type.likelihood = type.likelihood ,
             .flagAR = flagAR , .flagMA = flagMA ,
             .ARord = ARorder , .MAord = MAorder ))),
          
          
          
          
    #vfamily = c("vgtsff", "vgltsff-class", "ARff", "GARCHff"),
   vfamily = c("ARMAX.GARCHff", "GARCHff", 
               "ARMAX.GARCHvgltsmff", "ARMAXvgltsmff"),
          
          
          
   deriv = eval(substitute(expression({
     
     eta <- cbind(eta)
     NOS <- ncol(y)
     n   <- nrow(eta)
     nOrder <- extra$nOrder
     M1  <- 2 - .nodrift + .ARord + .MAord
     M   <- 2 + .ARord + .MAord #sum(M1) - nOrder - .OrderMA
     
     
     if ( ( .flagAR ) && ( .flagMA ) ) {
       arp.drMean <-  if ( .nodrift ) matrix(0.0, n, NOS) else
         eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
       
       if (.var.arg ) {
         arp.var <- eta2theta(eta[, 2 - ( .nodrift ), drop = FALSE],
                              .lvar , .evar )
         arp.sd  <- sqrt(arp.var)
       } else {
         arp.sd <- eta2theta(eta[, 2 - ( .nodrift ), drop = FALSE],
                             .lsd , .esd )
         arp.var  <- arp.sd^2
       }
       
       y.means <- y - arp.drMean
       
       dl.drfmean <- y.means / arp.var
       if ( .var.arg ) {
         dl.dvar <- y.means^2 / (2 * arp.var^2) - 1 / (2 * arp.var)
       } else {
         dl.dsd <- y.means^2 / arp.sd^3 - 1 / arp.sd
       }
       
       drfm.deta <- dtheta.deta(arp.drMean, .ldrift , earg = .edrift )
       if ( .var.arg ) {
         dvar.deta <- dtheta.deta(arp.var , .lvar , earg = .evar )
         dVarSD    <- dvar.deta
       } else {
         dsd.deta  <- dtheta.deta(arp.sd  , .lsd  ,  earg = .esd )
         dVarSD    <- dsd.deta
       }
       
       dl.deta1 <- c(w) * dl.drfmean * drfm.deta
       if ( .var.arg ){
         dl.deta2 <- c(w) * dl.dvar * dvar.deta
       } else {
         dl.deta2 <- c(w) * dl.dsd * dsd.deta
       }
       
       fDeriv <- cbind( if ( .nodrift ) NULL else dl.deta1, dl.deta2)
       
     } else {
       
       arp.drMean <- eta2theta(eta[, 1, drop = FALSE], .ldrift , .edrift )
       
       if ( .var.arg ) {
         arp.var <- eta2theta(eta[, 2, drop = FALSE], .lvar , .evar )
         arp.sd  <- sqrt(arp.var)
       } else {
         arp.sd  <- eta2theta(eta[, 2, drop = FALSE], .lsd , .esd )
         arp.var <- arp.sd^2
       }
       
       arp.The <- matrix(0.0, n, .Order )
       maq.The <- matrix(0.0, n, .OrderMA )
       kk <- 0
       
       if ( !( .flagAR ) ) {
         arp.The <- eta2theta(eta[, 3:(2 + .Order ), drop = FALSE],
                              .lARcoeff , .eARcoeff)
         kk <- .Order
       }
       
       if ( !( .flagMA ) )
         maq.The <- eta2theta(eta[, -(1:( 2 + kk )), drop = FALSE], 
                              .lMAcoeff , .eMAcoeff)
       y.lags  <- WN.lags(cbind(extra$y), lags = .Order )
       ma.lags <- WN.lags(cbind(extra$res), lags = .OrderMA )
       
       y.means <- y - (arp.drMean + rowSums(arp.The * y.lags) +
                         rowSums(maq.The * ma.lags))
       
       #y.means <- array(arp.The * y.lags, dim = c(n, max(nOrder), NOS))
       #y.means <- y - (arp.drMean + 
       #                  apply(y.means, c(1, 3), function(x) sum(x)))
       
       dl.drfmean <- y.means / arp.var
       if ( .var.arg ) {
         dl.dvar <- y.means^2 / (2 * arp.var^2) - 1 / (2 * arp.var)
       } else {
         dl.dsd <- y.means^2 / arp.sd^3 - 1 / arp.sd
       }
       
       dl.dThe <- y.lags * matrix(y.means, n, .Order ) / 
         matrix(arp.var, n, .Order )
       dl.dPhi <- ma.lags * matrix(y.means, n, .OrderMA ) /
         matrix( arp.var , n, .OrderMA)
       #y.lags <- array(y.lags, dim = c(n, max(nOrder), NOS))
       #dl.dThe <- array(NA_real_, dim = c(n, max(nOrder), NOS))
       #for (ii in 1:NOS)
       #  dl.dThe[, , ii] <- (array( y.means, dim = c(n, 1, NOS))[, , ii] *
       #                        y.lags[, , ii]) / arp.var[, ii]
       
       #dl.dThe <- interleaveArray.VGAMextra(dl.dThe)
       
       # Partial derivatives of theta wrt eta #
       
       drfm.deta <- dtheta.deta(arp.drMean, .ldrift , earg = .edrift )
       if ( .var.arg ) {
         dvar.deta <- dtheta.deta(arp.var , .lvar , earg = .evar )
         dVarSD    <- dvar.deta
       } else {
         dsd.deta  <- dtheta.deta(arp.sd  , .lsd  ,  earg = .esd )
         dVarSD    <- dsd.deta
       }
       
       #dThedEta <- interleaveArray.VGAMextra(array(arp.The,
       #                                        dim = c(n, ordMax, NOS)))
       dThedEta <- dtheta.deta(arp.The , .lARcoeff , earg = .eARcoeff)
       dPhidEta <- dtheta.deta(maq.The , .lMAcoeff , earg = .eMAcoeff)
       
       # The chain rule to obtain the total derivative of 'l' wrt
       # eta: (dl/deta_i) = (dl / dtheta_i) * (dtheta_i / deta_i)
       
       if ( .var.arg ){
         dl.deta2 <- c(w) * dl.dvar * dvar.deta
       } else {
         dl.deta2 <- c(w) * dl.dsd * dsd.deta
       }
       
       fDeriv <- cbind(c(w) * dl.drfmean * drfm.deta,
                       if ( .var.arg ) c(w) * dl.dvar * dvar.deta else
                         c(w) * dl.dsd * dsd.deta,
                       if ( .flagAR ) NULL else c(w) * dl.dThe * dThedEta,
                       if ( .flagMA ) NULL else c(w) * dl.dPhi * dPhidEta)
       
       #a.help  <- interleaveArray.VGAMextra(c(w) * dl.dThe * dThedEta,
       #                                     inverse = TRUE)
       
       #fDeriv <- numeric(0)
       #for (ii in 1:NOS) {
       # fDeriv <- cbind(fDeriv, if ( .nodrift ) NULL else  dl.deta1[, ii],
       #                 dl.deta2[, ii],
       #                 matrix(c(a.help[, 1:nOrder[ii] , ii]), nrow = n,
       #                        ncol = nOrder[ii]))
       #}
     }
     
     rownames(fDeriv) <- NULL
     
     fDeriv
     
   }), list( .ldrift = ldrift , .lvar = lvar , .lsd = lsd ,
             .lARcoeff  = lARcoeff , .lMAcoeff = lMAcoeff ,
             .edrift = edrift , .evar = evar , .esd = esd ,
             .eARcoeff = eARcoeff , .eMAcoeff = eMAcoeff ,
             .var.arg = var.arg , .nodrift = nodrift ,
             .Order = Order , .OrderMA = OrderMA ,
             .flagAR = flagAR , .flagMA = flagMA ,
             .ARord = ARorder , .MAord = MAorder ))),
          
          
          
          
    weight = eval(substitute(expression({
      
      if ( !( .flagAR ) && !( .flagMA ) ) {
        exact.eim <- ARMA.EIM.G2(y = cbind(extra$y),
                                 arCoe   = cbind(arp.The),
                                 maCoe   = cbind(maq.The),
                                 estRes  = cbind(extra$res),
                                 sdError = cbind(arp.sd),
                                 var.arg = .var.arg ,
                                 order   = c(.Order , .OrderMA ),
                                 nodrift = .nodrift )
        
        comb.wz <- combVGAMextra(1:M1)
        a.bind  <- cbind(if ( .nodrift ) NULL else drfm.deta ,
                         if (.var.arg ) dvar.deta else dsd.deta, 
                         dThedEta, dPhidEta)
        
        dthdeta <- apply(comb.wz, 1, function(x) {
          a.bind[, x[1]] * a.bind[, x[2]]
        })
        
        wz <- c(w) * exact.eim * dthdeta
        
      }
      
      if ( !( .flagAR ) && ( .flagMA ) ) {
        
        exact.eim <-  ARpEIM.G2(y = cbind(extra$y),
                                drift = arp.drMean,
                                sdError = cbind(arp.sd),
                                ARpcoeff = cbind(arp.The),
                                var.arg = .var.arg ,
                                order = .Order ,
                                nodrift = .nodrift )
        
        comb.wz <- combVGAMextra(x = 1:(2 + nOrder))
        a.bind  <- cbind(if ( .nodrift ) NULL else drfm.deta ,
                         if ( .var.arg ) dvar.deta else dsd.deta, 
                         dThedEta)
        
        dthdeta <- apply(comb.wz, 1, function(x) {
          a.bind[, x[1]] * a.bind[, x[2]]
        })
        
        wz <- c(w) * exact.eim * dthdeta
      }
      
      if ( ( .flagAR ) && !( .flagMA )) {
        exact.eim <- MAqEIM.G2(y = cbind(extra$y),
                               mean = arp.drMean,
                               sdError = cbind(arp.sd),
                               MAqcoeff = cbind(maq.The),
                               var.arg = .var.arg ,
                               order = .OrderMA ,
                               nomean = .nodrift )
        
        comb.wz <- combVGAMextra(1:(2 + .OrderMA ))
        a.bind  <- cbind(if ( .nodrift ) NULL else drfm.deta ,
                         if (.var.arg ) dvar.deta else dsd.deta, 
                         dPhidEta)
        
        dthdeta <- apply(comb.wz, 1, function(x) {
          a.bind[, x[1]] * a.bind[, x[2]]
        })
        
        wz <- c(w) * exact.eim * dthdeta
        
      }
      
      if ( ( .flagAR ) && ( .flagMA )) {
        ned2l.dsmn   <- 1 / arp.sd^2
        ned2l.dvarSD <- if ( .var.arg ) 1 / (2 * arp.var^2) else 2 /arp.var
        
        wz <- c(w) * cbind( if ( .nodrift ) NULL else 
          ned2l.dsmn * drfm.deta^2, ned2l.dvarSD * dvar.deta^2)
      }
      
      wz
      
    }), list( .var.arg  = var.arg , .idrift = idrift , 
              .Order   = Order ,  .OrderMA = OrderMA ,
              .nodrift  = nodrift , .eimARff = eimARff ,
              .flagAR = flagAR , .flagMA = flagMA ,
              .ivar = ivar , .isd = isd ,
              .iARcoeff = iARcoeff ))) ) # End of new("vgltsmff")
    alo
}  # End ARMAGARCH







## ARXff control function.
ARXff.control <- function(save.weights = TRUE,
                          summary.HDEtest = FALSE,...) { 
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}

ARMAX.GARCHff.control <- function(save.weights = TRUE,
                                  summary.HDEtest = FALSE,...) { 
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}




# ARffEIM Back up using matrix algebra
if (FALSE)
  ARpEIM.G2 <- function(y, drift, sdError, ARpcoeff,
                        var.arg = TRUE, order = 1,
                        nodrift = FALSE) {
    
    Order <- order; rm(order)
    y <- cbind(y); nn <- nrow(y)
    ARpcoeff <- matrix(ARpcoeff, nrow = nn, ncol = Order)
    
    RMat  <- diag(c(sdError^2))
    dRMat <- cbind( if(!nodrift) matrix(0, nrow = nn, ncol = 1) else NULL, 
                    if (var.arg) matrix(1, nrow = nn, ncol = 1) else
                      2 * sdError,  
                    matrix(0, nrow = nn, ncol = Order))
    
    y.lags <- WN.lags(y = y, lags = Order,
                      to.complete = 1 / (1 - ARpcoeff[1, ])^2)
    
    dmt0 <- matrix(1, nrow = nn, ncol = 1)
    dmt  <- cbind(if (!nodrift) dmt0 else NULL, rep_len(0, nn), y.lags)
    
    M <- 2 + Order  # For ARp
    try.comb <- combVGAMextra(1:M, nodrift = nodrift)
    
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

