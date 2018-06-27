###########################################################################
# These functions are
# Copyright (C) 2014-2018 V. Miranda ans T. W. Yee, University of Auckland
# All rights reserved.
#
# Function 'invweibull2mr'.
invweibull2mr    <- function(lscale  = "loge", 
                             lshape  = logoff(offset = -2),
                             iscale  = NULL, 
                             ishape  = NULL, 
                             imethod = 2,
                             lss     = TRUE,
                             gscale  = exp(-4:4), 
                             gshape  = exp(-4:4),
                             probs.y = c(0.25, 0.50, 0.75),
                             zero    = "shape" ) {
  
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("Bad input for argument 'iscale'.")
  
  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("Bad input for argument 'ishape'.")
  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE,
                  positive = TRUE) || imethod > 2)
    stop("Bad input for argument 'imethod'.")
  
  if (length(lss) && !is.logical(lss))
    stop("Bad input for argument 'lss'.")
  
  if (length(probs.y) < 2 || max(probs.y) > 1 || 
        !is.Numeric(probs.y, positive = TRUE))
    stop("Bad input for argument 'probs.y'.")
  
  if (!is.Numeric(gscale) || min(gscale) < 0)
    stop ("Bad input for argument 'gscale'.")
  
  if (!is.Numeric(gshape) || min(gshape) < 0)
    stop("Bad input for argument 'gshape'.")
  
  lscale    <- as.list(substitute(lscale))
  escale    <- link2list(lscale)
  lscale    <- attr(escale, "function.name")
  
  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")
  
  alo <- 
    new("vglmff", 
        blurb = c("2-parameter Inverse Weibull Distribution. \n", 
                  "Links:    ", 
                  ifelse(lss, 
                         namesof("scale", lscale, earg = escale),
                         namesof("shape", lshape, earg = eshape)), ", ",
                  ifelse(lss, 
                         namesof("shape", lshape, earg = eshape),
                         namesof("scale", lscale, earg = escale)), "\n ",
                  "Mean:     scale*[gamma(1 - 1/shape)],  \n",
                  "Variance: scale^2 * [gamma(1 - 2/shape) - ", 
                  "(gamma(1 - 1/shape)]^2), if shape > 2."),
        
    constraints = eval(substitute(expression({
      M1  <- 2
      constraints <- 
        cm.zero.VGAM(constraints, x = x, zero = .zero , M = M,
                     # old predictor.names 20170226
                     predictors.names = parameter.names,
                     M1 = M1)
      
    }), list( .zero = zero ))),
    
    infos = eval(substitute(function(...) {
      
      list(M1 = 2, 
           Q1 = 1,
           expected = TRUE, 
           multipleResponses = TRUE, 
           lscale = .lscale , 
           lshape = .lshape , 
           zero = .zero )
      
    }, list( .lscale = lscale, .lshape = lshape, 
             .zero = zero ))),
   
    initialize = eval(substitute(expression({
      M1 <- 2  
      check <- w.y.check(w = w, y =y, 
                         Is.positive.y = TRUE, 
                         ncol.w.max = Inf, 
                         ncol.y.max = Inf, 
                         out.wy = TRUE, 
                         colsyperw = 1, 
                         maximize = TRUE)
      
      w <- check$w
      y <- check$y
      
      M     <- M1 * ncol(y)
      ncoly <- ncol(y)
      NOS   <- ncoly
      
      if ( length( .iscale ) && (length( .iscale ) > NOS) )
        stop("Invalid number of initial values for the scale parameter.")
      if ( length( .ishape ) && (length( .ishape ) > NOS) )
        stop("Invalid number of initial values for the shape parameter.")
      
      if ( length( .iscale ) )
        myIsc <- rep( .iscale, length = NOS )[1:NOS]
      if ( length( .ishape ) )
        myIsh <- rep( .ishape, length = NOS )[1:NOS]
      
      sca.names <- 
        if (NOS == 1) "scale" else paste("scale", 1:NOS, sep = "")
      sha.names <- 
        if (NOS == 1) "shape" else paste("shape", 1:NOS, sep = "")
      
      if ( .lss ) {
        # 20170226. Allow invweibullMeanlink to adequately work.
        parameter.names <- c(sca.names, sha.names)
        
        predictors.names <- 
          c(namesof(sca.names, .lscale , earg = .escale , tag = FALSE), 
            namesof(sha.names, .lshape , earg = .eshape , tag = FALSE))
      } else {
        # 20170226. Allow invweibullMeanlink to adequately work.
        parameter.names <- c(sha.names, sca.names)
        
        predictors.names <- 
          c(namesof(sha.names, .lshape , earg = .eshape , tag = FALSE), 
            namesof(sca.names, .lscale , earg = .escale , tag = FALSE))
      }
      
      parameter.names <- 
        parameter.names[interleave.VGAM( M1 * NOS , M1 = M1)]
      
      predictors.names <- 
        predictors.names[interleave.VGAM( M1 * NOS , M1 = M1)]
      
      if (!length(etastart)) {
        init.scale <- matrix(if (length( .iscale )) myIsc else 0.0 + NA, 
                             n, NOS, byrow = TRUE)
        init.shape <- matrix(if (length( .ishape )) myIsh else 0.0 + NA, 
                             n, NOS, byrow = TRUE)
        
        scale.grid <- .gscale
        shape.grid <- .gshape
        
        if ( .imethod == 1 ) 
          
          for (k in 1:ncoly) {
            yvec <- y[, k]
            wvec <- w[, k]
            invWeiLoglink <- function(scaleval, y, x, w, extraargs) {
              ans <- sum (c(w) * dinvweibull(x = y,
                                             scale = scaleval, 
                                             shape = extraargs$Shape, 
                                             log = TRUE))
              ans
            }
            
            localMatrix <- matrix(-1.0, length(shape.grid), 3)
            for (jj in 1:length(shape.grid)) {
              localMatrix[jj, 1] <- shape.grid[jj]
              localMatrix[jj, 2:3] <- 
                grid.search(scale.grid, objfun = invWeiLoglink, 
                            y = yvec, x = x, w = wvec,
                            ret.objfun = TRUE,
                            extraargs = list(Shape = shape.grid[jj]))
            }
            
            index.Maxshape <- 
              which(localMatrix[, 3] == max(localMatrix[, 3]))[1]
            if ( !length( .iscale ) )
              init.scale[, k] <- localMatrix[index.Maxshape, 2]
            if ( !length( .ishape ) )
              init.shape[, k] <- localMatrix[index.Maxshape, 1]
          }
       
        if ( .imethod == 2 ) 
          for (j in 1:NOS) {
            probs.y <- .probs.y
            xvec    <- log(-log(probs.y))
            yvec    <- log(quantile(y[, j], probs = probs.y))
            locfit  <- lsfit(x = xvec, y = yvec)
            if ( !length( .iscale ) )
              init.scale[, j] <- exp(locfit$coefficients["Intercept"])
            if ( !length( .ishape ) ) 
              init.shape[, j] <- -1/locfit$coefficients["X"]
          }    
        
        if ( .lss ) 
          etastart <- 
            cbind(theta2eta(init.scale, .lscale , earg = .escale), 
                  theta2eta(init.shape, .lshape , earg = .eshape)) else 
          etastart <- 
            cbind(theta2eta(init.shape, .lshape , earg = .eshape), 
                  theta2eta(init.scale, .lscale , earg = .escale))
      }
      
      etastart <- 
        etastart[, interleave.VGAM(M1 * NOS , M1 = M1), drop = FALSE]
      
    }), list( .lscale = lscale, .lshape = lshape, 
              .escale = escale, .eshape = eshape,
              .iscale = iscale, .ishape = ishape, 
              .gshape = gshape, .gscale = gscale,
              .imethod = imethod, .probs.y = probs.y, 
              .lss = lss ))),  
   
   linkinv = eval(substitute(function(eta, extra = NULL) {
     
     M1 <- 2
     NOS <- ncol(eta)/M1
     
     if ( .lss ) {
       iwscale <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                            .lscale , earg = .escale)
       iwshape <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                            .lshape , earg = .eshape)
     } else {
       iwshape <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                            .lshape , earg = .eshape)
       iwscale <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                            .lscale , earg = .escale)
     }
     ans <- cbind(iwscale * gamma(1 - 1/iwshape))
     
     ans
     
   }, list( .lscale = lscale , .lshape = lshape,
            .escale = escale , .eshape = eshape, 
            .lss = lss ))),
   
   last = eval(substitute(expression({ 
      
      if ( .lss ) {
        iwscale <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                             .lscale , earg = .escale)
        iwshape <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                             .lshape , earg = .eshape)
      } else {
        iwshape <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                             .lshape , earg = .eshape)
        iwscale <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                             .lscale , earg = .escale)
      }
      
      if (any(iwshape < 2)) 
        warning("MLE regularity conditions are violated (shape < 2) ",
                "at the ", "\n", "final iteration: Either MLE's are ",
                "not consistent, or MLE's ", "\n", "exist but are not ", 
                "asymptotically normal.")
      if (any( abs(iwshape - 2) < 1e-2  )) 
        warning("MLE regularity conditions are violated (shape == 2) ",
                "at the ", "\n", "final iteration: MLEs exist and are ",
                "normal and asymptotically ", "\n", "efficient but with ",
                "a slower convergence rate than when shape > 2.")
      
      misc$M1  <- M1
      misc$lss <- .lss 
      misc$expected <- TRUE
      misc$multipleResponses <- TRUE
      
      mynames <- if ( .lss )
        c(sca.names, sha.names) else
          c(sha.names, sca.names)
      mynames <- mynames[interleave.VGAM( M1 * NOS, M1 = M1 )]
      
      auxvec <- c( rep( if ( .lss ) .lscale else .lshape, length = NOS ),
                   rep( if ( .lss ) .lshape else .lscale, length = NOS))
      misc$link <- auxvec[ interleave.VGAM( M1 * NOS, M1 = M1) ]
      names(misc$link) <- mynames
      
      misc$earg <- vector("list", length = M )
      for (ii in 1:NOS) {
        misc$earg[[ M1 * ii - 1]] <- if ( .lss ) .escale else .eshape
        misc$earg[[ M1 * ii  ]]   <- if ( .lss ) .eshape else .escale
      } 
      names(misc$earg) <- mynames
      
    }), list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape ,
              .lss = lss ))),
   
   loglikelihood = eval(substitute(
     function(mu, y, w, residuals = FALSE, eta, 
              extra = NULL, summation = TRUE) {
       
       M1  <- 2
       NOS <- ncol(eta)/M1
       
       if( .lss ) {
         iwscale <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                              .lscale , earg = .escale)
         iwshape <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                              .lshape , earg = .eshape)
       } else {
         iwshape <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                              .lshape , earg = .eshape)
         iwscale <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                              .lscale , earg = .escale)
       }
       
       if (residuals) {
         stop("Loglikelihood not implemented yet",
              "to handle residuals.")
       } else {
         loglik.terms <- c(w)*dinvweibull(x = y, 
                                          scale = iwscale,
                                          shape = iwshape,
                                          log = TRUE)
         if (summation) 
           sum(loglik.terms)
         else 
           loglik.terms
       } 
       
     }, list( .lscale = lscale, .lshape = lshape,
              .escale = escale, .eshape = eshape,
              .lss = lss ))),
   
   vfamily = c("invweibull2mr"),
   
   deriv = eval(substitute(expression({
     
     M1 <- 2
     NOS <- ncol(eta)/M1
     if ( .lss ) {
       iwscale <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                            .lscale , earg = .escale)
       iwshape <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                            .lshape , earg = .eshape)
     } else {
       iwshape <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                            .lshape , earg = .eshape)
       iwscale <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                            .lscale , earg = .escale)
     }
     
     dl.dscale <- (iwshape/iwscale)*(1 - (y/iwscale)^(-iwshape))
     dl.dshape <- 1/iwshape - log(y/iwscale) + log(y/iwscale) * 
       (y/iwscale)^(-iwshape)
     
     dscale.deta <- dtheta.deta(iwscale, .lscale, earg = .escale)
     dshape.deta <- dtheta.deta(iwshape, .lshape, earg = .eshape)
     
     if ( .lss ) 
       myderiv <- c(w) * cbind(dl.dscale*dscale.deta, 
                               dl.dshape*dshape.deta) else 
              myderiv <- c(w) * cbind(dl.dshape*dshape.deta,
                                      dl.dscale*dscale.deta)
     
     myderiv[, interleave.VGAM(M1 * NOS , M1 = M1)]
     
   }), list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape,
             .lss = lss ))),
   
   weight = eval(substitute(expression({
     
     d1gamma.ev2 <- gamma(2) * digamma(2)
     d2gamma.ev2 <- (trigamma(2) * gamma(2)^2 + d1gamma.ev2^2)/gamma(2)
     
     ned2l.dscale2 <- (iwshape/iwscale)^2
     ned2l.dshape2 <- (1 + d2gamma.ev2)/(iwshape^2)
     ned2l.dscaledshape <- d1gamma.ev2/iwscale
     
     wz <- matrix(as.numeric(0.0), n, M + (M - 1))
     if (.lss ) {
       wz[, M1*(1:NOS) - 1]    <- ned2l.dscale2*(dscale.deta^2)
       wz[, M1*(1:NOS)]        <- ned2l.dshape2*(dshape.deta^2)
       wz[, M1*(1:NOS)+(M -1)] <- ned2l.dscaledshape * 
         (dscale.deta)*(dshape.deta)
     } else {
       wz[, M1*(1:NOS) - 1]    <- ned2l.dshape2*(dshape.deta^2)              
       wz[, M1*(1:NOS)]        <- ned2l.dscale2*(dscale.deta^2)
       wz[, M1*(1:NOS)+(M -1)] <- ned2l.dscaledshape * 
         (dscale.deta)*(dshape.deta)
     }
     
     w.wz.merge(w = w, wz = wz, n = n, M = 2*M - 1, ndepy = NOS)
     
   }), list( .lscale = lscale, .lshape = lshape,
             .escale = escale, .eshape = eshape, 
             .lss = lss ))),  
   
   simslot = eval(substitute(function(object, nsim) {
     
     pwts <- 
       if (length(pwts <- object@prior.weights) > 0) pwts else
         weights(object, type = "prior")
     
     if (any(pwts != 1))  warning ("Ignoring prior weights")
     
     eta <- predict(object)
     if ( .lss ) {
       iwscale <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                            .lscale , earg = .escale)
       iwscale <- eta2theta(eta[, M1*(1:NOS), drop = FALSE], 
                            .lshape , earg = .eshape)
     } else {
       iwscale <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                            .lshape , earg = .eshape)
       iwscale <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                            .lscale , earg = .escale)
     }
     rinvweibull(nsim*length(iwscale),
                 scale = iwscale,
                 shape = iwshape)
     
   }, list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape, .lss = lss )))
   
   ) # 'alo'
  
  alo
  
} 
