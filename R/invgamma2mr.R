##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0

# Function 'invgamma2mr'.
# This ff was written based upon gamma2(), from Thomas W. Yee (2006).



invgamma2mr <-function(lmu      = "loglink", 
                       lshape   = logofflink(offset = -2), 
                       parallel = FALSE,  
                       ishape   = NULL, 
                       imethod  = 1, 
                       zero     = "shape") {
  
  if (is.logical(parallel) && parallel && length(zero))
    stop("If 'parallel = TRUE', then set 'zero = NULL'.")
  
  if (!is.Numeric(imethod, 
                  length.arg = 1, 
                  integer.valued = TRUE, 
                  positive = TRUE) || imethod > 2)
    stop("Bad input for argument 'imethod'. It must be either 1 or 2.")
  
  if (length(ishape) && !is.Numeric( ishape, positive = TRUE ))
    stop("Bad input for argument 'ishape'.")
  
  apply.parint <- FALSE
  
  lmu    <- as.list(substitute(lmu))
  emu    <- link2list(lmu)
  lmu    <- attr(emu, "function.name")
  
  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")
  lss    <-  TRUE # Might be an argument...
  
  alo <- 
    new("vglmff", 
        blurb = c("2-parameter Inverse Gamma Distribution. \n", 
                  "Links : ", 
                  namesof("mu"   , lmu  ,  earg = emu, tag = FALSE), 
                  ", ",  
                  namesof("shape", lshape, earg = eshape, tag = FALSE), 
                  "\n", 
                  "Mean:      mu  \n", 
                  "Variance:  ((mu)^2)/(shape - 2) , shape > 2."),
        
   constraints = eval(substitute(expression({   
      M1  <- 2
      constraints <- 
        cm.zero.VGAM(constraints, x = x, zero = .zero , M = M,
                     predictors.names = predictors.names,
                     M1 = M1)
      
    }), list ( .zero         = zero, 
               .parallel     = parallel,
               .apply.parint = apply.parint ))),
        
    infos = eval(substitute(function(...){
      
      list(M1 = 2, 
           Q1 = 1, 
           expected = TRUE, 
           multipleResponses = TRUE, 
           lmu = .lmu , 
           lshape = .lshape , 
           zero = .zero )
      
    }, list ( .lmu = lmu, .lshape = lshape, .zero = zero ))),
   
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
      ncoly <- NOS <- ncol(y)
      M     <- M1 * ncol(y)
      
      if ( length( .ishape ) && (length( .ishape ) > NOS)  )
        stop(" Invalid number of initial values for the shape parameter.")
      
      if( length( .ishape ) && (length( .ishape ) <= NOS) )
        myIsh <- rep( .ishape , length = NOS )[1:NOS]
      
      mu.names <- 
        if (NOS == 1) "mu"    else paste("mu"  , 1:NOS   , sep = "")
      sha.names <- 
        if (NOS == 1) "shape" else paste("shape", 1:NOS  , sep = "")
      
      predictors.names <- 
        c(namesof(mu.names , .lmu    , earg = .emu ,    tag = FALSE),
          namesof(sha.names, .lshape , earg = .eshape , tag = FALSE))
      
      predictors.names <-
             predictors.names[interleave.VGAM( M1 * ncol(y) , M1 = M1)]
      
      if (is.logical( .parallel )  &  .parallel  &  ncoly > 1)
        warning("The constraint matrices may not be correct with multiple",
                " responses.")
      
      if (!length(etastart)) {
        
        init.shape <- 
           matrix( if ( length( .ishape ) ) myIsh else 1.0, 
                   nrow = n, ncol = NOS, byrow = TRUE)
        mymu       <- matrix( 1.0 , n , NOS )
        
        if (.imethod == 1) 
          for (ii in 1:ncol(y)) 
            mymu[, ii]  <- weighted.mean(y[, ii], w = w[, ii])
          
        if ( .imethod == 2 ) 
          for (spp. in 1:ncol(y)) 
            mymu[, spp.] <- (0.5)*(y[, spp.] + matrix(mean(y[, spp.]),
                                                      nrow = nrow(y),
                                                      ncol = 1,
                                                      byrow = TRUE))  
        
        for (spp in 1:NOS) {
          notdsr <- lsfit(x, y[, spp], wt = w[, spp], intercept = FALSE)
          var.y.est <- sum( w[, spp] * (notdsr$resid^2)) / 
                                  (n - length(notdsr$coef))
          
          if ( !length( .ishape ) )
            init.shape[, spp] <- ((mymu[, spp]^2)/var.y.est) + 2
        }
        
        etastart <- cbind(theta2eta(mymu      , .lmu   , earg = .emu   ),
                          theta2eta(init.shape, .lshape, earg = .eshape))
        etastart <- 
          etastart[, interleave.VGAM( M1 * NOS, M1 = M1), drop = FALSE]
      }
    }), list ( .lmu = lmu , .lshape = lshape , 
               .emu = emu , .eshape = eshape , 
               .ishape = ishape , 
               .parallel = parallel , 
               .apply.parint = apply.parint , 
               .imethod = imethod ))),
   
    linkinv = eval(substitute(function(eta, extra = NULL) {
      
      M1  <- 2
      NOS <- ncol(eta)/M1
      
      eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                .lmu , earg = .emu )
      
    }, list( .lmu = lmu, .emu = emu ))),
    
   
    last = eval(substitute(expression({
      
      if (any(Shape < 2)) 
        warning("MLE regularity conditions are violated ",
                "(shape < 2) at the final ", "\n",
                "iteration: Either MLE's are not consistent ",
                "or MLE's exist but these ", "\n", 
                "are not asymptotically normal.")
      
      if (any( abs(Shape - 2) < 1e-2  )) 
        warning("MLE regularity conditions are violated (shape == 2) ",
                "at the final ", "\n", "iteration: In this case, MLEs ",
                "exist and are normal and asymptotically ", "\n", 
                "efficient but with a slower convergence rate than ", 
                "when shape > 2.")
      
      misc$M1  <- M1
      misc$lss <- .lss
      misc$expected <- TRUE
      misc$parallel <- .parallel
      misc$multipleResponses <- TRUE
      misc$apply.parint      <- .apply.parint
      
      mynames <- 
        c(mu.names, sha.names)[interleave.VGAM( M1 * NOS, M1 = M1 )]
      misc$link <- c( rep( .lmu , M ) )
      auxvec <- c(rep( .lmu , length = NOS), rep( .lshape , length = NOS))
      misc$link <- auxvec[ interleave.VGAM( M1 * NOS, M1 = M1) ]
      names(misc$link) <- mynames
      
      misc$earg <- vector("list", length = M )
      for (ii in 1:NOS) {
        misc$earg[[ M1 * ii - 1]] <- if ( .lss ) .emu else .eshape
        misc$earg[[ M1 * ii    ]] <- if ( .lss ) .eshape else .emu
      } 
      names(misc$earg) <- mynames
      
    }), list( .lmu = lmu, .lshape = lshape,
              .emu = emu, .eshape = eshape,
              .parallel = parallel, .lss = lss ,
              .apply.parint = apply.parint ))),
         
    loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, eta,
               extra = NULL, summation = TRUE) {
        
        M1       <- 2
        NOS      <- ncol(eta)/M1
        mymu     <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                              .lmu    , earg = .emu)
        Shape    <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                              .lshape , earg = .eshape )
        if (residuals) {
          stop("Loglikelihood residuals not implemented yet.")
        } else {
          ll.elts <- c(w)*dinvgamma(x = y, 
                                    scale = c(mymu*(Shape -1)), 
                                    shape = c(Shape),
                                    log = TRUE)
          if (summation) 
            sum(ll.elts) else 
              ll.elts
        }
        
      }, list ( .lmu = lmu, .lshape = lshape,
                .emu = emu, .eshape = eshape ))), 
    
    vfamily = c("invgamma2mr"),
    
    simslot = eval(substitute(function(object, nsim) {
      
      pwts <- if (length(pwts <- object@prior.weights) > 0) {
        pwts
      } else {
        weights(object, type = "prior")
      }
      if (any(pwts != 1)) warning("Ignoring prior weights")
      eta   <- predict(object)  # Predictions
      mymu  <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                         .lmu , earg = .emu )
      Shape <- eta2theta(eta[, M1*(1:NOS), drop = FALSE], 
                         .lshape, earg = .eshape )
      rinvgamma(nsim*length(shape), 
                scale = c(mymu*(Shape - 1)), 
                shape = c(Shape))
      
    }, list( .lmu = lmu, .lshape = lshape, 
             .emu = emu, .eshape = eshape ))), 
    
    deriv = eval(substitute(expression({
      
      M1  <- 2
      NOS <- ncol(eta)/M1  
      
      mymu  <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                         .lmu ,    earg = .emu )
      Shape <- eta2theta(eta[, M1*(1:NOS), drop = FALSE], 
                         .lshape , earg = .eshape )
      
      dl.dmu <- (Shape/mymu) - ((Shape - 1)/y)
      dl.dshape <- (Shape/(Shape-1)) + log(Shape - 1) + log(mymu) +
        - log(y) - digamma(Shape) -(mymu/y)
      
      dmu.deta    <- dtheta.deta(mymu , .lmu   , earg = .emu   )
      dshape.deta <- dtheta.deta(Shape, .lshape, earg = .eshape)
      
      myderiv <- c(w)*cbind( dl.dmu*dmu.deta, dl.dshape*dshape.deta )
      myderiv[, interleave.VGAM( M1 * NOS , M1 = M1)]
      
    }), list ( .lmu = lmu, .lshape = lshape, 
               .emu = emu, .eshape = eshape ))), 
    
    weight = eval(substitute(expression({
      
      M1  <- 2
      NOS <- ncol(eta)/M1
      
      ned2l.dmu2    <- Shape/(mymu^2)
      ned2l.dshape2 <- trigamma(Shape) - ((Shape - 2)/(Shape - 1)^2)
      ned2l.dmymudshape <- 1/(mymu*(Shape - 1))
      
      wz <- matrix(as.numeric(0.0), n, M + (M - 1)) 
      wz[, M1*(1:NOS) - 1   ] <- ned2l.dmu2 * (dmu.deta^2)
      wz[, M1*(1:NOS)       ] <- ned2l.dshape2 * (dshape.deta^2)
      wz[, M1*(1:NOS)+(M -1)] <- ned2l.dmymudshape * 
        (dmu.deta) * (dshape.deta)
      w.wz.merge(w = w, wz = wz, n = n, M = M + (M - 1), ndepy = NOS)
      
      
    }), list( .lmu = lmu, .lshape = lshape,
              .emu = emu, .eshape = eshape )))
    
    ) # End of alo
  
  alo
 } 
