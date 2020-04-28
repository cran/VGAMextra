##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

# This function is a 'updated version' of the original script written 
# by Thomas W. Yee, University of Auckland.
# Modified on: 20150123, 20150913 (Renamed as 'gen.BetaIImr')

gen.betaIImr    <- function(lscale    = "loglink", 
                            lshape1.a = "loglink", 
                            lshape2.p = "loglink", 
                            lshape3.q = "loglink", 
                            iscale    = NULL, 
                            ishape1.a = NULL, 
                            ishape2.p = NULL, 
                            ishape3.q = NULL,
                            imethod   = 1, 
                            lss       = TRUE,
                            gscale    = exp(-5:5),
                            gshape1.a = exp(-5:5),
                            gshape2.p = exp(-5:5),
                            gshape3.q = exp(-5:5),
                            probs.y   = c(0.25, 0.50, 0.75), 
                            zero      = "shape" ) {
  
  if (length(lss) && !is.logical(lss))
    stop("Argument 'lss' not specified correctly. ",
         "See online help for important information.")
  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, 
                  positive = TRUE) || imethod > 1)
    stop("Bad input for argument 'imethod'.")
  
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE)) 
    stop("Bad input for argument 'iscale'.")
  
  if (length(ishape1.a) && !is.Numeric(ishape1.a, positive = TRUE))
    stop("Bad input for argument 'ishape1.a'.")
  
  if (length(ishape2.p) && !is.Numeric(ishape2.p, positive = TRUE))
    stop("Bad input for argument 'ishape2.p'.")
  
  if (length(ishape3.q) && !is.Numeric(ishape3.q, positive = TRUE))
    stop("Bad input for argument 'ishape3.q'.")
  
  if (length(probs.y) < 2 || max(probs.y) > 1 || 
        !is.Numeric(probs.y, positive = TRUE))
    stop("Bad input for argument 'probs.y'.")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")
  
  lshape1.a <- as.list(substitute(lshape1.a))
  eshape1.a <- link2list(lshape1.a)
  lshape1.a <- attr(eshape1.a, "function.name")
  
  lshape2.p <- as.list(substitute(lshape2.p))
  eshape2.p <- link2list(lshape2.p)
  lshape2.p <- attr(eshape2.p, "function.name")
  
  lshape3.q <- as.list(substitute(lshape3.q))
  eshape3.q <- link2list(lshape3.q)
  lshape3.q <- attr(eshape3.q, "function.name")
  
  new("vglmff", 
      blurb = 
        c("Generalized Beta II distribution \n\n", 
          "Links:    ", 
          ifelse (lss, 
                  namesof("scale"   , lscale   , earg = escale), 
                  namesof("shape1.a", lshape1.a, earg = eshape1.a)), ", ", 
          ifelse (lss, 
                  namesof("shape1.a", lshape1.a, earg = eshape1.a), 
                  namesof("scale"   , lscale   , earg = escale)), ", ", 
          namesof("shape2.p" , lshape2.p, earg = eshape2.p), ", ", 
          namesof("shape3.q" , lshape3.q, earg = eshape3.q), "\n", 
          "Mean:     scale * gamma(shape2.p + 1/shape1.a) * ", 
          "gamma(shape3.q - 1/shape1.a) / ", 
          "(gamma(shape2.p) * gamma(shape3.q))"),
      
    constraints = eval(substitute(expression({
      M1  <- 4
      constraints <- 
        cm.zero.VGAM(constraints, x = x, zero = .zero , M = M,
                     predictors.names = predictors.names,
                     M1 = M1)
      
    }), list( .zero = zero ))),
    
    
    infos = eval(substitute(function(...) {
      
      list(M1 = 4, 
           Q1 = 1, 
           expected = TRUE, 
           multipleResponses = TRUE, 
           lscale = .lscale , 
           lshape1.a = .lshape1.a , 
           lshape2.p = .lshape2.p , 
           lshape3.q = .lshape3.q,
           .zero = zero )
      
    }, list( .lscale = lscale      , .lshape1.a = lshape1.a,
             .lshape1.a = lshape1.a, .lshape2.p = lshape2.p, 
             .lshape3.q = lshape3.q, .zero = zero ))),
    
    initialize = eval(substitute(expression({
      
      check <- w.y.check(w = w, y = y, 
                         Is.positive.y = TRUE, 
                         ncol.w.max = Inf, 
                         ncol.y.max = Inf, 
                         out.wy = TRUE, 
                         colsyperw = 1, 
                         maximize = TRUE)
      y    <- check$y
      w    <- check$w
      M1   <- 4   
      NOS  <- ncoly <- ncol(y)
      M    <- M1*ncol(y)
      
      scale.names <- 
        if (NOS == 1) "scale" else paste("scale", 1:NOS, sep = "")
      sha1.names <- 
        if (NOS == 1) "shape1.a" else paste("shape1.a", 1:NOS, sep = "")
      sha2.names <-
        if (NOS == 1) "shape2.p" else paste("shape2.p", 1:NOS, sep = "")
      sha3.names <-
        if (NOS == 1) "shape3.q" else paste("shape3.q", 1:NOS, sep = "")
      
      if ( length( .iscale ) && ( length( .iscale ) > NOS ) )
        stop("Conflicting number of initial values for the scale param.")
      if ( length( .ishape1.a ) && ( length( .ishape1.a ) > NOS ) )
        stop("Conflicting number of initial values for the shape1 param.")
      if ( length( .ishape2.p ) && ( length( .ishape2.p ) > NOS ) )
        stop("Conflicting number of initial values for the shape2 param.")
      if ( length( .ishape3.q ) && ( length( .ishape3.q ) > NOS ) )
        stop("Conflicting number of initial values for the shape3 param.")
      
      if ( length( .iscale ) )
        iniSca <- rep( .iscale , length = NOS)[1:NOS]
      if ( length( .ishape1.a ) )
        iniSh1 <- rep( .ishape1.a , length = NOS)[1:NOS]
      if ( length( .ishape2.p ) )
        iniSh2 <- rep( .ishape2.p , length = NOS)[1:NOS]
      if ( length( .ishape3.q ) )
        iniSh3 <- rep( .ishape3.q , length = NOS)[1:NOS]
      
      predictors.names <- 
        if ( .lss ) {
          c(namesof(scale.names, .lscale    , 
                    earg = .escale    , tag = FALSE), 
            namesof(sha1.names , .lshape1.a , 
                    earg = .eshape1.a , tag = FALSE),
            namesof(sha2.names , .lshape2.p , 
                    earg = .eshape2.p , tag = FALSE), 
            namesof(sha3.names , .lshape3.q , 
                    earg = .eshape3.q , tag = FALSE))
        } else {
          c(namesof(sha1.names , .lshape1.a , 
                    earg = .eshape1.a , tag = FALSE), 
            namesof(scale.names, .lscale    , 
                    earg = .escale    , tag = FALSE),
            namesof(sha2.names , .lshape2.p , 
                    earg = .eshape2.p , tag = FALSE), 
            namesof(sha3.names , .lshape3.q , 
                    earg = .eshape3.q , tag = FALSE))
        }
      predictors.names <- 
        predictors.names[interleave.VGAM( M1 *NOS, M1 = M1)]
      
      if (!length(etastart)) {
        sca.init <- matrix(if (length( .iscale   ))  iniSca   else 1.0, 
                           n, NOS, byrow = TRUE)
        aa       <- matrix(if (length( .ishape1.a )) iniSh1 else 1.0, 
                           n, NOS, byrow = TRUE)
        parg     <- matrix(if (length( .ishape2.p )) iniSh2 else 1.0, 
                           n, NOS, byrow = TRUE)
        qq       <- matrix(if (length( .ishape3.q )) iniSh3 else 1.0, 
                           n, NOS, byrow = TRUE)
        # grid.search (2 parameters) + lsfit (2 parameters) #
        if ( .imethod == 1 ) {
          for(kk in 1:NOS) {
            probs.y <- .probs.y
            ishape3.q <- if (length( .ishape3.q )) .ishape3.q else 1.0
            xvec <- log( (1-probs.y)^(-1/ ishape3.q ) - 1 )
            yvec <- quantile(y[, kk], probs = probs.y)
            fit0 <- lsfit(x = xvec, y = log(yvec))
            
            if (!length( .iscale )) 
              sca.init[, kk] <- exp(fit0$coef[1])
            if (!length( .ishape1.a ))
              aa[, kk]       <- abs(1 / fit0$coef[2])
          }
          
          gshape2.p  <- .gshape2.p
          gshape3.q  <- .gshape3.q
          
          for (spp. in 1:NOS) {
            
            yvec <- y[, spp.]
            wvec <- w[, spp.]
            mymat <- matrix(-1.0, nrow = length(gshape3.q), ncol = 2) 
            localfun <- function(shape2eval, x = x, y = y, w = w, 
                                 extraargs) {
              ans <- sum(c(w) * dgen.betaII(x = y,
                                            scale    = extraargs$Scale, 
                                            shape1.a = extraargs$Shape1.a,
                                            shape2.p = shape2eval, 
                                            shape3.q = extraargs$Shape3.q,
                                            log = TRUE))
              ans
            }
            # Grid.search applied to 'gshape3.q' and 'gshape2.p' #
            for(jj in 1:length(gshape3.q)) {
              mymat[jj, ] <- 
                grid.search(gshape2.p, objfun = localfun,
                            y = yvec, x = x, w = wvec,
                            ret.objfun = TRUE,
                            extraargs = list(Scale    = sca.init[jj, spp.],
                                             Shape1.a = aa[jj, spp.],
                                             Shape3.q = gshape3.q[jj]))
            }
            index.shamax <- which(mymat[, 2] == max(mymat[, 2]))[1]
            
            if (!length( .ishape2.p )) 
              parg[, spp.] <- mymat[index.shamax, 1]
            if (!length( .ishape3.q  ))
              qq[, spp.]   <- gshape3.q[index.shamax]
            
          } 
        } # end of .imethod == 1
        
        etastart <- if ( .lss ) 
          cbind(theta2eta(sca.init, .lscale    , earg = .escale    ), 
                theta2eta(aa      , .lshape1.a , earg = .eshape1.a ), 
                theta2eta(parg    , .lshape2.p , earg = .eshape2.p ), 
                theta2eta(qq      , .lshape3.q , earg = .eshape3.q)) else 
                  
              cbind(theta2eta(aa      , .lshape1.a , earg = .eshape1.a ), 
                    theta2eta(sca.init, .lscale    , earg = .escale    ),
                    theta2eta(parg    , .lshape2.p , earg = .eshape2.p ), 
                    theta2eta(qq      , .lshape3.q , earg = .eshape3.q))
      } # End of etastart
      
      etastart <- etastart[, interleave.VGAM(M1 * NOS, M1 = M1)]
      
    }), list( .lscale    = lscale    , .escale    = escale,
              .lshape1.a = lshape1.a , .eshape1.a = eshape1.a,
              .lshape2.p = lshape2.p , .eshape2.p = eshape2.p,
              .lshape3.q = lshape3.q , .eshape3.q = eshape3.q,
              .iscale    = iscale    , .ishape1.a = ishape1.a,
              .ishape2.p = ishape2.p , .ishape3.q = ishape3.q,
              .imethod   = imethod   , .probs.y = probs.y, .lss = lss, 
              .gscale    = gscale    , .gshape1.a = gshape1.a , 
              .gshape2.p = gshape2.p , .gshape3.q = gshape3.q ))),
    
    linkinv = eval(substitute(function(eta, extra = NULL) {
      M1 <- 4
      NOS <- ncol(eta)/M1
      
      if ( .lss ) {
        Scale <- eta2theta(eta[, M1*(1:NOS) - 3, drop = FALSE], 
                           .lscale ,  earg = .escale )
        aa    <- eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE], 
                           .lshape1.a , earg = .eshape1.a )
      } else {
        aa    <- eta2theta(eta[, M1*(1:NOS) - 3, drop = FALSE], 
                           .lshape1.a , earg = .eshape1.a )
        Scale <- eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE], 
                           .lscale    ,  earg = .escale )
      }
      parg <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                        .lshape2.p , earg = .eshape2.p )
      qq   <- eta2theta(eta[, M1*(1:NOS), drop = FALSE],
                        .lshape3.q , earg = .eshape3.q )
      
      ans <- cbind(Scale * exp(lgamma(parg + 1/aa) +
                      lgamma(qq   - 1/aa) - lgamma(parg) - lgamma(qq)))
      
      ans
      
    }, list( .lscale    = lscale   , .lshape1.a = lshape1.a,
             .lshape2.p = lshape2.p, .lshape3.q = lshape3.q, 
             .escale    = escale   , .eshape1.a = eshape1.a, 
             .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
             .lss = lss ))),
    
    last = eval(substitute(expression({
       
      mynames <- c( if ( .lss ) scale.names else sha1.names,
                    if ( .lss ) sha1.names else scale.names,
                    sha2.names, sha3.names)
      mynames <- mynames[interleave.VGAM( M1 * NOS, M1 = M1 )]
      
      auxvec <- c( rep( if ( .lss ) .lscale else .lshape1.a, len = NOS),
                   rep( if ( .lss ) .lshape1.a else .lscale , len = NOS),
                   rep( .lshape2.p, len = NOS ) ,
                   rep( .lshape3.q, len = NOS ))
      misc$link <- auxvec[ interleave.VGAM( M1 * NOS , M1 = M1)]
      names(misc$link) <- mynames
      
      misc$eaerg <- vector("list", length = M)
      for ( ii in 1:NOS ) {
        misc$earg[[ M1 * ii - 3]] <- if ( .lss ) .escale else .eshape1.a
        misc$earg[[ M1 * ii - 2]] <- if ( .lss ) .eshape1.a else .escale
        misc$earg[[ M1 * ii - 1]] <- .eshape2.p
        misc$earg[[ M1 * ii    ]] <- .eshape3.q
      }
      names(misc$earg) <- mynames
      
      misc$M1  <- M1
      misc$expected          <- TRUE
      misc$multipleResponses <- TRUE
      
    }), list( .lscale    = lscale   , .lshape1.a = lshape1.a,
              .lshape2.p = lshape2.p, .lshape3.q = lshape3.q, 
              .escale    = escale   , .eshape1.a = eshape1.a,
              .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
              .lss = lss ))), 
    
    loglikelihood = eval(substitute(
      function(mu, y, w, residuals = FALSE, 
               eta, extra = NULL, summation = TRUE) {
        
        M1  <- 4
        NOS <- ncol(eta)/M1
        if ( .lss ) {
          scale <- eta2theta(eta[, M1*(1:NOS) - 3, drop = FALSE], 
                             .lscale    , earg = .escale    )
          aa    <- eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE], 
                             .lshape1.a , earg = .eshape1.a )
        } else {
          aa    <- eta2theta(eta[, M1*(1:NOS) - 3, drop = FALSE], 
                             .lshape1.a , earg = .eshape1.a )
          scale <- eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE], 
                             .lscale    , earg = .escale    )
        }
        parg   <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                            .lshape2.p , earg = .eshape2.p )
        qq     <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                            .lshape3.q , earg = .eshape3.q )
        if (residuals) {
          stop("loglikelihood residuals not implemented yet")
        } else {
          ll.elts <-
            c(w) * dgen.betaII(x = y, 
                               scale    = scale, 
                               shape1.a = aa,
                               shape2.p = parg,
                               shape3.q = qq,
                               log      = TRUE)
          if (summation) 
            sum(ll.elts) else 
              ll.elts
        }
        
      }, list( .lscale    = lscale   , .lshape1.a = lshape1.a, 
               .lshape2.p = lshape2.p, .lshape3.q = lshape3.q,
               .escale    = escale   , .eshape1.a = eshape1.a, 
               .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
               .lss = lss ))),
    
    vfamily = c("gen.betaIImr"),
    
    deriv = eval(substitute(expression({
      M1  <- 4
      NOS <- ncol(eta)/M1
      if ( .lss ) { 
        scale <- eta2theta(eta[, M1*(1:NOS) - 3, drop = FALSE], 
                           .lscale    , earg = .escale   )
        aa    <- eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE], 
                           .lshape1.a , earg = .eshape1.a)
      } else {
        aa    <- eta2theta(eta[, M1*(1:NOS) - 3, drop = FALSE], 
                           .lshape1.a , earg = .eshape1.a)
        scale <- eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE], 
                           .lscale    , earg = .escale )
      }  
      parg <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE], 
                        .lshape2.p , earg = .eshape2.p)
      qq   <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE], 
                        .lshape3.q , earg = .eshape3.q)
      temp1  <- log(y/scale)
      temp2  <- (y/scale)^aa
      temp3  <- digamma(parg + qq)
      temp3a <- digamma(parg)
      temp3b <- digamma(qq)
      temp4  <- log1p(temp2)
      
      dl.dscale <- (aa/scale) * (-parg + (parg+qq) / (1+1/temp2))
      dl.da <- 1/aa + parg * temp1 - (parg+qq) * temp1 / (1+1/temp2)
      dl.dp <- aa * temp1 + temp3 - temp3a - temp4
      dl.dq <- temp3 - temp3b - temp4
      
      dscale.deta <- dtheta.deta(scale, .lscale ,    earg = .escale )
      da.deta     <- dtheta.deta(aa,    .lshape1.a , earg = .eshape1.a )
      dp.deta     <- dtheta.deta(parg,  .lshape2.p , earg = .eshape2.p )
      dq.deta     <- dtheta.deta(qq,    .lshape3.q , earg = .eshape3.q )
      
      if ( .lss ) {
        myderiv <- c(w) * cbind(dl.dscale * dscale.deta, 
                                dl.da * da.deta, 
                                dl.dp * dp.deta, 
                                dl.dq * dq.deta)
      } else {
        myderiv <- c(w) * cbind(dl.da * da.deta, 
                                dl.dscale * dscale.deta, 
                                dl.dp * dp.deta, 
                                dl.dq * dq.deta)
      }
      
      myderiv[, interleave.VGAM(M1 * NOS, M1 = M1)]
      
    }), list(  .lscale    = lscale   , .lshape1.a = lshape1.a, 
               .lshape2.p = lshape2.p, .lshape3.q = lshape3.q, 
               .escale    = escale   , .eshape1.a = eshape1.a, 
               .eshape2.p = eshape2.p, .eshape3.q = eshape3.q, 
               .lss = lss ))),
    
    weight = eval(substitute(expression({
      
      temp5  <- trigamma(parg + qq)
      temp5a <- trigamma(parg)
      temp5b <- trigamma(qq)
      
      ned2l.da <- (1 + parg + qq + parg * qq * (temp5a + temp5b + 
                  (temp3b - temp3a + (parg-qq)/(parg*qq))^2 - 
                  (parg^2 + qq^2) / (parg*qq)^2)) / (aa^2 * (1+parg+qq))
      ned2l.dscale  <- (aa^2) * parg * qq / (scale^2 * (1+parg+qq))
      ned2l.dp      <- temp5a - temp5
      ned2l.dq      <- temp5b - temp5
      ned2l.dascale <- (parg - qq - parg * qq * (temp3a -temp3b)) / 
                       (scale*(1 + parg+qq))
      ned2l.dap     <- -(qq   * (temp3a -temp3b) -1) / (aa*(parg+qq))
      ned2l.daq     <- -(parg * (temp3b -temp3a) -1) / (aa*(parg+qq))
      ned2l.dscalep <-  aa * qq   / (scale*(parg+qq))
      ned2l.dscaleq <- -aa * parg / (scale*(parg+qq))
      ned2l.dpq     <- -temp5
      wz <- matrix(as.numeric(0.0), n, M + (M - 1) + (M - 2) + (M -3))
      wz <- if ( .lss ) {
        array(c(ned2l.dscale * dscale.deta^2, 
                ned2l.da * da.deta^2, 
                ned2l.dp * dp.deta^2, 
                ned2l.dq * dq.deta^2, 
                ned2l.dascale * da.deta * dscale.deta,
                ned2l.dap * da.deta * dp.deta,
                ned2l.dpq * dp.deta * dq.deta,
                ned2l.dscalep * dscale.deta * dp.deta, 
                ned2l.daq * da.deta * dq.deta,
                ned2l.dscaleq * dscale.deta * dq.deta),
              dim = c(n, M/M1, 10))
      } else {
        array(c(ned2l.da * da.deta^2,
                ned2l.dscale * dscale.deta^2,
                ned2l.dp * dp.deta^2,
                ned2l.dq * dq.deta^2,
                ned2l.dascale * da.deta * dscale.deta,
                ned2l.dscalep * dscale.deta * dp.deta,
                ned2l.dpq * dp.deta * dq.deta,
                ned2l.dap * da.deta * dp.deta,
                ned2l.dscaleq * dscale.deta * dq.deta,
                ned2l.daq * da.deta * dq.deta),
              dim = c(n, M/M1, 10))  # could be M1 * (M1 + 2) / 2
      }
      wz <- c(w) * arwz2wz(wz, M = M, M1 = M1)
      
      wz
      
    }), list( .lscale    = lscale   , .lshape1.a = lshape1.a, 
              .lshape2.p = lshape2.p, .lshape3.q = lshape3.q, 
              .eshape1.a = eshape1.a, .escale = escale      , 
              .eshape2.p = eshape2.p, .eshape3.q = eshape3.q,
              .lss = lss ))) 
    ) #End of new object....
}
