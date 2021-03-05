##########################################################################
# These functions are
# Copyright (C) 2014-2021 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

 weibullRff <-
    function(link1 = c("loglink", "weibullMlink", 
                       "weibullQlink")[1],
             lshape = "loglink", percentile = 50,
             imu = NULL, iscale = NULL, ishape = NULL,
             lss = TRUE, nrfs = 1,
             probs.y = c(0.2, 0.5, 0.8),
             imethod = 1, zero = "shape") {

  ## zzz - Discuss with Thomas
  test.link <- c("loglink", "weibullMlink", "weibullQlink")
  lscale <- match.arg(link1, test.link)[1]; rm(link1)
  Reliab <- 0.5 # not needed yet
    
  if (lscale == "weibullRLifelink"){
    if (!length(Reliab))
      stop("Enter a valid value for 'Reliab', between 0 and 1.")
    
    if (length(Reliab) && (log(Reliab) >= 0))
      stop("Wrong input for 'Reliab. Must be beetween 0 and 1.")
  }
  
  if (any(percentile < 1))
  warning("Some 'percentiles' are < 1? Usually they lie between 10 and 95")
  
  if (lscale == "weibullQlink") {
    percentile <- sort(percentile)
    if (!is.Numeric(percentile, positive = TRUE) || 
        (any(percentile > 100)))
      stop("Wrong input for 'percentile'.")
  }
  
  
  
  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")
  
  if (lscale == "weibullQlink") {
    lscale <- as.list(substitute(lscale(percentile = percentile))) 
  } else {
    if (lscale == "weibullRLifelink") {
      lscale <- as.list(substitute(lscale(Reliab = Reliab))) 
    } else {
      lscale <- as.list(substitute(lscale))
    }
  }

  
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(probs.y, positive  = TRUE) ||
      length(probs.y) < 2 ||
      max(probs.y) >= 1)
    stop("bad input for argument 'probs.y'")


  if (!is.Numeric(nrfs, length.arg = 1) ||
      nrfs < 0 ||
      nrfs > 1)
    stop("bad input for argument 'nrfs'")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' values must be positive")

  if (length(imu))
    if (!is.Numeric(imu, positive = TRUE))
      stop("argument 'imu' values must be positive")

  scale.TF <- if (lss) c(TRUE, FALSE) else c(FALSE, TRUE)
  scale.12 <- if (lss) 1:2 else 2:1
  blurb.vec <- c(namesof("scale", lscale, earg = escale),
                 namesof("shape", lshape, earg = eshape))
  blurb.vec <- blurb.vec[scale.12]

  new("vglmff",
  blurb = c("Weibull distribution\n\n",
            "Links:    ",
            blurb.vec[1], ", ",
            blurb.vec[2], "\n",
            "Mean:     scale * gamma(1 + 1/shape)\n",
            "Variance: scale^2 * (gamma(1 + 2/shape) - ",
                      "gamma(1 + 1/shape)^2)"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                          bool = FALSE ,
                           constraints = constraints )
   #  zzz Discuss with Thomas/ must be parameters.names.
   #  Otherwise, grep() will detect 'shape' twice therefore no cols.
   #  predictors.names = parameters.names instead of parameters.names
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                 predictors.names = parameters.names,
                                 M1 = 2)
   }), list( .zero = zero,
             .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = if ( .lss )
           c("scale", "shape") else
           c("shape", "scale"),
         lss = .lss ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .scale.12 = scale.12, .scale.TF = scale.TF,
           .lscale = lscale ,
           .lshape = lshape ,
           .lss = lss
         ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    ncoly <- ncol(y)
    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly
    
    ## zzz 
    if ( NCOL(y) != length( .percentile )) {
      toret <- "Response and 'percentile' dimensions do not match."
      toret <- paste(toret, "Use Q.reg() to set up the response.")
      stop(toret)
    }

    if (is.SurvS4(y))
      stop("only uncensored observations are allowed; ",
           "don't use SurvS4()")

    if ( .lss ) {
      mynames1 <- param.names("scale", ncoly, skip1 = TRUE)
      mynames2 <- param.names("shape", ncoly, skip1 = TRUE)
      predictors.names <-
          c(namesof(mynames1, .lscale , earg = .escale , tag = FALSE),
            namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))

    } else {
      mynames1 <- param.names("shape", ncoly, skip1 = TRUE)
      mynames2 <- param.names("scale", ncoly, skip1 = TRUE)
      predictors.names <-
          c(namesof(mynames1, .lshape , earg = .eshape , tag = FALSE),
            namesof(mynames2, .lscale , earg = .escale , tag = FALSE))
    }
    
    ### zzz
    parameters.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    ##
    predictors.names <- predictors.names[
          interleave.VGAM(M, M1 = M1)]

    ## zzz
    if (length( .imu )) {
      .iscale  <-  .imu / gamma(1 + 1 / ( .ishape) ) 
    }
    ###
    
    Shape.init <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                         n, ncoly, byrow = TRUE)
    Scale.init <- matrix(if (length( .iscale )) .iscale else 0 + NA,
                         n, ncoly, byrow = TRUE)

    if (!length(etastart)) {
      if (!length( .ishape ) ||
          !length( .iscale )) {
        for (ilocal in 1:ncoly) {

          anyc <- FALSE  # extra$leftcensored | extra$rightcensored
          i11 <- if ( .imethod == 1) anyc else FALSE  # Can be all data
          probs.y <- .probs.y
          xvec <- log(-log1p(-probs.y))
          fit0 <- lsfit(x  = xvec,
                        y  = log(quantile(y[!i11, ilocal],
                                 probs = probs.y )))
          
          if (!is.Numeric(Shape.init[, ilocal]))
            Shape.init[, ilocal] <- 1 / fit0$coef["X"]
          if (!is.Numeric(Scale.init[, ilocal]))
            Scale.init[, ilocal] <- exp(fit0$coef["Intercept"])
        }  # ilocal
        
        
        ## zzz
        extra$y <- y
        newescale <- .escale
       # if ( ( .lscale == "weibullMlink") ||
      #       ( .lscale == "weibullRMedianlink") ||
      #       ( .lscale == "weibullFRFlink") ||
      #       ( .lscale == "weibullRLifelink") ||
      #       ( .lscale == "weibullQlink"))  {
       if (!( .lscale == "loglink" )) {
           newescale$shape <- Shape.init
          if  ( .lscale == "weibullFRFlink")
            newescale$ycons <- y
          if  ( .lscale == "weibullRLifelink")
            newescale$Reliab <- .Reliab
          if  ( .lscale == "weibullQlink")
            newescale$percentile <- .percentile
        } 
        ###

        etastart <- if ( .lss )
          cbind(theta2eta(Scale.init, .lscale , earg = newescale ),
                theta2eta(Shape.init, .lshape , earg = .eshape ))[,
                interleave.VGAM(M, M1 = M1)] else
          cbind(theta2eta(Shape.init, .lshape , earg = .eshape ),
                theta2eta(Scale.init, .lscale , earg = newescale ))[,
                interleave.VGAM(M, M1 = M1)]
      }
    }
    
    
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .iscale = iscale, .ishape = ishape, .imu = imu,
            .probs.y = probs.y, .Reliab = Reliab ,
            .percentile = percentile,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss,
            .imethod = imethod ) )),
 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Shape <- eta2theta(eta[, !( .scale.TF )], .lshape , earg = .eshape )
    
    ## zzz
 #   if ( ( .lscale == "weibullMlink") ||
#         ( .lscale == "weibullRMedianlink") ||
#         ( .lscale == "weibullFRFlink") ||
#         ( .lscale == "weibullRLifelink") ||
#         ( .lscale == "weibullQlink")) {

   if (!( .lscale == "loglink" )) {
     
      newescale <- .escale
      newescale$shape <- Shape
      if  ( .lscale == "weibullRLifelink")
        newescale$Reliab <- .Reliab
      
      if  ( .lscale == "weibullFRFlink")
        newescale$ycons <- extra$y
      
      if  ( .lscale == "weibullQlink")
        newescale$percentile <- .percentile
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = newescale)
      
    } else {
      
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = .escale )  
    }
    ###
    
   ## zzz
    if ( .lscale == "weibullRMedianlink") {
      toret <- Scale * (log(2))^(1 / Shape)
    } else {
      
      if ( .lscale == "weibullFRFlink") {
        toret <- (Shape * extra$y^(Shape - 1)) / Scale^Shape
      } else {
        
        if ( .lscale == "weibullRLifelink") {
          toret <- Scale * (-log( .Reliab ))^(1 / Shape)
        } else {
          toret <-  Scale * gamma(1 + 1 / Shape) 
        }
      }
    } 
    
    if ( .lscale == "weibullQlink") {
      NOS <- NCOL(eta) / 2
      myperc <- matrix( .percentile / 1e2, nrow(eta), NOS, byrow = TRUE)
      toret <- Scale * (-log(1 - myperc ))^(1 / Shape)
    }

    toret
    
  }, list( .lscale = lscale, .lshape = lshape,
           .percentile = percentile ,
           .escale = escale, .eshape = eshape, .Reliab = Reliab ,
           .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ) )),
 
 
 
  last = eval(substitute(expression({
    regnotok <- any(Shape <= 2)
    if (any(Shape <= 1)) {
      warning("MLE regularity conditions are violated",
              "(shape <= 1) at the final iteration: ",
              "MLEs are not consistent")
    } else if (any(1 < Shape & Shape < 2)) {
      warning("MLE regularity conditions are violated",
              "(1 < shape < 2) at the final iteration: ",
              "MLEs exist but are not asymptotically normal")
    } else if (any(2 == Shape)) {
      warning("MLE regularity conditions are violated",
              "(shape == 2) at the final iteration: ",
              "MLEs exist and are normal and asymptotically ",
              "efficient but with a slower convergence rate than when ",
              "shape > 2")
    }



    M1 <- extra$M1
    avector <- if ( .lss ) c(rep_len( .lscale , ncoly),
                             rep_len( .lshape , ncoly)) else
                           c(rep_len( .lshape , ncoly),
                             rep_len( .lscale , ncoly))
    misc$link <- avector[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    
    ### zzz
    #if ( ( .lscale == "weibullMlink") ||
    #     ( .lscale == "weibullRMedianlink") ||
    #     ( .lscale == "weibullFRFlink" ) ||
    #     ( .lscale == "weibullRLifelink") ||
    #     ( .lscale == "weibullQlink"))
  
   if (!( .lscale == "loglink" )) {
      newescale <- .escale
      
      if ( .lscale == "weibullFRFlink") {
        newescale$ycons <- extra$y
      }
      
      if ( .lscale == "weibullRLifelink") {
        newescale$Reliab <- .Reliab
      }
      
      if  ( .lscale == "weibullQlink")
        newescale$percentile <- .percentile
      
      newescale$shape  <- if ( .lss )
        eta2theta(coef(fit)[1], .lshape , earg = .eshape ) else
          eta2theta(coef(fit)[2], .lshape , earg = .eshape )
    }
    
    for (ii in 1:ncoly) {   ### zzz
      misc$earg[[M1*ii-1]] <- if ( .lss ) newescale else .eshape
      misc$earg[[M1*ii  ]] <- if ( .lss ) .eshape else newescale 
    }
    
    ## zzz
    #if (.lscale == "weibullRMeanlink") {
    #  fit$fittedvalues <- 
    #    weibullRMeanlink(theta = Scale , shape = Shape, deriv = 0,
    #                     inverse  = FALSE )
    #}
    ###
    if ( .lscale == "loglink" ) {
      cat("\nNOTE: Same as VGAM::weibullR() as link1 = 'loglink'.\n\n")
    }
    if ( .lscale == "weibullQlink" ) {
      cat("\nLink = 'weibullQlink', percentiles =", 
          noquote(paste( .percentile, collapse = ", ")),
          ".\n\n")
    }
    
    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE


    misc$nrfs <- .nrfs
    misc$RegCondOK <- !regnotok # Save this for later
    
    
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .imethod = imethod, .Reliab = Reliab ,
            .percentile = percentile ,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss,
            .nrfs = nrfs ) )),
 
 
  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Shape <- eta2theta(eta[, !( .scale.TF )], .lshape , earg = .eshape )
    
    ## zzz
    #if ( ( .lscale == "weibullMlink") ||
    #     ( .lscale == "weibullRMedianlink") ||
    #     ( .lscale == "weibullFRFlink") ||
    #     ( .lscale == "weibullRLifelink") ||
    #     ( .lscale == "weibullQlink")) {
    
    if (!( .lscale == "loglink" )) {
      
      newescale <- .escale
      newescale$shape <- Shape
      
      if ( .lscale == "weibullRLifelink") {
        newescale$Reliab <- .Reliab
      }
      if ( .lscale == "weibullFRFlink") {
        newescale$ycons <- extra$y
      }
      
      if  ( .lscale == "weibullQlink")
        newescale$percentile <- .percentile
      
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = newescale) 
      
    } else {
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = .escale ) 
    }
    ###

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dweibull(y, shape = Shape, scale = Scale, log = TRUE))
    }
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape, .Reliab = Reliab ,
           .percentile = percentile ,
           .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ) )),
 
  vfamily = c("weibullRff"),
 
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    
    ## zzz
    Shape <- eta2theta(eta[, !( .scale.TF )], .lshape , earg = .eshape )
    
   # if ( ( .lscale == "weibullMlink") ||
  #       ( .lscale == "weibullRMedianlink") ||
  #       ( .lscale == "weibullFRFlink") ||
  #       ( .lscale == "weibullRLifelink") ||
  #       ( .lscale == "weibullQlink")) {
      
    if (!( .lscale == "loglink" )) {
      
      newescale <- .escale
      newescale$shape <- Shape
      
      if ( .lscale == "weibullRLifelink") {
        newescale$Reliab <- .Reliab
      }
      
      if ( .lscale == "weibullFRFlink") {
        newescale$ycons <- extra$y
      }
      if  ( .lscale == "weibullQlink")
        newescale$percentile <- .percentile
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = newescale) 
      
    } else {
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = .escale ) 
    }
    ###

    okay1 <- all(is.finite(Scale)) && all(0 < Scale) &&
             all(is.finite(Shape)) && all(0 < Shape)
    okay1
  }, list( .lscale = lscale, .lshape = lshape,
           .percentile = percentile ,
           .escale = escale, .eshape = eshape, .Reliab = Reliab ,
           .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ) )),
  
 deriv = eval(substitute(expression({
    M1 <- 2
    
    ## zzz
    NOS <- extra$ncoly
    Shape <- eta2theta(eta[, !( .scale.TF )], .lshape , earg = .eshape )
    
    # zzz 
    #if ( ( .lscale == "weibullMlink") ||
    #     ( .lscale == "weibullRMedianlink") ||
    #     ( .lscale == "weibullFRFlink") ||
    #     ( .lscale == "weibullRLifelink") || 
    #     ( .lscale == "weibullQlink")) {
    
    if (!( .lscale == "loglink" )) { 
      newescale <- .escale
      newescale$shape <- Shape
    
      if ( .lscale == "weibullRLifelink") {
        newescale$Reliab <- .Reliab
      }
      
      if ( .lscale == "weibullFRFlink") {
        newescale$ycons <- extra$y
      }
      
      if  ( .lscale == "weibullQlink")
        newescale$percentile <- .percentile
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = newescale) 

      # dscale.dshape <- dscale.dshape
      if ( .lscale == "weibullMlink") 
        dscale.dshape <- -Scale * digamma(1 + 1/ Shape) * (- 1/ Shape^2)
      if ( .lscale == "weibullRMedianlink")
        dscale.dshape <- log(log(2)) * Scale / Shape^2
      if ( .lscale == "weibullFRFlink") 
        dscale.dshape <- -Shape * log(Scale) +
          (Scale /Shape^2) * (1 + Shape * log(extra$y))
      if ( .lscale == "weibullRLifelink")
        dscale.dshape <- log(-log( .Reliab )) * Scale / Shape^2
      if ( .lscale == "weibullQlink") {
        myperc <- matrix( .percentile / 1e2 , nrow(eta), NOS, byrow = TRUE)
        dscale.dshape <- log(-log( 1 - myperc )) * Scale / Shape^2
      }
 

      ## zzz CHANGED  dl.dshape
      dz.da <- ((y/Scale)^Shape) * log(y/Scale)
      dz.db <- -(Shape / Scale) *  ((y/Scale)^Shape)
      dl.dshape <- 1 / Shape + log(y / Scale) - dz.da -
                 (dz.db + Shape / Scale) * dscale.dshape
      dl.dscale <- (Shape / Scale) * (-1.0 + (y / Scale)^Shape) # SAME
     
      newescale$wrt.param <- 1
      newescale$deriv <- 1
      newescale$inverse <- TRUE
      
      dscale.deta1 <- dtheta.deta(Scale , .lscale , 
                        earg = newescale)[, 1:NOS, drop = FALSE]
      dshape.deta1 <- dtheta.deta(Shape, .lscale ,
                       earg = newescale )[, -(1:NOS), drop = FALSE]
      
      newescale$wrt.param <- 2
      dscale.deta2 <- dtheta.deta(Scale, .lscale ,
                       earg = newescale)[, 1:NOS, drop = FALSE]
      dshape.deta2 <- dtheta.deta(Shape, .lshape , earg = .eshape )
      
      myderiv <-  if ( .lss )
        c(w) * cbind(dl.dscale * dscale.deta1 +
                       dl.dshape * dshape.deta1,
                     dl.dscale * dscale.deta2 +
                       dl.dshape * dshape.deta2) else 
                         c(w) * cbind(dl.dscale * dscale.deta2 +
                                        dl.dshape * dshape.deta2,
                                      dl.dscale * dscale.deta1 +
                                        dl.dshape * dshape.deta1)
      ####### end of zzz
      
    } else {
      # Ordinary weibullR()
      Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = .escale ) 
      
      dl.dshape <- 1 / Shape + log(y / Scale) -
                    log(y / Scale) * (y / Scale)^Shape
      dl.dscale <- (Shape / Scale) * (-1.0 + (y / Scale)^Shape)
      
      dshape.deta <- dtheta.deta(Shape, .lshape, earg = .eshape )
      dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
      
      myderiv <- if ( .lss )
        c(w) * cbind(dl.dscale * dscale.deta,
                     dl.dshape * dshape.deta) else
                       c(w) * cbind(dl.dshape * dshape.deta,
                                    dl.dscale * dscale.deta)  
    }
    
     myderiv[, interleave.VGAM(M, M1 = M1)]

  }), list( .lscale = lscale, .lshape = lshape,
            .percentile = percentile ,
            .escale = escale, .eshape = eshape, .Reliab = Reliab ,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ) )),
  
 weight = eval(substitute(expression({
    
    ### zzz
  # if ( ( .lscale == "weibullMlink") ||
  #      ( .lscale == "weibullRMedianlink") ||
  #      ( .lscale == "weibullFRFlink") || 
  #      ( .lscale == "weibullRLifelink") ||
  #      ( .lscale == "weibullQlink")) {
   
   if (!( .lscale == "loglink" )) {   
      EulerM <- -digamma(1.0)
      
      ned2l.dscale2 <- (Shape / Scale)^2
      ned2l.dshape.part <- (6*(EulerM - 1)^2 + pi^2)/(6*Shape^2) # KK(2003)
     
      ## zzz
      ned2l.dshape2 <- ned2l.dshape.part - 
        (2 / Scale) * digamma(2) * dscale.dshape +
        ( (Shape / Scale) * dscale.dshape)^2
      ### zzz
      ned2l.dshapescale <- (EulerM-1) / Scale + 
                      dscale.dshape * (Shape / Scale)^2
      
     term1 <- ned2l.dshape2 * dshape.deta2^2 + 
       2 * ned2l.dshapescale  * dshape.deta2 * dscale.deta2 +
       ned2l.dscale2 * dscale.deta2^2
     
     term2 <- ned2l.dshape2 * dshape.deta1^2 + 
       2 * ned2l.dshapescale  * dshape.deta1 * dscale.deta1 +
       ned2l.dscale2 * dscale.deta1^2
     
     term3 <- ned2l.dshape2 * dshape.deta1 * dshape.deta2 +
       ned2l.dshapescale  * ( dshape.deta1 * dscale.deta2 +
                               dshape.deta2 * dscale.deta1) +
       ned2l.dscale2 * dscale.deta1 * dscale.deta2
     
     if ( .lss ) {
       wz <- array(c(c(w) * term2, c(w) * term1, c(w) * term3 ),
                   dim = c(n, M / M1, 3)) 
       
     } else {
       wz <- array(c(c(w) * term1, c(w) * term2, c(w) * term3 ),
                   dim = c(n, M / M1, 3)) 
     }
      ### end of zzz 
    } else {
      # Ordinary weibullR()
      
      EulerM <- -digamma(1.0)  
      ned2l.dscale <- (Shape / Scale)^2
      ned2l.dshape <- (6*(EulerM - 1)^2 + pi^2)/(6*Shape^2)  # KK (2003)
      ned2l.dshapescale <- (EulerM-1) / Scale
      
      wz <- if ( .lss )
        array(c(c(w) * ned2l.dscale * dscale.deta^2,
                c(w) * ned2l.dshape * dshape.deta^2,
                c(w) * ned2l.dshapescale * dscale.deta * dshape.deta),
              dim = c(n, M / M1, 3)) else
            array(c(c(w) * ned2l.dshape * dshape.deta^2,
                  c(w) * ned2l.dscale * dscale.deta^2,
                  c(w) * ned2l.dshapescale * dscale.deta * dshape.deta),
                dim = c(n, M / M1, 3))
    }
    
    
    wz <- arwz2wz(wz, M = M, M1 = M1)


    wz
  }), list( .lscale = lscale , .lshape = lshape ,
            .eshape = eshape, .nrfs = nrfs,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ))))
}  # weibullRff
