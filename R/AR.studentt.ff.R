##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.

AR.studentt.ff <- function(order = 1,
                           zero = c("scale", "df"),
                           llocation = "identitylink",
                           lscale    = "loglink",
                           ldf       = "logloglink",
                           ilocation = NULL,
                           iscale = NULL, idf = NULL,
                           imethod = 1) {
  
  
  
  lloc <- as.list(substitute(llocation))
  eloc <- link2list(lloc)
  lloc <- attr(eloc, "function.name")
  
  lsca <- as.list(substitute(lscale))
  esca <- link2list(lsca)
  lsca <- attr(esca, "function.name")
  
  ldof <- as.list(substitute(ldf))
  edof <- link2list(ldof)
  ldof <- attr(edof, "function.name")
  
  
  iloc <- ilocation
  isca <- iscale
  idof <- idf
  ord.int <- order; rm(order)
  
  
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  
  if (length(iloc))
    if (!is.Numeric(iloc))
      stop("bad input in argument 'ilocation'")
  if (length(isca))
    if (!is.Numeric(isca, positive = TRUE))
      stop("argument 'iscale' should be positive")
  if (length(idof))
    if (!is.Numeric(idof) || any(idof <= 1))
      stop("argument 'idf' should be > 1")
  
  
  
  new("vglmff",
      blurb = c("Order--",ord.int," AR model with Student-t errors\n\n",
                "Link:     ",
                namesof("location", lloc, earg = eloc), ", ",
                namesof("scale",    lsca, earg = esca), ", ",
                namesof("df",       ldof, earg = edof), "\n",
                "Variance: scale^2 * df / (df - 2) if df > 2\n"),
      constraints = eval(substitute(expression({
        
        constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                    predictors.names = predictors.names,
                                    M1 = 3)
      }), list( .zero = zero ))),
      infos = eval(substitute(function(...) {
        list(M1 = 3,
             Q1 = 1,
             expected = TRUE,
             multipleResponses = TRUE,
             parameters.names = c("location", "scale", "df"),
             zero = .zero)
      }, list( .zero = zero ))),
      
      
      first = eval(substitute(expression({
        
        y.int <- extra$y.int <- cbind(y)
        
        if (NCOL(y) > 1)
          stop("Currently, only univariate TS handled.")
        
        if ( ( .ord.int ) ) {
          mat.int <- NULL
          for (jj in 1:NCOL(y.int)) {
            temp1 <-  WN.lags(cbind(y.int[, jj]),  .ord.int )
            colnames(temp1) <- paste(colnames(y.int[, jj, drop = FALSE]),
                                     "Lag", 1:( .ord.int ), sep = "")
            mat.int <- cbind(mat.int, temp1)
          }
          extra$coeff <- coef(lm(y ~ mat.int))
          names(extra$coeff) <- c(c("(Intercept)"), 
                                     paste("AR", 1:( .ord.int ), sep = ""))
          mat.int <- residuals(lm(y ~ mat.int))
          y <- mat.int
        }
        
        if (FALSE) {
          x.matrix <- cbind(x[, 1], mat.int)
          colnames(x.matrix) <- c("(Intercept)", "Stut-residuals")
          list.names <- vector("list", NCOL(x.matrix) )
          names(list.names) <- colnames(x.matrix)
          for (ii in 1:NCOL(x.matrix)) 
            list.names[[ii]] <- ii
          attr(x.matrix, "assign") <- list.names
          #x <- x.matrix
          #print(head(mat.int))
        }
        
      }), list( .ord.int  = ord.int ))),
      
      
      initialize = eval(substitute(expression({
        M1 <- 3
        
        
        temp5 <-
          w.y.check(w = w, y = y,
                    ncol.w.max = Inf,
                    ncol.y.max = Inf,
                    out.wy = TRUE,
                    maximize = TRUE)
        w <- temp5$w
        y <- temp5$y
        
        extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
        extra$M1 <- M1
        M <- M1 * ncoly #
        
        mynames1 <- param.names("location", NOS, skip1 = TRUE)
        mynames2 <- param.names("scale",    NOS, skip1 = TRUE)
        mynames3 <- param.names("df",       NOS, skip1 = TRUE)
        predictors.names <-
          c(namesof(mynames1, .lloc , earg = .eloc , tag = FALSE),
            namesof(mynames2, .lsca , earg = .esca , tag = FALSE),
            namesof(mynames3, .ldof , earg = .edof , tag = FALSE))
        predictors.names <-
          predictors.names[interleave.VGAM(M1 * NOS, M1 = M1)]
        
        if (!length(etastart)) {
          init.loc <- if (length( .iloc )) .iloc else {
            if ( .imethod == 2) apply(y, 2, median) else
              if ( .imethod == 3) (colMeans(y) + t(y)) / 2 else {
                colSums(w * y) / colSums(w)
              }
          }
          
          sdvec <- apply(y, 2, sd)
          init.sca <- if (length( .isca )) .isca else
            sdvec / 2.3
          
          sdvec    <- rep_len(sdvec,    max(length(sdvec), length(init.sca)))
          init.sca <- rep_len(init.sca, max(length(sdvec), length(init.sca)))
          ind9 <- (sdvec / init.sca <= (1 + 0.12))
          sdvec[ind9] <- sqrt(1.12) * init.sca[ind9]
          init.dof <- if (length( .idof )) .idof else
            (2 * (sdvec / init.sca)^2) / ((sdvec / init.sca)^2  - 1)
          if (!is.Numeric(init.dof) || init.dof <= 1)
            init.dof <- rep_len(3, ncoly)
          
          mat1 <- matrix(theta2eta(init.loc, .lloc , earg = .eloc ), n, NOS,
                         byrow = TRUE)
          mat2 <- matrix(theta2eta(init.sca, .lsca , earg = .esca ), n, NOS,
                         byrow = TRUE)
          mat3 <- matrix(theta2eta(init.dof, .ldof , earg = .edof ), n, NOS,
                         byrow = TRUE)
          etastart <- cbind(mat1, mat2, mat3)
          etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
        }
      }), list( .lloc = lloc, .eloc = eloc, .iloc = iloc,
                .lsca = lsca, .esca = esca, .isca = isca,
                .ldof = ldof, .edof = edof, .idof = idof,
                .imethod = imethod ))),
      linkinv = eval(substitute(function(eta, extra = NULL) {
        NOS    <- extra$NOS
        y.int  <- extra$y.int
        M1 <- extra$M1
        Loc <-  eta2theta(eta[, M1*(1:NOS)-2], .lloc , earg = .eloc )
        Dof <-  eta2theta(eta[, M1*(1:NOS)-0], .ldof , earg = .edof )
        Loc[Dof <= 1] <- NA
        #print(head(cbind(extra$y.int, Loc)))
        Loc  + y.int
        
      }, list( .lloc = lloc, .eloc = eloc,
               .lsca = lsca, .esca = esca,
               .ldof = ldof, .edof = edof ))),
      last = eval(substitute(expression({
        M1 <- extra$M1
        misc$link <- c(rep_len( .lloc , NOS),
                       rep_len( .lsca , NOS),
                       rep_len( .ldof , NOS))
        misc$link <- misc$link[interleave.VGAM(M1 * NOS, M1 = M1)]
        temp.names <- c(mynames1, mynames2, mynames3)
        temp.names <- temp.names[interleave.VGAM(M1 * NOS, M1 = M1)]
        names(misc$link) <- temp.names
        
        misc$earg <- vector("list", M1 * NOS)
        names(misc$earg) <- temp.names
        for (ii in 1:NOS) {
          misc$earg[[M1*ii-2]] <- .eloc
          misc$earg[[M1*ii-1]] <- .esca
          misc$earg[[M1*ii  ]] <- .edof
        }
        
        misc$M1 <- M1
        misc$imethod <- .imethod
        misc$expected <- TRUE
        misc$multipleResponses <- TRUE
        cat("\nEstimated coefficients of the AR component:\n")
        print(extra$coef)
        cat("\n")
      }), list( .lloc = lloc, .eloc = eloc,
                .lsca = lsca, .esca = esca,
                .ldof = ldof, .edof = edof,
                .imethod = imethod ))),
      loglikelihood = eval(substitute(
        function(mu, y, w, residuals = FALSE, eta,
                 extra = NULL,
                 summation = TRUE) {
          NOS <- extra$NOS
          M1 <- extra$M1
          Loc <-  eta2theta(eta[, M1*(1:NOS)-2], .lloc , earg = .eloc )
          Sca <-  eta2theta(eta[, M1*(1:NOS)-1], .lsca , earg = .esca )
          Dof <-  eta2theta(eta[, M1*(1:NOS)-0], .ldof , earg = .edof )
          zedd <- (y - Loc) / Sca
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            ll.elts <- c(w) * (dt(x = zedd, df = Dof, log = TRUE) - log(Sca))
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
          }
        }, list(  .lloc = lloc, .eloc = eloc,
                  .lsca = lsca, .esca = esca,
                  .ldof = ldof, .edof = edof ))),
      vfamily = c("studentt3"),
      validparams = eval(substitute(function(eta, y, extra = NULL) {
        M1 <- extra$M1
        NOS <- extra$NOS
        Loc <- eta2theta(eta[, M1*(1:NOS)-2], .lloc , earg = .eloc )
        Sca <- eta2theta(eta[, M1*(1:NOS)-1], .lsca , earg = .esca )
        Dof <- eta2theta(eta[, M1*(1:NOS)-0], .ldof , earg = .edof )
        okay1 <- all(is.finite(Loc)) &&
          all(is.finite(Sca)) && all(0 < Sca) &&
          all(is.finite(Dof)) && all(0 < Dof)
        okay1
      }, list(  .lloc = lloc, .eloc = eloc,
                .lsca = lsca, .esca = esca,
                .ldof = ldof, .edof = edof ))),
      
      
      
      
      
      
      simslot = eval(substitute(
        function(object, nsim) {
          
          pwts <- if (length(pwts <- object@prior.weights) > 0)
            pwts else weights(object, type = "prior")
          if (any(pwts != 1))
            warning("ignoring prior weights")
          eta <- predict(object)
          Loc <-  eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lloc , earg = .eloc )
          Sca <-  eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lsca , earg = .esca )
          Dof <-  eta2theta(eta[, c(FALSE, FALSE, TRUE)], .ldof , earg = .edof )
          
          Loc + Sca * rt(nsim * length(Dof), df = Dof)
        }, list(  .lloc = lloc, .eloc = eloc,
                  .lsca = lsca, .esca = esca,
                  .ldof = ldof, .edof = edof ))),
      
      
      
      
      
      deriv = eval(substitute(expression({
        M1 <- extra$M1
        NOS <- extra$NOS
        Loc <- eta2theta(eta[, M1*(1:NOS)-2], .lloc , earg = .eloc )
        Sca <- eta2theta(eta[, M1*(1:NOS)-1], .lsca , earg = .esca )
        Dof <- eta2theta(eta[, M1*(1:NOS)-0], .ldof , earg = .edof )
        
        dloc.deta <- cbind(dtheta.deta(theta = Loc, .lloc , earg = .eloc ))
        dsca.deta <- cbind(dtheta.deta(theta = Sca, .lsca , earg = .esca ))
        ddof.deta <- cbind(dtheta.deta(theta = Dof, .ldof , earg = .edof ))
        
        zedd  <- (y - Loc) / Sca
        temp0 <- 1 / Dof
        temp1 <- temp0 * zedd^2
        dl.dloc <- (Dof + 1) * zedd / (Sca * (Dof + zedd^2))
        dl.dsca <- zedd * dl.dloc - 1 / Sca
        dl.ddof <- 0.5 * (-temp0 - log1p(temp1) +
                            (Dof+1) * zedd^2 / (Dof^2 * (1 + temp1)) +
                            digamma((Dof+1)/2) - digamma(Dof/2))
        
        ans <- c(w) * cbind(dl.dloc * dloc.deta,
                            dl.dsca * dsca.deta,
                            dl.ddof * ddof.deta)
        ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
        ans
      }), list( .lloc = lloc, .eloc = eloc,
                .lsca = lsca, .esca = esca,
                .ldof = ldof, .edof = edof ))),
      weight = eval(substitute(expression({
        
        const1 <- (Dof + 1) / (Dof + 3)
        const2 <- (Dof + 0) / (Dof + 3)
        const1[!is.finite(Dof)] <- 1  # Handles Inf
        const2[!is.finite(Dof)] <- 1  # Handles Inf
        
        const4 <- dnorm(0)
        ned2l.dlocat2 <-      const1 / (Sca * (Kayfun.studentt(Dof) / const4))^2
        ned2l.dscale2 <- 2  * const2 /  Sca^2
        
        DDS  <- function(df)          digamma((df + 1) / 2) -  digamma(df/2)
        DDSp <- function(df)  0.5 * (trigamma((df + 1) / 2) - trigamma(df/2))
        
        
        tmp6 <- DDS(Dof)
        edl2.dnu2 <- 0.5 * (tmp6 * (const2 * tmp6 - 2 / (Dof + 1)) - DDSp(Dof))
        ned2l.dshape2 <- cbind(edl2.dnu2)  # cosmetic name change
        
        ned2l.dshape.dlocat <- cbind(0 * Sca)
        ned2l.dshape.dscale <- cbind((-1 / (Dof + 1) + const2 * DDS(Dof))/Sca)
        
        
        
        wz <- array(c(c(w) * ned2l.dlocat2 * dloc.deta^2,
                      c(w) * ned2l.dscale2 * dsca.deta^2,
                      c(w) * ned2l.dshape2 * ddof.deta^2,
                      c(w) * ned2l.dshape2 * 0,
                      c(w) * ned2l.dshape.dscale  * dsca.deta * ddof.deta,
                      c(w) * ned2l.dshape.dlocat * dloc.deta * ddof.deta),
                    dim = c(n, M / M1, 6))
        wz <- arwz2wz(wz, M = M, M1 = M1)
        
        
        
        if (FALSE) {
          wz <- matrix(0.0, n, dimm(M))
          wz[, M1*(1:NOS) - 2] <- ned2l.dlocat2 * dloc.deta^2
          wz[, M1*(1:NOS) - 1] <- ned2l.dscale2 * dsca.deta^2
          wz[, M1*(1:NOS) - 0] <- ned2l.dshape2 * ddof.deta^2
          
          for (ii in ((1:NOS) - 1)) {
            ind3 <- 1 + ii
            wz[, iam(ii*M1 + 1, ii*M1 + 3, M = M)] <-
              ned2l.dshape.dlocat[, ind3] *
              dloc.deta[, ind3] * ddof.deta[, ind3]
            wz[, iam(ii*M1 + 2, ii*M1 + 3, M = M)] <-
              ned2l.dshape.dscale[, ind3] *
              dsca.deta[, ind3] * ddof.deta[, ind3]
          }
          
          while (all(wz[, ncol(wz)] == 0))
            wz <- wz[, -ncol(wz)]
        }
        
        
        
        wz
      }), list( .lloc = lloc, .eloc = eloc,
                .lsca = lsca, .esca = esca,
                .ldof = ldof, .edof = edof ))))
}