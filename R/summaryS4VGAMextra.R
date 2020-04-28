##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
##########################################################################

# S4 summary methods for vector generalized time series ff at VGAMextra.
# Methods here are included as extra output into the
# summary output for vglm's, i.e., summary(object of class 'vglm'.)
# Extra material: 
# A) Standard errors according to MLE asympt. theory.
# B) Checks on stationarity / invertibility. 
# Specific methods of function 'summaryvglmS4VGAM' for TSff must be defined

#setClass("ARff", contains = "vglmff")
#setClass("MAff", contains = "vglmff")
#setClass("ARMAff", contains = c("ARff", "MAff"))
##setClass("ARMA.GARCHff", contains = c("vglmff") )
#setClass("AR.GARCHff", contains = c("ARMAff"))
#setClass("vgltsmff", contains = c("AR.GARCHff"),
#         slots = c(typeTS = "character"))
#setClass("vgtsff", contains = c("vgltsff-class"))
#setClass("vgltsmff", contains = "vglmff", )


# Old: vgltsff-class
setClass("ARMAvgltsmff", contains = c("vglmff"))

setClass("vgltsmff", contains = c("ARMAvgltsmff"),
         slots = c(typeTS = "character"))


setMethod(summaryvglmS4VGAM, signature(VGAMff = "ARMAvgltsmff"),
          function(object, VGAMff = VGAMff,...){
            
      no.mean <- (object@misc$nomean)
      no.mean <- if (length(no.mean)) no.mean else TRUE
      m.data  <- object@y
      x.data  <- object@x
      NoS     <- object@misc$NOS
      nn      <- nrow(m.data)
      m.Coef  <- numeric(0)
      zero    <- object@misc$zero
      intOrd  <- object@misc$Order
      var.arg <- object@misc$var.arg
      m.links <- object@misc$link
      m.eargs <- object@misc$earg
      process <- object@misc$process
      res.obj <- object@misc$residuals  # 20170418. ARs is NULL.
      

      process <- match.arg(process,
                           c("AR", "MA", "ARMA", "ARIMA", "pureGARCH"))[1]
      bigFlag <- ( (object@misc$flagAR) &&  (object@misc$flagMA) ) 
      
      if (process == "ARIMA") {
        
        if (!object@misc$MAord) {
          process <- "AR"
          intOrd  <- rep(object@misc$ARord, NoS)
        }
          
        if (!object@misc$ARord) {
          process <- "MA"
          intOrd  <- rep(object@misc$MAord, NoS)
        }
        
        if (object@misc$ARord && object@misc$MAord) {
          process <- "ARMA"
          intOrd  <- c(object@misc$ARord, object@misc$MAord)
        }
      }
      
      m.predict <- object@predictors
      m.names   <- names(object@constraints)
      sd.err    <- matrix(0.0, nrow = NoS, ncol = length(m.names))
      if (length(m.names) > 1) {
        const  <- Coef(object, matrix = TRUE)
        m.coln <- grep("noise", colnames(m.predict))
        for (jk in 1:NoS) 
          sd.err[jk, ] <- const[, m.coln[jk]]
        flag11 <- (ncol(sd.err) > 1)
        colnames(sd.err) <- m.names
      } else {
        flag11 <- FALSE
      }
        
      ### Coeffs ###
      which.nC <- grep("coeff", colnames(m.predict))
      if (length(which.nC)) {
        
        predCoe <- m.predict[, c(which.nC), drop = FALSE]
        predlink <- m.links[which.nC]
        predearg <- m.eargs[which.nC]

        crossCh <- if (process == "ARMA" || process == "ARIMA")
                  NoS * sum(intOrd) else sum(intOrd)
        if (ncol(predCoe) != crossCh)
          warning("Conflicting number of parameters. Differs from ",
                  " the order entered.")
        
        m.Coef <- numeric(0)
        for (ii in 1:crossCh) {
          predCoe[, ii] <- eta2theta(predCoe[, ii],
                                     link = predlink[ii],
                                     earg = predearg[[ii]])
          m.Coef <- c(m.Coef, mean(predCoe[, ii]))
        }
      } else {
        m.Coef <- NULL
        if (!bigFlag)
          warning("NO coefficients found.")
      }      
      ### WN variances ###
      which.nV <- grep("noise", colnames(m.predict))
      if (length(which.nV)) {
        varSd  <- m.predict[, c(which.nV), drop = FALSE]
        m.link <- m.links[which.nV]
        m.earg <- m.eargs[which.nV]
        
        m.var <- numeric(0)  
        for(kk in 1:NoS) {
          varSd[, kk]  <- eta2theta(varSd[, kk],
                                    link = m.link[kk],
                                    earg = m.earg[[kk]])
          m.var <- c(m.var, mean(varSd[, kk])^(2 - var.arg))
        }
      } else {
        m.var <- NULL
        warning("No white noise variance found.")
      }
      ### Drift or mean ###
      if (!no.mean) {
        which.nD <- switch(object@misc$process,
                           "AR" = "drift",
                           "MA" = "Mean",
                           "ARMA"  = "drift.mean",
                           "ARIMA" = "drift.mean",
                           "pureGARCH" = "Mean")
        which.nD <- grep(which.nD, colnames(m.predict))
        driftP <- m.predict[, c(which.nD), drop = FALSE]
        m.link <- m.links[c(which.nD)]
        m.earg <- m.eargs[c(which.nD)]
        
        m.interc <- numeric(0)
        for (ii in 1:NoS) {
          driftP[, ii] <- eta2theta(driftP[, ii],
                                    link = m.link[ii],
                                    earg = m.earg[[ii]])
          
          m.interc <- c(m.interc, mean(driftP[, ii]))
        }
      } else {
        m.interc <- NULL
      }
      
      ###  Allocating the coefficients by row i.e. per response ###
      aux <- 0
      if (!bigFlag) {
        if ((process == "AR" || process == "MA" )) {
          m.thetas <- matrix(NA, 
                             nrow = NoS, 
                             ncol = max(intOrd) + !no.mean)
          
          for (kk in 1:NoS) {
            m.thetas[kk, 1:(intOrd[kk] + !no.mean)] <-
              c(m.Coef[(aux + 1):(aux + intOrd[kk])], m.interc[kk])
            aux <- aux + intOrd[kk]
          }
          
        } else {
          m.thetas <- matrix(NA, 
                             nrow = NoS, 
                             ncol = sum(intOrd) + !no.mean)
          m.thetas[, 1:sum(intOrd)] <- 
            matrix(m.Coef, ncol = sum(intOrd), byrow = TRUE)
          
          if (!no.mean)
            m.thetas[, sum(intOrd) + !no.mean] <- cbind(m.interc)
        }
      } else {
        m.thetas <- matrix( if (!no.mean) m.interc else 0.0, NoS,  1)
      }
      
      ### standard errors ###
      if (process == "AR" && !bigFlag) {
        Covs   <- matrix(NA,  nrow = NoS, ncol = max(intOrd))
        se.Mat <- matrix(NA,  nrow = NoS, 
                         ncol = max(intOrd) + !no.mean)
        for (spp. in 1:NoS) {
          aOrd    <- intOrd[spp. ]
          m.data[, spp. ] <- if (no.mean) m.data[, spp. ] else
                             m.data[, spp. ] - mean(m.data[, spp. ])
          
          x.Mat <- matrix(0.0, nrow = nn, ncol = aOrd + !no.mean)
          for (kk in 1:aOrd)
            x.Mat[(kk + 1):nn, kk] <- m.data[1:(nn - kk), spp. ]
          if (!no.mean) 
            x.Mat[, aOrd + 1] <- 1
          ToepMat <- t(x.Mat) %*% x.Mat /nn
          if (!no.mean)
          ToepMat[aOrd + 1, aOrd + 1] <- ToepMat[aOrd + 1, aOrd + 1] *
                            (1 - sum(m.thetas[spp. , 1:aOrd]))^(-2)
          if ( abs(det(ToepMat)) < 1e-2 )
            warning("Variance matrix is not invertible.") 
          se.Mat[spp. , 1:(aOrd + !no.mean)]  <- 
            round(sqrt(m.var[spp. ] * diag(solve(ToepMat))/nn), 3)
        }
      } else {
        se.Mat <- matrix(sqrt(m.var / nn), NoS, 1)
      }
      
      if (process == "MA") {
        Covs <- matrix(NA,  nrow = NoS, ncol = max(intOrd))
        u.ar <- matrix(NA_real_, nrow = nn , ncol = NoS)
        se.Mat <- matrix(NA,  nrow = NoS, 
                         ncol = max(intOrd) + !no.mean)
        resMat <- matrix(NA, nrow = nn, ncol = NoS)
        
        for (spp. in 1:NoS) {
          aOrd    <- intOrd[spp. ]
          m.data[, spp. ] <- if (no.mean) m.data[, spp. ] else
                                 m.data[, spp. ] - m.interc[spp. ] 
        #resMat[, spp.] <- WN.InitARMA(tsData = data.frame(m.data[,spp.]),
        #                              order = c(0, 0, aOrd),
        #                              whiteN = TRUE,
        #                              moreOrder = 1,
        #                              updateWN = TRUE)[[2]]
          resMat[, spp. ] <- res.obj[, spp. ]
          x.Mat <- matrix(0.0, nrow = nn, ncol = aOrd + !no.mean)
          for (kk in 1:aOrd)
            x.Mat[(kk + 1):nn, kk] <- m.data[1:(nn - kk), spp. ]
          if (!no.mean) 
            x.Mat[, aOrd + 1] <- 1
          ToepMat <- t(x.Mat) %*% x.Mat /nn

          if ( abs(det(ToepMat)) < 1e-2 )
            warning("Variace matrix is not invertible.") 
          
          se.Mat[spp. , 1:(aOrd + !no.mean)]  <- 
            round(sqrt( m.var[spp. ] * diag(solve(ToepMat)) / nn ), 3)
        }
      }
      
      if ((process == "ARMA") || (process == "ARIMA")) {  
        arOrd  <- object@misc$ARord
        maOrd  <- object@misc$MAord
        sOrd   <- sum(c(arOrd, maOrd)) 
        the.ar <- m.thetas[, 1:arOrd, drop = FALSE]
        the.ma <- m.thetas[, (arOrd + 1):(arOrd + maOrd), drop = FALSE]
        resMat <- matrix(NA, nrow = nn, ncol = NoS)
        se.Mat <- matrix(NA, nrow = NoS, ncol = sOrd + !no.mean)
        
        for (spp. in 1:NoS)  {
          m.data[, spp. ] <- if (no.mean) m.data[, spp. ] else
                                  m.data[, spp. ] - mean(m.data[, spp. ])
          #resMat[, spp.] <- WN.InitARMA(tsData = 
          #                                data.frame(m.data[, spp. ]),
          #                              order = c(arOrd, 0, maOrd),
          #                              whiteN = TRUE,
          #                              moreOrder = 1,
          #                              updateWN = TRUE)[[2]]
          resMat[, spp. ] <- res.obj[, spp. ]
          
          x.Mat <- matrix(0.0, nrow = nn, ncol = sOrd + !no.mean)
          for (kk in 1:arOrd)
            x.Mat[(kk + 1):nn, kk] <- m.data[1:(nn - kk), spp. ]
          for (kk in 1:maOrd)
            x.Mat[(kk + 1):nn, kk + arOrd] <- resMat[1:(nn - kk), spp. ]
          if (!no.mean) 
            x.Mat[, sOrd + 1] <- 1
          m.Mat <- t(x.Mat) %*% x.Mat /nn
          if (!no.mean)
            m.Mat[sOrd + 1, sOrd + 1] <- m.Mat[sOrd + 1, sOrd + 1] *
                                        (1 - sum(the.ar))^(-2 )
          if (any(eigen(m.Mat)$values <= 1e-5)) 
            for (kk in 1:sOrd) {
              aux.1.a <- m.Mat
              aux.1.a[col(aux.1.a) != row(aux.1.a)] <- (2/3) *
                             aux.1.a[col(aux.1.a) != row(aux.1.a)]
              m.Mat <- aux.1.a
              if (all(eigen(m.Mat)$values > 0))
                break()
            }
          
          se.Mat[spp. , 1:(sOrd + !no.mean)]  <- 
            round(sqrt(m.var[spp. ] * diag(solve(m.Mat))/nn), 3)
        }
      }
      
      object@post <- list(coeff   = round(m.thetas, 3),
                          varcov  = se.Mat,
                          no.mean = no.mean,
                          vglmobj = object@family@vfamily,
                          m.Coef  = m.Coef,
                          flag11  = flag11,
                          extraWN = if (flag11) sd.err else NULL,
                          var  = m.var,
                          NoS  = NoS,
                          Ord  = intOrd,
                          nmes = object@misc$theta.names,
                          proc = process,
                          bigFlag   = bigFlag,
                          m.varff   = object@misc$var.arg,
                          print.AIC = TRUE)
      
      object@post
       
    })






setMethod(showsummaryvglmS4VGAM, signature(VGAMff = "ARMAvgltsmff"),
          function(object, VGAMff = VGAMff,...) {
            
      options(digits = 4)
      NoS     <- object@post$NoS
      Ord     <- object@post$Ord
      m.name  <- paste(object@post$proc, "coeff", sep = "")
      matcoef <- object@post$coeff
      varcoef <- object@post$var
      varCov  <- object@post$varcov
      no.mean <- object@post$no.mean
      m.Coef  <- object@post$m.Coef
      cl.name <- vector("character")
      flag11  <- object@post$flag11
      sd.err  <- object@post$extraWN
      proc    <- object@post$proc
      proc2c  <- object@post$proc
      m.varff <- object@post$m.varff
      p.AIC   <- object@post$print.AIC
      bigFlag <- object@post$bigFlag
      
      if (p.AIC && NoS > 1) {
        loglik.tsff <- if (p.AIC ) object@family@loglikelihood else NULL
        matrix.AIC  <- object@predictors
        res.RespY <- extract.Residuals(object = object, TSprocess = proc)
      }
      
      if ((proc == "AR") || (proc == "MA")) {
        proc <- paste(object@post$proc, "coeff", sep = "")
      } else {
        arOrd <- object@post$Ord[1]
        maOrd <- object@post$Ord[2]
        sOrd  <- arOrd + maOrd
        proc  <- c(paste("ARcoeff", 1:arOrd, sep = ""),
                   paste("MAcoeff", 1:maOrd, sep = ""))
      }
      
      cat("\n-----\n")
      cat("\n** Standard errors based on the asymptotic",
          "\n distribution of the MLE estimates: \n\n")
      
      if (NoS > 1) {
        Ord <- if ((proc2c == "ARMA") || (proc2c == "ARIMA")) 
                          rep(sOrd, NoS) else Ord
        for (kk in 1:NoS) {
          
          ### A little change if an ARMA model is being analyzed. ###
          aux1 <- if ((proc2c == "ARMA") || (proc2c == "ARIMA")) proc else
                                  paste(proc, 1:Ord[kk], sep = "")
          cl.name2 <- c(paste(aux1, kk, sep = "" ), 
                        if (!no.mean) 
                          switch(object@misc$process, 
                                 "AR"    = "drift", 
                                 "MA"    = "Intercept (mean)", 
                                 "ARMA"  = "drift",
                                 "ARIMA" = "drift") else NULL)
          
          cl.name <- c(cl.name, cl.name2)
          preMat  <- matrix(NA, nrow = 2, ncol = Ord[kk] + !no.mean)
          preMat[1, ] <- matcoef[kk, 1:(Ord[kk] + !no.mean) ]
          preMat[2, ] <- varCov[kk, 1:(Ord[kk] + !no.mean)]
          colnames(preMat) <- cl.name2
          
          rownames(preMat) <- c("", "s.e.")
          cat("Response ",kk, ":\n")
          print(preMat, digits = 4)
          names(object@post$var) <- paste("WNvar", 1:NoS, sep = "")
          
          if (flag11) {
            sd.help <- sd.err[kk, ]
            names(sd.help) <- colnames(sd.err)
            if (m.varff) {
              cat("\nEstimated linear predictor for sigma^2:\n")
            } else {
              cat("\nEstimated linear predictor for sigma (errors SD):\n")
            }
            print(sd.help, digits = 4)
            cat("\n")
          } else {
            cat("\n  Estimated variance (sigma^2) of the errors: ", 
                object@post$var[kk], ".\n\n", sep = "")
          }
          
          # AIC-AICC-BIC
          if (p.AIC) {
            options(digits = 6)
            eta.respY <- matrix.AIC[ , 1:(Ord[kk] + 2 - no.mean),
                                     drop = FALSE]
            object@extra$NOS <- 1
            object@extra$y   <- depvar(object)[, kk, drop = FALSE]
            object@extra$M1  <- Ord[kk] + 2 - no.mean
            object@extra$nOrder <- Ord[kk]
            
            # Old : mu = fittedVGAMextra(object)[, kk, drop = FALSE]
            logLik.vgts <- 
              loglik.tsff(mu = NULL, # fitted values computed at @loglikel
                          y = NULL,  # Entered via @extra
                          w = weightsVGAMextra(object,
                                    type.w = "prior")[, kk, drop = FALSE],
                          residuals = FALSE,
                          eta   = eta.respY,
                          extra = object@extra,
                          summation = TRUE)
            
            cat("\nLoglikelihood:", logLik.vgts, "\n")
            
            n.n   <- nrow(object@y[, kk, drop = FALSE])
            AIC1  <- (-2) * logLik.vgts + 2 * (Ord[kk] + 2 - no.mean)
            k.k   <- Ord[kk] + 2 - !no.mean
            
            m.AIC <- c( AIC1,
                        AIC1 + 2 * k.k * (k.k + 1)/(n.n - k.k - 1),
                        (-2) * logLik.vgts + 
                          (Ord[kk] + !no.mean) * log(n.n) ) 
            cat("AIC: " , round(m.AIC[1], 6), ", ",  
                "AICc: ", round(m.AIC[2], 6), ", ",
                "BIC: " , round(m.AIC[3], 6), "\n\n", sep = "")
          }
          matrix.AIC <- matrix.AIC[, -(1:(Ord[kk] + 2 - no.mean)),
                                   drop = FALSE]
          options(digits = 4)
        }
    
      } else {
        
        if (!bigFlag) {
          ### A little change if an ARMA model is being analyzed. ###
          Ord <- if ( (proc2c == "ARMA") ||
                      (proc2c == "ARIMA")) sOrd else Ord
          nam.help <- if (!(proc2c == "AR" || proc2c == "MA")) proc else 
            paste(proc, 1:Ord, sep = "")
          cl.name <- c(nam.help,
                       if (!no.mean)
                         switch(object@misc$process, 
                                "AR"    = "drift", 
                                "MA"    = "Intercept (mean)", 
                                "ARMA"  = "drift",
                                "ARIMA" = "drift") else NULL)
          
          preMat <- matrix(NA, nrow = 2, ncol = Ord + !no.mean)
          preMat[1, ] <- matcoef[1, 1:(Ord + !no.mean)]
          preMat[2, ] <-  varCov[1, 1:(Ord + !no.mean)]
          
          colnames(preMat) <- cl.name # drift
          rownames(preMat) <- c("", "s.e.")
          print(preMat, digits = 4)
          names(object@post$var) <- c("") 
        } else {
          if (!no.mean) {
            cl.name <- switch(object@misc$process, 
                              "AR"    = "drift", 
                              "MA"    = "Intercept (mean)", 
                              "ARMA"  = "drift",
                              "ARIMA" = "drift",
                              "pureGARCH" = "Intercept (mean)")
            preMat <- rbind(matcoef, varCov)
            colnames(preMat) <- cl.name # drift
            rownames(preMat) <- c("", "s.e.")
            print(preMat, digits = 4)
            names(object@post$var) <- c("")
          }
        }
        
        if (flag11) {
          cat("\n  Estimated linear predictor of sigma^2 (SD errors):\n")
          print(sd.err[1, ], digits = 4)
          cat("\n")
        } else {
          cat("\n  Estimated variance (sigma^2) of the errors: ", 
              object@post$var, ".\n\n", sep = "")
        }
        
        # AIC-AICC-BIC
        if (p.AIC) {
          options(digits = 6)
          n.n   <- nrow(object@y)
          AIC1  <- (-2) * logLik.vlm(object, ...) + 
                                         2 * (sum(Ord) + 2 - no.mean)
          
          cat("\nLoglikelihood:", logLik.vlm(object, ...), "\n")
          
          k.k   <- sum(Ord) + 2 - no.mean
          m.AIC <- c( AIC1,
                      AIC1 + 2 * k.k * (k.k + 1)/(n.n - k.k -1),
                      (-2) * logLik.vlm(object, ...) + 
                        (sum(Ord) + 2 - no.mean) * log(n.n) ) 
          cat("AIC: " , round(m.AIC[1], 6), ", ",  
              "AICc: ", round(m.AIC[2], 6), ", ",
              "BIC: " , round(m.AIC[3], 6), "\n", sep = "")
          options(digits = 6)
        }
      }
      
      
      # IF NO pure GARCH process
      if (!bigFlag) {
        ### Checks on stationarity / invertibility 2016/March/02 ###
        m.Coef  <- object@post$m.Coef
        lastCoe <- numeric(0)
        
        if (object@post$proc == "ARMA" || object@post$proc == "ARIMA") {
          coefMa <- coefAr  <- numeric(0)
          
          for (kk in 1:NoS) {
            aux1  <- c(matcoef[kk, ])[1:(arOrd + maOrd)]
            coefAr <- c(coefAr, aux1[1:arOrd])
            coefMa <- c(coefMa, aux1[-(1:arOrd)]) 
          }
          
          cat("-----\n")
          cat("\n** Summary of checks on stationarity /",
              "invertibility:\n")
          checkAR <- checkTS.ffs(thetaEst = coefAr,
                                 tsclass  = "AR",
                                 NofS     = NoS, 
                                 chOrder  = arOrd , 
                                 pRoots   = TRUE, 
                                 retmod   = TRUE)
          checkMA <- checkTS.ffs(thetaEst = coefMa, 
                                 tsclass  = "MA",
                                 NofS     = NoS, 
                                 chOrder  = maOrd ,
                                 pRoots   = TRUE,
                                 retmod   = TRUE)
          
          if ( any(checkAR < 1 + 5e-3))  {   # 2015/11/04  
            if ( arOrd != 1)
              warning("Accuracy of parameter estimates may be low since",
                      " some root(s) of the polynomial ", "\n",
                       " 1 + ARcoeff1 * z + ... + ARcoeff",arOrd ," * z^",
                       arOrd, ", \n", 
                       " lie inside the unit circle.") else
               warning("Accuracy of parameter estimates may be low since",
                       " some root(s) of the polynomial ", "\n",
                       " 1 + ARcoeff",arOrd ," * z", "\n",
                       " lie inside the unit circle.") 
          } 
          
          if ( any(checkMA < 1 + 5e-3))  {  # 2015/11/04
            if ( maOrd != 1)
              warning("Accuracy of parameter estimates may be low since",
                      " some root(s) of the polynomial ", "\n",
                      " 1 + MAcoeff1 * z + ... + MAcoeff",maOrd," * z^",
                      maOrd , "\n",
                      " lie inside the unit circle.") else
               warning("Accuracy of parameter estimates may be low since",
                       " some root(s) of the polynomial ", "\n",
                       " 1 + MAcoeff",maOrd ," * z", "\n",
                       " lie inside the unit circle.") 
          }
          
        } else { ## AR or MA ##
          
          for (ii in 1:NoS) {
            aux5 <- m.Coef[1:Ord[ii]]
            if (length(aux5) < max(Ord))
              aux5 <- c(aux5, rep(0, max(Ord) - length(aux5)))
            
            lastCoe <- c(lastCoe, aux5)
            m.Coef <- m.Coef[-(1:Ord[ii])]
          }
          
          cat("-----\n")
          cat("\n** Summary of checks on stationarity /",
              "invertibility:\n")
          names(lastCoe) <- rep("", times = length(lastCoe))
          options(digits = 4)   
          flag    <- FALSE
          myroots <- 
            checkTS.ffs(thetaEst = lastCoe,
                        tsclass  = object@post$proc,
                        NofS     = NoS,
                        chOrder  = max(Ord),
                        pRoots   = TRUE,
                        retmod   = TRUE)
          for ( jj in 1:NoS ) 
            if (all( myroots[, jj ] > 0 ) && 
                any( myroots[, jj ] <= 1 + 5e-3)) {
              flag <- TRUE
              break()
            } 
          if ( flag && Ord[jj] == 1 )
            warning("\nAccuracy of parameter estimates may be low ",
                "since the root of the polynomial ", "\n",
                  " 1 + ARcoeff1 * z,", "\n",
                   "Response ", jj, ", ",
                    "lies inside the unit circle.") else
               if ( flag )
                  warning("\nAccuracy of parameter estimates may be low ",
                          "since some root(s) of the polynomial ",  "\n",
                          " 1 + ARcoeff1 * z + ... + ARcoeff", Ord[jj] ,
                          "* z^", Ord[jj] , ", \n",
                          "Response ", jj, " lie inside the unit circle.")
        }
        options(digits = 7)
      } else {
        # If PURE GARCH ...
        cat("\n** Summary of checks on stationarity /",
            "invertibility:\n")
        cat("\n No ARMA component involved.")
        
      }
      
  })

