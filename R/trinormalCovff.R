##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0


trinormalCovff.control <- function(save.weights = TRUE,
                              summary.HDEtest = FALSE,...) {
  list(save.weights = save.weights,
       summary.HDEtest = summary.HDEtest,...)
}

trinormalCovff <- function(zero = c("var", "cov"),
                  lmean = "identitylink",
                      lvar  = "loglink",
                      lcov  = "identitylink") {
  
  trinorm <-  TRUE ### If FALSE, then 'trinormalCovff()' mimics VGAM::binormal.
                   ### It may be and argument later (included at early stages for
                   ### testing purposes only).
  
  if (!is.logical(trinorm))
    stop("Wrong input for 'trinorm'.")
  
  loffd <- lcov; rm(lcov)
  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")
  
  lvar <- as.list(substitute(lvar))
  evar <- link2list(lvar)
  lvar <- attr(evar, "function.name")
  
  loffd <- as.list(substitute(loffd))
  eoffd <- link2list(loffd)
  loffd <- attr(eoffd, "function.name")
  
  alo <-
    new("vglmff",
      blurb = c("Trivariate normal distribution \n\n",
                "Links:   ",
                namesof("mean", lmean, earg = emean), ", ",
                namesof("variance", lvar, earg = evar), ", ",
                namesof("Error covariances", loffd, earg = eoffd), ".\n",
                "Mean:    c(mean1, mean2, mean3)^T"),
      
      
      
      constraints = eval(substitute(expression({
        M1 <- if ( .trinorm ) 9 else 5
        constraints <- 
          cm.zero.VGAM(constraints, x = x, zero = .zero , M = M ,
                       predictors.names = parameters.names, 
                       M1 = M1)
      }), list( .zero = zero, .trinorm = trinorm ))),
      
      
      
      infos = eval(substitute(function(...){
        
        my.names <- c(paste("mean", 1:(2 + .trinorm ), sep = ""),
                      paste("var", 1:(2 + .trinorm ), sep = ""),
                      paste("cov", 1:(2 + .trinorm ), sep = ""))
        
        list(M1 = if ( .trinorm ) 9 else 5,
             Q1 = 2 + .trinorm ,
             expected = TRUE,
             multipleResponses = FALSE,
             lmean = .lmean ,
             lvar  = .lvar ,
             loffd = .loffd ,
             parameters.names = my.names,
             zero  = .zero )
      }, list( .zero = zero , .lmean = lmean , .lvar = lvar ,
               .loffd = loffd , .trinorm = trinorm ))),
      
      
      
      initialize = eval(substitute(expression({
        
        Q1 <- 2 + .trinorm
        M <- M1 <- if (.trinorm ) 9 else 5
        check <- w.y.check(w = w, y =y, 
                           Is.positive.y = FALSE, 
                           ncol.w.max = 1, 
                           ncol.y.max = Q1, 
                           ncol.y.min = Q1,
                           out.wy = TRUE, 
                           colsyperw = Q1, 
                           maximize = TRUE)
        
        w <- check$w
        y <- check$y
        NOS.check   <- ncol(y)
        nn <- nrow(y)
        
        if ((NOS.check == 2) && ( .trinorm ))
          stop("Wrong input in the number of responses.")
        
        if ((NOS.check != 2) && (NOS.check != 3))
          stop("Wrong input in the response. " ,
               "Only bivariate and trivariate Normal handled.")
        
        NOS.b <- 2 + .trinorm
        parameters.names <- c(paste("mean", 1:NOS.b, sep = ""),
                              paste("var", 1:NOS.b, sep = ""))
        pre.name <- character(0)
        
        for (ii in 1:(NOS.b - 1)) 
          pre.name <- c(pre.name,
                        paste(paste("cov", ii, sep = ""), 
                              (ii + 1):NOS.b, sep = ""))
        parameters.names <- c(parameters.names, pre.name)
        
        
        predictors.names <- 
          c(namesof(paste("mean", 1:NOS.b, sep = ""),
                    link = .lmean , earg = .emean , tag = FALSE),
            namesof(paste("var", 1:NOS.b, sep = ""),
                    link = .lvar , earg = .evar , tag = FALSE))
        
        pre.name <- character(0)
        for (ii in 1:(NOS.b - 1)) 
          pre.name <- c(pre.name,
                        paste(paste("cov", ii, sep = ""), 
                              (ii + 1):NOS.b, sep = ""))
        rm(NOS.b)
        
        predictors.names <- c(predictors.names,
                              namesof(pre.name, .loffd , .eoffd ))
        
        
        extra$colnames.y <- colnames(y)
        
        ini.mean <- matrix(apply(y, 2, function(x) weighted.mean(x, w = w)),
                           nn, 2 + .trinorm , byrow = TRUE)
        ini.var  <- cov(y)
        
        real.var <- if (NOS.check == 2) matrix(ini.var[1, 2], nn, 1) else
                matrix(c(ini.var[1, 2], ini.var[1, 3], ini.var[2, 3]),
                       nn, 2 + .trinorm , byrow = TRUE)
        ini.var <- matrix(diag(ini.var), nn, 2 + .trinorm , byrow = TRUE)
        
        etastart <- cbind(theta2eta(ini.mean, .lmean , earg = .emean ),
                      theta2eta(ini.var, .lvar , earg = .evar ),
                      theta2eta(real.var, .loffd , earg = .eoffd ))
        etastart
        
      }), list( .trinorm = trinorm , 
                .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                .emean = emean , .evar = evar , .eoffd = eoffd ))),
      
      
      
      
      linkinv = eval(substitute(function(eta, extra = NULL) {
        
        NOS <- 1
        eta2theta(eta[, 1:(2 + .trinorm ), drop = FALSE],
                              .lmean , earg = .emean )
        
      }, list( .trinorm = trinorm ,
                .lmean =lmean , .emean = emean ))),
      
      
      
      
      last = eval(substitute(expression({
        
        if ( .trinorm ) {
          misc$link <- c("mean1" = .lmean , "mean2" = .lmean ,
                         "mean3" = .lmean ,
                         "var1"  = .lvar , "var2"  = .lvar , 
                         "var3"  = .lvar ,
                         "cov12" = .loffd , "cov13" = .loffd ,
                         "cov23" = .loffd )
          
          misc$earg <- c("mean1" = .emean , "mean2" = .emean ,
                         "mean3" = .emean ,
                         "var1"  = .evar , "var2"  = .evar , 
                         "var3"  = .evar ,
                         "cov12" = .eoffd , "cov13" = .eoffd ,
                         "cov23" = .eoffd )
        } else {
          misc$link <- c("mean1" = .lmean , "mean2" = .lmean ,
                         "var1"  = .lvar , "var2"  = .lvar , 
                         "cov12" = .loffd )
          
          misc$earg <- c("mean1" = .emean , "mean2" = .emean ,
                         "var1"  = .evar , "var2"  = .evar , 
                         "cov12" = .eoffd )
        }
        
        misc$expected <- TRUE
        misc$multipleResponses <- FALSE
        
      }), list( .trinorm = trinorm , 
                .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                .emean = emean , .evar = evar , .eoffd = eoffd ))),
      
      
      
      loglikelihood = eval(substitute(
        function(mu, y , w,
                 residuals = FALSE, eta,
                 extra = NULL,
                 summation = TRUE) {
          
          nn <- nrow(y)
          
          if (.trinorm ) {
            vec.mean <- eta2theta(eta[, 1:3], .lmean , earg = .emean )
            var.only  <- eta2theta(eta[, 4:6], .lvar  , earg = .evar )
            offds <- eta2theta(eta[, 7:9], .loffd , earg = .eoffd )
          } else {
            vec.mean <- cbind(eta2theta(eta[, 1], .lmean , earg = .emean ),
                              eta2theta(eta[, 2], .lmean , earg = .emean ))
            var.only <- cbind(eta2theta(eta[, 3], .lvar  , earg = .evar ),
                              eta2theta(eta[, 4], .lvar  , earg = .evar ))
            offds <- eta2theta(eta[, 5], .loffd , earg = .eoffd )
          }
          
          if (residuals) {
            stop("loglikelihood residuals not implemented yet")
          } else {
            
            ll.elts <- c(w) * dmultinorm(vec.x = y, vec.mean = vec.mean,
                                mat.cov = cbind(var.only, offds),
                                log = TRUE)
            if (summation) {
              sum(ll.elts)
            } else {
              ll.elts
            }
          }
        }, list( .trinorm = trinorm , 
                 .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                 .emean = emean , .evar = evar , .eoffd = eoffd  ))),
      
      
      
      vfamily = c("trinormalCovff"),
      
      
      
      deriv = eval(substitute(expression({
        
        nn <- nrow(eta)
        counter <- cbind(1:nn)
        
        if (.trinorm ) {
          vec.mean <- eta2theta(eta[, 1:3], .lmean , earg = .emean )
          var.only  <- eta2theta(eta[, 4:6], .lvar  , earg = .evar )
          offds <- eta2theta(eta[, 7:9], .loffd , earg = .eoffd )
          
          
          bigMat.der <- matrix(NA_real_, nn, 9)
          mat11 <- cbind(var.only, offds)
          temp1 <- apply(mat11, 1 , function(x) {
            matrix(c(x[1], x[4], x[5],
                     x[4], x[2], x[6],
                     x[5], x[6], x[3]), 3, 3, byrow = TRUE)
          })
          
          
          attr(temp1, "dim") <- c(3, 3, nn)
          # R Inverse
          temp1 <- apply(temp1, 3, function(x) solve(x))
          attr(temp1, "dim") <- c(3, 3, nn)
          
          dV.ds21 <- matrix(0, 3, 3); dV.ds21[1, 1] <- 1
          dV.ds22 <- matrix(0, 3, 3); dV.ds22[2, 2] <- 1
          dV.ds23 <- matrix(0, 3, 3); dV.ds23[3, 3] <- 1
          
          dV.dthe1 <- matrix(0, 3, 3); dV.dthe1[1, 2] <- dV.dthe1[2, 1] <- 1
          dV.dthe2 <- matrix(0, 3, 3); dV.dthe2[1, 3] <- dV.dthe2[3, 1] <- 1
          dV.dthe3 <- matrix(0, 3, 3); dV.dthe3[2, 3] <- dV.dthe3[3, 2] <- 1
          
          ###### For sigma^21
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*% dV.ds21 )))
          })
          #dim(temp2) <- c(3, 3, 10)
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.ds21 %*% x 
          })
          
          dim(temp3) <- c(3, 3, nn)
          #temp3
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          bigMat.der[, 4] <- temp2 + temp3
          
          
          ###### For sigma^22
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*% dV.ds22 )))
          })
          
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.ds22 %*% x 
          })
          dim(temp3) <- c(3, 3, nn)
          
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          
          bigMat.der[, 5] <- temp2 + temp3
          
          
          ###### For sigma^23
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*% dV.ds23 )))
          })
          
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.ds23 %*% x 
          })
          dim(temp3) <- c(3, 3, nn)
          
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          
          bigMat.der[, 6] <- temp2 + temp3
          
          
          ###### For theta1
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*%  dV.dthe1 )))
          })
          
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.dthe1 %*% x 
          })
          dim(temp3) <- c(3, 3, nn)
          
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          
          bigMat.der[, 7] <- temp2 + temp3
          
          
          
          ###### For theta2
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*%  dV.dthe2 )))
          })
          
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.dthe2 %*% x 
          })
          dim(temp3) <- c(3, 3, nn)
          
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          
          bigMat.der[, 8] <- temp2 + temp3
          
          
          ###### For theta3
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*%  dV.dthe3 )))
          })
          
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.dthe3 %*% x 
          })
          dim(temp3) <- c(3, 3, nn)
          
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          
          bigMat.der[, 9] <- temp2 + temp3
          
          
          ### For the means
          temp3 <- apply(counter, 1, function(x) {
            temp1[, , x] %*% cbind((y - vec.mean))[x, ]
          })
          
          bigMat.der[, 1:3] <- t(temp3)
          
          
        } else {
          
          vec.mean <- cbind(eta2theta(eta[, 1], .lmean , earg = .emean ),
                            eta2theta(eta[, 2], .lmean , earg = .emean ))
          var.only <- cbind(eta2theta(eta[, 3], .lvar  , earg = .evar ),
                            eta2theta(eta[, 4], .lvar  , earg = .evar ))
          offds <- eta2theta(eta[, 5], .loffd , earg = .eoffd )
          
          bigMat.der <- matrix(NA_real_, nn, 5)
          mat11 <- cbind(var.only, offds)
          temp1 <- apply(mat11, 1 , function(x) {
            matrix(c(x[1], x[3],
                     x[3], x[2]), 2, 2, byrow = TRUE)
          })
          
          attr(temp1, "dim") <- c(2, 2, nn)
          
          # R Inverse
          temp1 <- apply(temp1, 3, function(x) solve(x))
          attr(temp1, "dim") <- c(2, 2, nn)
          
          dV.ds21 <- matrix(0, 2, 2); dV.ds21[1, 1] <- 1
          dV.ds22 <- matrix(0, 2, 2); dV.ds22[2, 2] <- 1
          
          dV.dthe1 <- matrix(0, 2, 2)
          dV.dthe1[1, 2] <- dV.dthe1[2, 1] <- 1
          
          ###### For sigma^21
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*% dV.ds21 )))
          })
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.ds21 %*% x 
          })
          
          dim(temp3) <- c(2, 2, nn)
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          bigMat.der[, 3] <- temp2 + temp3
          
          
          
          ###### For sigma^22
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*% dV.ds22 )))
          })
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.ds22 %*% x 
          })
          
          dim(temp3) <- c(2, 2, nn)
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          bigMat.der[, 4] <- temp2 + temp3
          
          
          ###### For theta1
          temp2 <- apply(temp1, 3, function(x) {
            (-0.5) *  sum(diag(( x %*% dV.dthe1 )))
          })
          temp3 <- apply(temp1, 3, function(x) {
            x %*% dV.dthe1 %*% x 
          })
          
          dim(temp3) <- c(2, 2, nn)
          temp3 <- apply(counter, 1, function(x) {
            (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
          })
          
          bigMat.der[, 5] <- temp2 + temp3
          
          
          ### For the means
          temp3 <- apply(counter, 1, function(x) {
            temp1[, , x] %*% cbind((y - vec.mean))[x, ]
          })
          
          bigMat.der[, 1:2] <- t(temp3)
          
          
        }
        
        
        ### Der theta / d.eta
        dmeans.deta <- dtheta.deta(vec.mean, .lmean , .emean )
        dvars.deta  <- dtheta.deta(var.only, .lvar, .evar )
        doffds.deta <- dtheta.deta(offds, .loffd , .eoffd )
       
        c(w) * bigMat.der * cbind(dmeans.deta,
                                  dvars.deta,
                                  doffds.deta)
        
      }), list( .trinorm = trinorm , 
                .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                .emean = emean , .evar = evar , .eoffd = eoffd  ) )),
      
      
      weight = eval(substitute(expression({
        
        
        if ( .trinorm ) {
          
          wz <- matrix(0.0, nn, dimm(9))
          wz[, 1:3] <- cbind(temp1[1, 1, ],
                             temp1[2, 2, ],
                             temp1[3, 3, ]) * dmeans.deta^2
          
          ned2l.ds21 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.ds21 %*% 
                               temp1[, , x] %*% dV.ds21))
          })
          
          ned2l.ds22 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.ds22 %*% 
                               temp1[, , x] %*% dV.ds22))
          })
          
          ned2l.ds23 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.ds23 %*% 
                               temp1[, , x] %*% dV.ds23))
          })
          
          wz[, 4:6] <- cbind(ned2l.ds21,
                             ned2l.ds22,
                             ned2l.ds23) * dvars.deta^2
          
          
          ned2l.dthe1 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.dthe1 %*% 
                               temp1[, , x] %*% dV.dthe1))
          })
          
          ned2l.dthe2 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.dthe2 %*% 
                               temp1[, , x] %*% dV.dthe2))
          })
          
          ned2l.dthe3 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.dthe3 %*% 
                               temp1[, , x] %*% dV.dthe3))
          })
          
          wz[, 7:9] <- cbind(ned2l.dthe1,
                             ned2l.dthe2,
                             ned2l.dthe3) * doffds.deta^2
          
          
        } else {
          
          wz <- matrix(0.0, nn, dimm(5))
          wz[, 1:2] <- cbind(temp1[1, 1, ],
                             temp1[2, 2, ]) * dmeans.deta^2
          
          ned2l.ds21 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.ds21 %*% 
                               temp1[, , x] %*% dV.ds21))
          })
          
          ned2l.ds22 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.ds22 %*% 
                               temp1[, , x] %*% dV.ds22))
          })
          
          wz[, 3:4] <- cbind(ned2l.ds21,
                             ned2l.ds22) * dvars.deta^2
          
          ned2l.dthe1 <- apply(counter, 1, function(x) {
            (0.5) * sum(diag(temp1[, , x] %*% dV.dthe1 %*% 
                               temp1[, , x] %*% dV.dthe1))
          })
          
          wz[, 5] <- cbind(ned2l.dthe1) * doffds.deta^2
        }
        
        c(w) * wz
        
      }), list( .trinorm = trinorm , 
                .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                .emean = emean , .evar = evar , .eoffd = eoffd  ))) )

  alo   
}


