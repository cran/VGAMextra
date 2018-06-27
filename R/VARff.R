##########################################################################
# These functions are 
# Copyrigth (C) 2014-2018 V. Miranda and T. W. Yee. University of Auckland
# All rights reserved.

VARff <- function(VAR.order = 1,
                  zero = c("var", "cov"),
                  lmean = "identitylink",
                  lvar  = "loge",
                  lcov  = "identitylink") {
  
  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")
  
  lvar <- as.list(substitute(lvar))
  evar <- link2list(lvar)
  lvar <- attr(evar, "function.name")
  
  loffd <- lcov; rm(lcov)
  loffd <- as.list(substitute(loffd))
  eoffd <- link2list(loffd)
  loffd <- attr(eoffd, "function.name")
  
  ecm.order <- VAR.order; rm(VAR.order)
  V.class  <- c("ECM", "VAR")[2]
  trinorm  <- FALSE
  noRemove <- FALSE
  ordtsDyn <- 0
  lag.res  <- 1
  if (!Is.Numeric(ordtsDyn, length.arg =  1) || ordtsDyn < 0)
    stop("Bad input for 'ordtsDyn'")
  
  resids.pattern <- match.arg("intercept",  c("neither", "intercept",
                                              "trend", "both"))[1]
  r.pat <- resids.pattern; rm(resids.pattern)
  V.class <- match.arg(V.class, c("ECM", "VAR"))[1]
  
  if (!Is.Numeric(lag.res, length.arg = 1, 
                  isInteger = TRUE) || lag.res < 1)
    stop("Bad enter for argument 'lag.res'.")
  
  if (!is.logical(noRemove))
    stop("Wrong input for argument 'noRemove'.")
  
  if ((length(my.ord <- ecm.order ) != 1 ) || any(ecm.order <= 0))
    stop("Wrong input for argument 'VAR.order'.")
  
  temp1 <- switch(V.class,
                  ECM = "Error Correction Model (Engle -- Granger)\n",
                  VAR = c("Vector Autoregressive Model of order ", 
                          ecm.order, "\n"))
  alo <-
    new("vglmff",
        blurb = c(temp1,
                  "Links:   ",
                  namesof("Cond.Mean", lmean, earg = emean), ", ",
                  namesof("variance", lvar, earg = evar), ", ",
                  namesof("ErrorCovariances", loffd, earg = eoffd)),
        
        
        constraints = eval(substitute(expression({
          
          M <- M1 <- extra$n.means + ( extra$n.means *
                                         (extra$n.means + 1) )/2
          constraints <- 
            cm.zero.VGAM(constraints, x = x, zero = .zero , M = M ,
                         predictors.names = parameters.names, 
                         M1 = M1)
        }), list( .zero = zero, .trinorm = trinorm ))),
        
        
        
        
        first = eval(substitute(expression({
          if (!length(ncol(y)))
            stop("\n Refer to uninormal() from VGAM to fit the",
                 " univariate Normal.")
          
          class.in <- .V.class; res.p <- .r.pat
          NOS.f <- NCOL(y); nn <- NROW(y)
          int.ord <- rep ( .my.ord , NOS.f )
          T.trend <- 1:nn
          
          ### FALSE at all times for VARff
          if (class.in %in% c("ECM")) {
            my.form <- switch(res.p,
                              "neither" = y[, 2] ~ y[, 1] - 1,
                              "intercept" = y[, 2] ~ y[, 1],
                              "trend" = y[, 2] ~ y[, 1] + T.trend - 1,
                              "both"  = y[, 2] ~ y[, 1] + T.trend )
            int.regre  <- lm(my.form)
            int.resids <- cbind(residuals(int.regre))
            extra$coint.vec <- coef(int.regre)
            names(extra$coint.vec) <- switch(res.p,
                         "neither" = "betaY2",
                         "intercept" = c("(Intercept)", "betaY2"),
                         "trend" = c("betaY2", "T.trend"),
                          "both"  =  c("(Intercept)", "betaY2", "T.trend"))
            rm(T.trend, int.regre)
            
            rem.res <- .lag.res - 1
            temp1   <- if ( !int.ord[1] ) NULL else 
              embed(diff(y[, 1]), 1 + max(int.ord))
            temp2   <- if ( !int.ord[2] ) NULL else
              embed(diff(y[, 2]), 1 + max(int.ord))
            x.temp  <- embed(x, 2)[-c(1:max(int.ord)), 1:NCOL(x),
                                   drop = FALSE]
            
            if ( .ordtsDyn ) {
              temp3 <- int.resids[c((.ordtsDyn + 1):
                                      (nn - .lag.res - 1 + .ordtsDyn))]
            } else {
              # As Pfaff(2011)
              temp3 <-   int.resids[c(1:(nn - .lag.res - 1))]
            }
            
            if (max(int.ord) - .lag.res > 0) {
              temp3 <- temp3[-c(1:(max(int.ord) - .lag.res))]
            } else {
              if (max(int.ord) - .lag.res < 0) {
                temp1 <- temp1[-c(1:abs(max(int.ord) - .lag.res)), ,
                               drop = FALSE]
                temp2 <- temp2[-c(1:abs(max(int.ord) - .lag.res)), , 
                               drop = FALSE]
                x.temp <- x.temp[-c(1:abs(max(int.ord) - .lag.res)), , 
                                 drop = FALSE]
              } else {
                
              }
            }
            
            x.matrix <- cbind(x.temp, temp3,
                              if (!int.ord[1]) NULL else 
                                temp1[, 2:(1 + int.ord[1]), drop = FALSE],
                              if (!int.ord[2]) NULL else 
                                temp2[,  2:(1 + int.ord[2]), drop = FALSE])
            
            colnames(x.matrix) <- 
              c(colnames(x), paste("ErrorsLag", .lag.res , sep = ""),
                #paste("dify", 1:NOS.f, sep = ""), 
                if ( !int.ord[1] ) NULL else
                  paste("diffy1Lag", 1:int.ord[1], sep = ""),
                if ( !int.ord[2] ) NULL else
                  paste("diffy2Lag", 1:int.ord[2], sep = ""))
            
            list.names <- vector("list", NCOL(x) + 1 + sum(int.ord))
            names(list.names) <- colnames(x.matrix)
            for (ii in 1:( NCOL(x) + 1 + sum(int.ord)))
              list.names[[ii]] <- ii
            attr(x.matrix, "assign") <- list.names
            y <- cbind(temp1[, 1], temp2[, 1])
            colnames(y) <- paste("diffy", 1:NOS.f, sep = "")
            
            x <- x.matrix
            w <- w[c(1:NROW(y))]
            
          } else {
            
                  ### TRUE at all times ### 
            #x.temp <- if ( !int.ord[1] ) NULL else   20180227
            #  embed(y[, 1], 1 + max(int.ord))[, 2:(1 + int.ord[1])]
            
            x.temp <- WN.lags(cbind(y[, 1]), int.ord[1], 
                              to.complete = y[c(1:int.ord[1]), 1])
            
            for (ii in 2:NOS.f)
              x.temp <- cbind(x.temp, if (!int.ord[ii]) NULL else 
                      WN.lags(cbind(y[, ii]), int.ord[ii],
                              to.complete = y[c(1:int.ord[ii]), ii]))
            temp.names <- c(colnames(x))
            for (ii in 1:NOS.f)
              temp.names <- c(temp.names,
                              paste(paste(paste("Y", ii, sep = ""),
                                          "Lag", sep = ""),
                                    1:int.ord[ii], sep = ""))
            
            x.matrix <- cbind(x, x.temp)
            colnames(x.matrix) <- temp.names 
            
            list.names <- vector("list", NCOL(x) + sum(int.ord))
            names(list.names) <- colnames(x.matrix)
            for (ii in 1:( NCOL(x) + sum(int.ord)))
              list.names[[ii]] <- ii
            attr(x.matrix, "assign") <- list.names
            x <- x.matrix
          }
          
        }), list( .my.ord = my.ord , .lag.res = lag.res ,
                  .V.class = V.class , .r.pat = r.pat ,
                  .ordtsDyn = ordtsDyn ))),
        
        
        
        
        infos = eval(substitute(function(...){
          list(M1 = NULL, # Can't set the expression for M1, Q1. No variables
               Q1 = NULL, # are allowed to be read in this slot,
               expected = TRUE,
               multipleResponses = FALSE,
               lmean = .lmean ,
               lvar  = .lvar ,
               loffd = .loffd ,
               #parameters.names = parameters.names,
               zero  = .zero )
          
        }, list( .zero = zero , .lmean = lmean , .lvar = lvar ,
                 .loffd = loffd , .trinorm = trinorm ))),
        
        
        
        
        
        initialize = eval(substitute(expression({
          
          Q1 <- ncol(y)
          M <- M1 <- Q1 + Q1 * (Q1 + 1) / 2
          
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
          NOS.b <- extra$n.means  <- extra$NOS.b <- Q1
          nn <- nrow(y)
          extra$M <- M
          
          der.list <- vector("list", extra$M)
          der.MUs <- vector("list", extra$M)
          der.Vars <- vector("list", extra$M)
          
          dMUs <- array(0, dim = c(extra$NOS.b, 1, extra$NOS.b))
          for (ii in 1:extra$NOS.b) 
            der.list[[ii]] <- dMUs[, , ii] <- 
            diag(extra$NOS.b)[, ii, drop = FALSE]
          
          dVs <- array(0, dim = c(extra$NOS.b, extra$NOS.b, extra$NOS.b))
          
          for(ii in 1:extra$NOS.b) {
            dVs[, , ii][ii, ii] <- 1
            der.list[[ii + extra$NOS.b]] <- dVs[, ,  ii]
          }
          
          dTs <- array(0, dim = c(extra$NOS.b, extra$NOS.b,
                                  M1 - 2 * extra$NOS.b))
          
          hlcount <- 1
          for (ii in 1:(extra$NOS.b - 1)) {
            for (jj in (ii + 1):extra$NOS.b) {
              dTs[, , hlcount][ii, jj] <- 1
              dTs[, , hlcount][jj, ii] <- 1
              hlcount <- hlcount + 1
            }
          }; rm(hlcount)
          
          for (ii in 1:dim(dTs)[3])
            der.list[[ii + 2 * extra$NOS.b]] <- dTs[, , ii]
          
          for (ii in 1:extra$NOS.b)
            der.MUs[[ii]] <- cbind(dMUs[, , ii])
          for (ii in (1 + extra$NOS.b):M)
            der.MUs[[ii]] <- matrix(0, extra$NOS.b, 1)
          
          for (ii in 1:extra$NOS.b)
            der.Vars[[ii]] <- matrix(0, extra$NOS.b, extra$NOS.b)
          
          for (ii in (1 + extra$NOS.b):(2 * extra$NOS.b))
            der.Vars[[ii]] <- der.list[[ii]]
          
          for (ii in (2 * extra$NOS.b + 1):M)
            der.Vars[[ii]] <- der.list[[ii]]
          
          extra$derVs  <- dVs
          extra$derTs  <- dTs
          extra$derMUs <- dMUs
          rm(dMUs)
          extra$der.MUs <- der.MUs
          extra$der.Vars <- der.Vars
          
          parameters.names <-  c( if ( .V.class %in% c("ECM") )
                         paste("Diff", 1:NOS.b, sep = "") else 
                                paste("Mean", 1:NOS.b, sep = ""),
                                      paste("var", 1:NOS.b, sep = ""))
          
          pre.name <- character(0)
          
          for (ii in 1:(NOS.b - 1)) 
            pre.name <- c(pre.name, paste(paste("cov", ii, sep = ""), 
                                  (ii + 1):NOS.b, sep = ""))
          
          parameters.names <- c(parameters.names, pre.name)
          
          predictors.names <- 
            c(namesof(c( if ( .V.class %in% c("ECM") )
              paste("Diff", 1:NOS.b, sep = "") else 
                paste("Mean", 1:NOS.b, sep = "")),
              link = .lmean , earg = .emean , tag = FALSE),
              namesof(paste("var", 1:NOS.b, sep = ""),
                      link = .lvar , earg = .evar , tag = FALSE))
          
          pre.name <- character(0)
          for (ii in 1:(NOS.b - 1)) 
            pre.name <- c(pre.name,
                          paste(paste("cov", ii, sep = ""), 
                                (ii + 1):NOS.b, sep = ""))
          
          predictors.names <- c(predictors.names,
                                namesof(pre.name, .loffd , .eoffd ))
          extra$colnames.y <- colnames(y)
          
          ini.mean <- matrix(apply(y, 2, function(x) weighted.mean(x, w = w)),
                             nn, NOS.b , byrow = TRUE)
          ini.var <- cov(y)
          
          real.var <- ini.var[row(ini.var) > col(ini.var)]
          real.var <- matrix(real.var, nn, M - 2 * Q1, byrow = TRUE)
          ini.var <- matrix(diag(ini.var), nn, NOS.b, byrow = TRUE)
          
          etastart <- cbind(theta2eta(ini.mean, .lmean , earg = .emean ),
                            theta2eta(ini.var, .lvar , earg = .evar ),
                            theta2eta(real.var, .loffd , earg = .eoffd ))
          etastart
          
        }), list( .trinorm = trinorm , 
                  .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                  .emean = emean , .evar = evar , .eoffd = eoffd ,
                  .V.class = V.class ))),
        
        
        
        
        linkinv = eval(substitute(function(eta, extra = NULL) {
          eta2theta(eta[, 1:extra$NOS.b, drop = FALSE],
                    .lmean , earg = .emean )
          
        }, list( .trinorm = trinorm ,
                 .lmean =lmean , .emean = emean ))),
        
        
        
        
        last = eval(substitute(expression({
          
          misc$link <- c(rep( .lmean , extra$NOS.b),
                         rep( .lvar , extra$NOS.b),
                         rep( .loffd , M - 2*extra$NOS.b))
          names(misc$link) <- parameters.names
          
          misc$earg <- vector("list", length = M)
          for (ii in 1:extra$NOS.b) {
            misc$earg[[ii]] <- .emean
            misc$earg[[extra$NOS.b + ii]] <- .evar
          }
          
          for (ii in (2 * extra$NOS.b):(M - extra$NOS.b))
            misc$earg[[ii]] <- .eoffd
          names(misc$earg) <- parameters.names
          
          misc$expected <- TRUE
          misc$multipleResponses <- FALSE
          misc$order <-.my.ord
          
          if ( .V.class %in% c("ECM")) {
            cat("\nCo-integrated vector:\n")
            print(extra$coint.vec)
          }
          
          
        }), list( .trinorm = trinorm , .my.ord = my.ord ,
                  .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                  .emean = emean , .evar = evar , .eoffd = eoffd ,
                  .V.class = V.class , .r.pat = r.pat ))),
        
        
        
        loglikelihood = eval(substitute(
          function(mu, y , w,
                   residuals = FALSE, eta,
                   extra = NULL,
                   summation = TRUE) {
            
            nn <- nrow(y)
            
            vec.mean <- eta2theta(eta[, 1:extra$NOS.b],
                                  .lmean , earg = .emean )
            var.only <- eta2theta(eta[, (1 + extra$NOS.b):(2 * extra$NOS.b)],
                                  .lvar  , earg = .evar )
            offds <- eta2theta(eta[, -c(1:(2 * extra$NOS))],
                               .loffd , earg = .eoffd )
            
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
        
        
        
        vfamily = c("VARff"),
        
        
        
        deriv = eval(substitute(expression({
          nn <- nrow(eta)
          counter <- cbind(1:nn)
          M1 <- extra$n.means + ( extra$n.means *
                                    (extra$n.means + 1) )/2
          
          dVs  <- extra$derVs 
          dTs <- extra$derTs
          dMUs <- extra$derMUs
          der.MUs  <- extra$der.MUs
          der.Vars <- extra$der.Vars
          
          vec.mean <- eta2theta(eta[, 1:extra$NOS.b, drop = FALSE],
                                .lmean , earg = .emean )
          var.only <- eta2theta(eta[, (1 + extra$NOS.b):(2 * extra$NOS.b),
                                    drop = FALSE], .lvar  , earg = .evar )
          offds <- eta2theta(eta[, -c(1:(2 * extra$NOS)), drop = FALSE],
                             .loffd , earg = .eoffd )
          
          bigMat.der <- matrix(NA_real_, nn, M)
          nos.int <- extra$NOS.b
          
          mat11 <- cbind(var.only, offds)
          temp1 <-  apply(mat11, 1 , function(x) {
            int.mat <- matrix(0, nos.int, nos.int)
            diag(int.mat) <- x[1:nos.int]
            int.mat[row(int.mat) > col(int.mat)] <- x[-c(1:nos.int)]
            int.mat <- t(int.mat)
            int.mat[row(int.mat) > col(int.mat)] <- x[-c(1:nos.int)]
            int.mat
          }); rm(nos.int)
          
          attr(temp1, "dim") <- c(extra$NOS.b, extra$NOS.b, nn)
          
          temp1 <- apply(temp1, 3, function(x) solve(x))
          attr(temp1, "dim") <- c(extra$NOS.b, extra$NOS.b, nn)
          
          # From @initialize
          dVs <- extra$derVs
          dTs <- extra$derTs
          
          for (ii in 1:(extra$NOS.b)) {
            temp2 <- apply(temp1, 3, function(x) {
              (-0.5) *  sum(diag(( x %*% dVs[, , ii] )))
            })
            
            temp3 <- apply(temp1, 3, function(x) {
              x %*% dVs[, , ii] %*% x 
            })
            
            dim(temp3) <- c(extra$NOS.b, extra$NOS.b, nn)
            
            temp3 <- apply(counter, 1, function(x) {
              (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
            })
            
            bigMat.der[, extra$NOS.b + ii] <- temp2 + temp3
            
          }
          
          for (ii in 1:(M1 - 2 * extra$NOS.b)) {
            
            temp2 <- apply(temp1, 3, function(x) {
              (-0.5) *  sum(diag(( x %*% dTs[, , ii] )))
            })
            
            temp3 <- apply(temp1, 3, function(x) {
              x %*% dTs[, , ii] %*% x 
            })
            
            dim(temp3) <- c(extra$NOS.b, extra$NOS.b, nn)
            
            temp3 <- apply(counter, 1, function(x) {
              (0.5) * (y - mu)[x, ] %*% temp3[, , x] %*% cbind((y - mu)[x, ])
            })
            
            bigMat.der[, 2 * extra$NOS.b + ii] <- temp2 + temp3
            
            
          }
          
          temp3 <- apply(counter, 1, function(x) {
            temp1[, , x] %*% cbind((y - mu))[x, ]
          })
          
          bigMat.der[, 1:(extra$NOS.b)] <- t(temp3)
          
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
          
          wz <- matrix(0.0, nn, dimm(extra$M))
          
          
          for (ii in 1:(extra$NOS.b)) 
            wz[, ii] <- temp1[ii, ii, ] * dmeans.deta[, ii]^2
          
          for (ii in 1:extra$NOS.b)
            wz[, (ii + extra$NOS.b)] <- (0.5) * temp1[ii, ii, ]^2 *
              dvars.deta[, ii]^2
          
          for (ii in 1:(M1 - 2 * extra$NOS.b))
            wz[, (2 * extra$NOS.b + ii)] <- apply(counter, 1, function(x) {
              (0.5) * sum(diag(temp1[, , x] %*% dTs[, , ii] %*% 
                                 temp1[, , x] %*% dTs[, , ii])) 
            }) * doffds.deta[, ii]^2
          
          pos.thetas <- combVGAMextra(x = 1:extra$M)[-c(1:extra$M), ]
          Der.dETAS <- cbind(dmeans.deta, dvars.deta, doffds.deta)
          
          for (ii in 1:nrow(pos.thetas)) {
            x.row <- pos.thetas[ii, ]
            wz[, extra$M + ii] <- apply(counter, 1, function(x){
              ((0.5) * sum(diag(temp1[, , 1] %*% der.Vars[[x.row[1]]] %*% 
                                  temp1[, , 1] %*% der.Vars[[x.row[2]]])) +
                 t(der.MUs[[x.row[1]]]) %*%
                 temp1[, , 1] %*% der.MUs[[x.row[2]]]) 
            })
            wz[, extra$M + ii] <- wz[, extra$M + ii] * 
              Der.dETAS[, x.row[1]] * Der.dETAS[, x.row[2]] 
          }
          
          c(w) * wz
          
        }), list( .trinorm = trinorm , 
                  .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                  .emean = emean , .evar = evar , .eoffd = eoffd  ))) )
  
  alo   
}


VARff.control <- function(save.weights = TRUE,
                          summary.HDEtest = FALSE, ...) {
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}