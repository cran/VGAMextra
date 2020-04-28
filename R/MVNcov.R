##########################################################################
# These functions are
# Copyright (C) 2014-2020 V. Miranda & T. Yee
# Auckland University of Technology & University of Auckland
# All rights reserved.
#
# Links renamed on Jan-2019 conforming with VGAM_1.1-0

MVNcov <- function(zero = c("var", "cov"),
                   lmean = "identitylink",
                   lvar  = "loglink",
                   lcov = "identitylink") {
  
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
  trinorm <- FALSE
  
  alo <-
    new("vglmff",
      blurb = c("Multivariate Normal distribution \n",
                "Links:   ",
                namesof("mean", lmean, earg = emean), ", ",
                namesof("variance", lvar, earg = evar), ", ",
                namesof("Error covariances", loffd, earg = eoffd)),
      
      
      
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
      }))),
      
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
        
        
        extra$derVs    <- dVs
        extra$derTs    <- dTs
        #extra$derMUs <- dMUs
        rm(dMUs)
        extra$der.MUs  <- der.MUs
        extra$der.Vars <- der.Vars
        
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
                .emean = emean , .evar = evar , .eoffd = eoffd ))),
      
      
      
      
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
        
      }), list( .trinorm = trinorm , 
                .lmean = lmean , .lvar = lvar , .loffd = loffd ,
                .emean = emean , .evar = evar , .eoffd = eoffd ))),
      
      
      
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
      
      
      
      vfamily = c("MVNcov"),
      
      
      
      deriv = eval(substitute(expression({
        nn <- nrow(eta)
        counter <- cbind(1:nn)
        M1 <- extra$n.means + ( extra$n.means *
                                  (extra$n.means + 1) )/2
        
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
        der.Vars  <- extra$der.Vars 
        dMUs      <- extra$derMUs
        der.MUs   <- extra$der.MUs
        
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



dmultinorm <- function(vec.x, vec.mean = c(0, 0),
                     mat.cov = c(1, 1, 0), log = FALSE) {
  
  if (!(is.vector(vec.x) || is.matrix(vec.x) ))
    stop("Entries must be vectors or matrices. Here, vec.x is a ",
         class(vec.x))
  
  
  if (!(is.vector(vec.mean) || is.matrix(vec.mean) ))
    stop("Entries must be vectors or matrices. Here, vec.mean is a ",
         class(vec.mean))
  
  if (!(is.vector(mat.cov) || is.matrix(mat.cov) ))
    stop("Entries must be vectors or matrices. Here, mat.cov is a ",
         class(mat.cov))

  if (length(c(mat.cov)) <= 2)
    stop("At least two responses handled. Else, refer to dnorm().")
  
  if (!length(nrow(vec.x))) {
    vec.x    <- rbind(vec.x)
    vec.mean <- matrix(vec.mean, nrow(vec.x), ncol(vec.x), byrow =  TRUE)
    mat.cov  <- matrix(mat.cov, nrow(vec.x),
                       ncol(vec.x) * (ncol(vec.x) + 1) / 2, byrow =  TRUE)
  } else {
    vec.mean <- matrix(c(vec.mean), nrow(vec.x), ncol(vec.x),
                       byrow = ( !length(dim(vec.mean)) ))
    mat.cov  <- matrix(c(mat.cov), nrow(vec.x),
                       ncol(vec.x) * (ncol(vec.x) + 1) / 2, 
                       byrow = ( !length(dim(mat.cov)) ))
  }
  
  nn <- nrow(vec.x)
  nos.int <- ncol(vec.x)
  count <- cbind(1:nn)
  
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  mat.cov <- apply(mat.cov, 1 , function(x) {
    int.mat <- matrix(0, nos.int, nos.int)
    diag(int.mat) <- x[1:nos.int]
    int.mat[row(int.mat) > col(int.mat)] <- x[-c(1:nos.int)]
    int.mat <- t(int.mat)
    int.mat[row(int.mat) > col(int.mat)] <- x[-c(1:nos.int)]
    int.mat
  })
  
  attr(mat.cov, "dim") <- c(nos.int, nos.int, nn)
  
  solve.mat.cov <- apply(count, 1, function(x) solve(mat.cov[, , x]))
  attr(solve.mat.cov, "dim") <- c(nos.int, nos.int, nn)
  
  logpdf <- apply(count, 1, function(x) {
    log( (det(2 * pi * mat.cov[, , x]) )^(-0.5) ) + 
       (-0.5) * (vec.x - vec.mean)[x, ] %*%
      solve.mat.cov[, , x] %*%  cbind(vec.x - vec.mean)[x, ]
  })
  
  if (log.arg) logpdf else exp(logpdf)
  
}


MVNcov.control <- function(save.weights = TRUE, 
                           summary.HDEtest = FALSE,...) {
  list(save.weights = save.weights, 
       summary.HDEtest = summary.HDEtest,...)
}
