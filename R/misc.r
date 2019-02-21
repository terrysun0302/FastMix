## misc useful functions

## a robust matrix inverse function. Here min.cond.num refers to the
## threshold of the minimum condition number of X
rsolve <- function(X, min.cond.num=1e-6){
  o <- svd(X)
  ## replace smaller singular values with min.cond.num * the largest
  ## singular value of X
  lambda <- min.cond.num*max(o$d)
  dd <- pmax(o$d, lambda)
  ## below are some tricks to ensure rsolve() works for 1x1 matrix (a number).
  if(length(dd)==1){
    dd.inv <- 1/dd
  } else {
    dd.inv <- diag(1/dd)
  }
  return(o$v %*% dd.inv %*% t(o$u))
}

## a robust version of half inverse of a *symmetric, semi-definite* matrix
rhalfinv <- function(X, min.cond.num=1e-6) {
  X <- as.matrix(X)
  ## make sure X is a squared matrix
  if (nrow(X) != ncol(X)) stop("X must be a squared matrix!")
  if (!isSymmetric(X, tol=1e-6)) stop("X must be a symmetric matrix!")
  ## if X is a 1x1 matrix (number), just return sqrt(X)
  if (all(dim(X)==c(1,1))) {
    Xhalfinv <- 1/sqrt(pmax(X,0))
  } else {
    o <- svd(X); u <- o$u
    lambda <- min.cond.num*max(o$d)
    dd <- pmax(o$d, lambda)
    Xhalfinv <- u %*% diag(1/sqrt(dd)) %*% t(u)
  }
  return(Xhalfinv)
}


data_gen <- function(m, n, seed = 1, outlier = T, cor = 0, balance = T){
  ## generate some new simulated data
  n <- n
  set.seed(seed)
  CellProp <- cbind(Cell1=runif(n, .2, .3),
                    Cell2=runif(n, .2, .3),
                    Cell3=runif(n, .2, .3))
  rownames(CellProp) <- paste0("Subj", 1:n)
  Demo <- cbind(Severity=round(runif(n, 0, 10), 1),
                Sex=c(rep(0, n/2), rep(1, n/2)))
  rownames(Demo) <- paste0("Subj", 1:n)

  ## generate some true signals

  m <- m                               #genes
  bb <- c(1.5, 0, 0)                     #assoc. w. CellProp
  sigma2.u <- 1
  sigma2.e <- 0.6
  sigma2.alpha <- 1.2
  sigma2.err <- 1/16
  #sigma2.err <- 0.5 #i.i.d. noise
  ## mean interactions
  aa <- rep(0, 6); aa[1] <- .75

  Bmat <- t(replicate(m, c(bb, 0, 0, aa)))
  #Bmat <- t(replicate(m, c(1.5, 0.5, 0.5, -0.5, -0.5, 0.75, 0.5, 0.5, 0.5, -0.5, -0.5)))
  L <- ncol(Bmat)
  sdvec <- (c(rep(sigma2.u, 3), rep(sigma2.e,2), rep(sigma2.alpha, 6)))

  if(outlier == T){
    ### the location shift
    idx = c(1,4,7,10)
    str = sdvec[idx] * 3.5
    OutSig = matrix(NA, ncol = 4, nrow = m)
    for(i in 1:4){
      OutSig[,i] <- c(rep(0, m * 0.05 * (i-1)), rep(str[i], m * 0.05), rep(0, m*(1 - i * 0.05)))
    }
    if(balance == T) OutSig <- OutSig * (2*rbinom(length(OutSig), 1, .5) - 1) # approximately 0 mean
    else if (balance == F) OutSig <- OutSig * (2*rbinom(length(OutSig), 1, 0.9) - 1) # approximately 0 mean

    Bmat[,idx] = Bmat[,idx] + OutSig
    if(cor == 0){
      Gmat <- matrix(rnorm(L*m),m,L) %*% diag(sdvec)
      sigmamat = diag(sdvec^2)
      for(i in 1:4){
        Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]] <- Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]]
      }
    }
    else{
      cormat <- matrix(cor, ncol = 11, nrow = 11)
      diag(cormat) = 1
      sigmamat <- sweep(sweep(cormat, 1, sdvec, "*"), 2, sdvec, "*")
      a = chol(sigmamat)
      Gmat <- matrix(rnorm(L*m),m,L) %*% a
      for(i in 1:4){
        Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]] <- Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]]
      }
    }
  }
  else{
    Bmat = Bmat
    if(cor == 0){
      Gmat <- matrix(rnorm(L*m),m,L) %*% diag(sdvec)
      sigmamat = diag(sdvec^2)
    }
    else{
      cormat <- matrix(cor, ncol = 11, nrow = 11)
      diag(cormat) = 1
      sigmamat <- sweep(sweep(cormat, 1, sdvec, "*"), 2, sdvec, "*")
      a = chol(sigmamat)
      Gmat <- matrix(rnorm(L*m),m,L) %*% a
    }
  }

  ## generate the gene expression
  Err <- matrix(rnorm(m*n, sd=sqrt(sigma2.err)), nrow=m, ncol=n)
  colnames(Err) <- paste0("Subj", 1:n); rownames(Err) <- paste0("Gene", 1:m)

  ## Use DataPrep to get Xmat; Y is irrelevant
  Xmat <- t(DataPrep(Err, CellProp, Demo)$X[1:n, -1])

  ## the gene expression matrix
  GeneExp <- (Bmat + Gmat) %*% Xmat + Err

  return(list(GeneExp = GeneExp, CellProp = CellProp, Demo = Demo, Bmat = Bmat, Gmat = Gmat, sigmamat = sigmamat))
}

.cov.est <- function(bmat, var.epsilon, xx, m, coef){
  #m <- nrow(bmat)
  #vc <- diag(sqrt(coef)) %*% cov(bmat) %*% diag(sqrt(coef)) - 1/m * var.epsilon * xx
  vc <- cov(bmat)
  diag(vc) = diag(vc) * coef
  vc <- vc - 1/m * var.epsilon * xx
  ## We may force it to be positive semi-definite. Right now just return vc
  return(vc)
}

robust.cov.est <- function(bmat, var.epsilon, xx, m, robust){
  #m <- nrow(bmat)
  vc <- covRob(bmat, estim = robust)$cov - 1/m * var.epsilon * xx
  ## We may force it to be positive semi-definite. Right now just return vc
  return(vc)
}

## A convenience function to generate the combined covariate matrix
DataPrep <- function(GeneExp, CellProp, Demo, include.demo=TRUE){
  m <- nrow(GeneExp); n <- nrow(Demo)
  ## assign gene names if empty
  if (is.null(rownames(GeneExp))) rownames(GeneExp) <- paste0("Gene", 1:m)
  ## turn vectors into matrix
  CellProp <- as.matrix(CellProp); Demo <- as.matrix(Demo)
  ## assign colnames (variable names) if empty
  if (is.null(colnames(CellProp))) colnames(CellProp) <- paste0("Cell", 1:ncol(CellProp))
  if (is.null(colnames(Demo))) colnames(Demo) <- paste0("Var", 1:ncol(Demo))
  ## standardize all clinical covariates
  Demo0 <- scale(Demo)
  #Demo0 = Demo
  ## create interaction terms
  Crossterms.idx <- cbind(rep(1:ncol(CellProp), ncol(Demo0)),
                          rep(1:ncol(Demo0), each=ncol(CellProp)))
  rownames(Crossterms.idx) <- paste(colnames(CellProp)[Crossterms.idx[,1]],
                                    colnames(Demo0)[Crossterms.idx[,2]], sep=".")
  Crossterms <- sapply(1:nrow(Crossterms.idx), function(l) {
    k <- Crossterms.idx[l,1]; p <- Crossterms.idx[l,2]
    CellProp[,k] * Demo0[,p]
  }); colnames(Crossterms) <- rownames(Crossterms.idx)
  ## combine all variables
  if (include.demo==TRUE){
    X0 <- cbind(CellProp, Demo0, Crossterms)
  } else {
    X0 <- cbind(CellProp, Crossterms)
  }
  ## check the SVD and provide some warning messages
  ss <- svd(X0)$d
  if (min(ss)<1e-7) stop("The design matrix is numerically singular. Please consider: (a) include less cell types in the model, or (b) use 'include.demo=FALSE' to exclude the main effects of demographic covariates from the model.")
  ## create the long table for all covariates
  X <- cbind("ID"=rep(1:m, each=n),
             do.call(rbind, replicate(m,X0,simplify=FALSE)))
  ## create the long vector of Y
  Y <- as.vector(t(GeneExp))
  return(list(X=X,Y=Y))
}


#----------------------------------------------------------------------------#
# the function to calculate the chi-square type value. Similar as ols.eblup  #
#----------------------------------------------------------------------------#
hy.ols.blup.wrapper <- function(Des, Y, var.epsilon, number, random = random, vc, independent = F, trim.idx = NULL) {
  if(length(random) == 1){
    A = sqrt(max(vc, var.epsilon/100))
    Des.prime <- Des[,c(1+random)] * A # transformt the design matrix
    DZ.prime <- Des.prime
    vc.hat <- A^2

    lambda.hat <- 1/var.epsilon
    ########## estimation of the covariance matrix and the vakue of beta
    ZZ <- lapply(1:m, function(i) sum(DZ.prime[(number[i]+1):(number[i+1])]^2) )
    cap <-  lapply(1:m, function(i) {solve(diag(length(random)) + lambda.hat * ZZ[[i]])})
    XZ <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]),-1]) %*%  DZ.prime[(number[i]+1):(number[i+1])])
    XX <-  lapply(1:m, function(i) {t(Des[(number[i]+1):(number[i+1]),-1]) %*% Des[(number[i]+1):(number[i+1]),-1]})

    if(is.null(trim.idx)){
      sigmabeta_i <-  lapply((1:m), function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
      sigmabeta <- var.epsilon * solve(Reduce("+", sigmabeta_i))

      ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1])]) %*% Y[(number[i]+1):(number[i+1])])
      XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])

      #######################################################################################
      ##### modify the trim,idx case at 6/14/2018: estimate fixed effect with trimed sunjects
      #######################################################################################
      betai <- lapply((1:m), function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
      betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
    }
    else{
      sigmabeta_i <-  lapply((1:m)[trim.idx], function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
      sigmabeta <- var.epsilon * solve(Reduce("+", sigmabeta_i))

      ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1])]) %*% Y[(number[i]+1):(number[i+1])])
      XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])

      #######################################################################################
      ##### modify the trim,idx case at 6/14/2018: estimate fixed effect with trimed subjects
      #######################################################################################
      betai <- lapply((1:m)[trim.idx], function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
      betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
    }

    ######EBLUP:
    gamma.hat <- lapply(1:m, function(i) lambda.hat * A %*% (ZY[[i]] - t(XZ[[i]]) %*% betahat - lambda.hat * ZZ[[i]] %*% cap[[i]] %*%
                                                               ZY[[i]] + lambda.hat * ZZ[[i]] %*% cap[[i]] %*% t(XZ[[i]]) %*% betahat))
    blup <- t(do.call(cbind, gamma.hat))
    ## assign colnames to blup
    colnames(blup) <- colnames(Des)[-1][random]

    ################################################################################################
    ################## DEFINE SOME USEFUL QUANTITIES TO COMPUTE THE VARIANCE OF EBLUP ##############
    ################################################################################################
    yvar <- lapply(1:m, function(i) solve(ZZ[[i]] + diag(var.epsilon, length(random))))
    #beta.inverse <- lapply(1:m, function(i) yvar[[i]] %*% ZZ[[i]])
    #beta.inverse <- lapply(1:m, function(i) yvar[[i]] %*% XX[[i]])
    #beta.inverse2 <- lapply(1:m, function(i) yvar[[i]] %*% ZZ[[i]])
    ############# this is actually another way tp compute sigmabeta at line 88 #####################
    #varbeta <- solve(solve(A) %*% Reduce("+", beta.inverse) %*% solve(A))

    var.part1 <- lapply(1:m, function(i) A %*%  ZZ[[i]] %*% yvar[[i]] %*% A)
    var.part2 <- lapply(1:m, function(i) A %*% yvar[[i]] %*% t(XZ[[i]]) %*% sigmabeta %*%
                          XZ[[i]] %*% yvar[[i]] %*% A)
    var2.part1 <- lapply(1:m, function(i) ZZ[[i]] %*% yvar[[i]])
    var2.part2 <- lapply(1:m, function(i) yvar[[i]] %*% t(XZ[[i]]) %*% sigmabeta %*%
                           XZ[[i]] %*% yvar[[i]])
    var2.eblup <- lapply(1:m, function(i) var2.part1[[i]] - var2.part2[[i]])
    #var.part3 = -2*var.part2
    var.eblup <- lapply(1:m, function(i) var.part1[[i]] - var.part2[[i]])

    ######################## refit the B matrix  ####################################################
    ######################## the eblup values #######################################################
    eta.stat <- unlist(lapply(1:m, function(i) {t(blup[i,]) %*% solve(var.eblup[[i]]) %*% blup[i,]}))

    # ### recover the covariance matrix
    # eta.stat3 <- lapply(1:m, function(i) {a = eigen(solve(var2.eblup[[i]])); solve(A) %*% a$vectors %*% diag(sqrt(a$values)) %*% t(a$vectors) %*% blup[i,]})
    # eta.stat3 <- t(do.call("cbind",eta.stat3))
    #
    eta.stat2 <- lapply(1:m, function(i) {a = eigen(solve(var.eblup[[i]])); a$vectors * sqrt(a$values) * t(a$vectors) * blup[i,]})
    eta.stat2 <- t(do.call("cbind",eta.stat2))

    eta.stat3 <- lapply(1:m, function(i) {1/sqrt(diag(var.eblup[[i]])) * blup[i,]})
    eta.stat3 <- t(do.call("cbind",eta.stat3))

    eta.test <- lapply(1:m, function(i) {blup[i,]/sqrt((var.eblup[[i]]))})
    eta.test <- t(do.call("cbind",eta.test))
    ## the covariance estimation
    cov = vc.hat*var.epsilon * lambda.hat
    #rownames(cov) <- colnames(cov) <- colnames(Des)[-1][random]
    return(list(eta.stat = eta.stat, eta.stat2 = eta.stat2, eta.test = eta.test,
                blup = blup, betahat = betahat, sigmabeta = sigmabeta,
                cov = cov, lambda.hat = lambda.hat))

  }
  else{
    if(independent == F){
      a <- eigen(vc)
      a$values <- pmax(a$values, var.epsilon/100)
      A <- (a$vectors %*% diag(sqrt(a$values)) %*% solve(a$vectors))
      vc.hat <- (a$vectors %*% diag(a$values) %*% solve(a$vectors))

      Des.prime <- Des[,c(1+random)] %*% A # transformt the design matrix
      DZ.prime <- Des.prime
    }
    else{
      a <- pmax(diag(vc), var.epsilon/100)
      A <- diag(sqrt(a))
      vc.hat <- diag(a)

      Des.prime <- Des[,c(1+random)] %*% A # transformt the design matrix
      DZ.prime <- Des.prime
    }
  }
  ##########Step 5: estimate \lambda
  lambda.hat <- 1/var.epsilon

  ########## estimation of the covariance matrix and the vakue of beta
  ZZ <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% DZ.prime[(number[i]+1):(number[i+1]),] )
  cap <-  lapply(1:m, function(i) {solve(diag(length(random)) + lambda.hat * ZZ[[i]])})
  XZ <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]),-1]) %*%  DZ.prime[(number[i]+1):(number[i+1]),])
  XX <-  lapply(1:m, function(i) {t(Des[(number[i]+1):(number[i+1]),-1]) %*% Des[(number[i]+1):(number[i+1]),-1]})

  if(is.null(trim.idx)){
    sigmabeta_i <-  lapply((1:m), function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
    sigmabeta <- var.epsilon * solve(Reduce("+", sigmabeta_i))

    ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% Y[(number[i]+1):(number[i+1])])
    XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])

    #######################################################################################
    ##### modify the trim,idx case at 6/14/2018: estimate fixed effect with trimed sunjects
    #######################################################################################
    betai <- lapply((1:m), function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
    betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
  }
  else{
    sigmabeta_i <-  lapply((1:m)[trim.idx], function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
    sigmabeta <- var.epsilon * solve(Reduce("+", sigmabeta_i))

    ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% Y[(number[i]+1):(number[i+1])])
    XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])

    #######################################################################################
    ##### modify the trim,idx case at 6/14/2018: estimate fixed effect with trimed subjects
    #######################################################################################
    betai <- lapply((1:m)[trim.idx], function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
    betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
  }
  #betahat.back <- A %*% betahat

  ######EBLUP:
  gamma.hat <- lapply(1:m, function(i) lambda.hat * A %*% (ZY[[i]] - t(XZ[[i]]) %*% betahat - lambda.hat * ZZ[[i]] %*% cap[[i]] %*%
                                                             ZY[[i]] + lambda.hat * ZZ[[i]] %*% cap[[i]] %*% t(XZ[[i]]) %*% betahat))
  blup <- t(do.call(cbind, gamma.hat))
  ## assign colnames to blup
  colnames(blup) <- colnames(Des)[-1][random]

  ################################################################################################
  ################## DEFINE SOME USEFUL QUANTITIES TO COMPUTE THE VARIANCE OF EBLUP ##############
  ################################################################################################
  yvar <- lapply(1:m, function(i) solve(ZZ[[i]] + diag(var.epsilon, length(random))))
  #beta.inverse <- lapply(1:m, function(i) yvar[[i]] %*% ZZ[[i]])
  #beta.inverse <- lapply(1:m, function(i) yvar[[i]] %*% XX[[i]])
  #beta.inverse2 <- lapply(1:m, function(i) yvar[[i]] %*% ZZ[[i]])
  ############# this is actually another way tp compute sigmabeta at line 88 #####################
  #varbeta <- solve(solve(A) %*% Reduce("+", beta.inverse) %*% solve(A))

  var.part1 <- lapply(1:m, function(i) A %*%  ZZ[[i]] %*% yvar[[i]] %*% A)
  var.part2 <- lapply(1:m, function(i) A %*% yvar[[i]] %*% t(XZ[[i]]) %*% sigmabeta %*%
                                         XZ[[i]] %*% yvar[[i]] %*% A)
  var2.part1 <- lapply(1:m, function(i) ZZ[[i]] %*% yvar[[i]])
  var2.part2 <- lapply(1:m, function(i) yvar[[i]] %*% t(XZ[[i]]) %*% sigmabeta %*%
                        XZ[[i]] %*% yvar[[i]])
  var2.eblup <- lapply(1:m, function(i) var2.part1[[i]] - var2.part2[[i]])
  #var.part3 = -2*var.part2
  var.eblup <- lapply(1:m, function(i) var.part1[[i]] - var.part2[[i]])

  ######################## refit the B matrix  ####################################################
  ######################## the eblup values #######################################################
  eta.stat <- unlist(lapply(1:m, function(i) {t(blup[i,]) %*% solve(var.eblup[[i]]) %*% blup[i,]}))

  # ### recover the covariance matrix
  # eta.stat3 <- lapply(1:m, function(i) {a = eigen(solve(var2.eblup[[i]])); solve(A) %*% a$vectors %*% diag(sqrt(a$values)) %*% t(a$vectors) %*% blup[i,]})
  # eta.stat3 <- t(do.call("cbind",eta.stat3))
  #
  eta.stat2 <- lapply(1:m, function(i) {a = eigen(solve(var.eblup[[i]])); a$vectors %*% diag(sqrt(a$values)) %*% t(a$vectors) %*% blup[i,]})
  eta.stat2 <- t(do.call("cbind",eta.stat2))

  eta.stat3 <- lapply(1:m, function(i) {1/sqrt(diag(var.eblup[[i]])) * blup[i,]})
  eta.stat3 <- t(do.call("cbind",eta.stat3))

  eta.test <- lapply(1:m, function(i) {blup[i,]/sqrt(diag(var.eblup[[i]]))})
  eta.test <- t(do.call("cbind",eta.test))
  ## the covariance estimation
  cov = vc.hat*var.epsilon * lambda.hat
  rownames(cov) <- colnames(cov) <- colnames(Des)[-1][random]
  return(list(eta.stat = eta.stat, eta.stat2 = eta.stat2, eta.stat3 = eta.stat3, eta.test = eta.test,
              blup = blup, betahat = betahat, sigmabeta = sigmabeta,
              cov = cov, lambda.hat = lambda.hat))
}
