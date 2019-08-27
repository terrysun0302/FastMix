## misc useful functions

#============================================================================================#
#========================= some functions for simulation study ==============================#
#============================================================================================#

# the data generation for main simulation
data_gen <- function(m, n, seed = 1, outlier = T, cor = 0, balance = T, weight = F){
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
  if(typeof(weight) == "logical"){
    if(weight == F){
      Err <- matrix(rnorm(m*n, sd=sqrt(sigma2.err)), nrow=m, ncol=n)
    }
  }
  else{ # i.d but not i.i.d
    Err <- matrix(rnorm(m*n, sd=sqrt(sigma2.err)), nrow=m, ncol=n)
    ### the cor structure
    weight = weight
    e = eigen(weight)
    w = e$vectors %*% sqrt(diag(e$values)) %*% t(e$vectors)
    Err = Err %*% w
  }
  colnames(Err) <- paste0("Subj", 1:n); rownames(Err) <- paste0("Gene", 1:m)

  ## Use DataPrep to get Xmat; Y is irrelevant
  Xmat <- t(DataPrep(Err, CellProp, Demo, w = diag(rep(1, n)))$X[1:n, -1])

  ## the gene expression matrix
  GeneExp <- (Bmat + Gmat) %*% Xmat + Err

  return(list(GeneExp = GeneExp, CellProp = CellProp, Demo = Demo, Bmat = Bmat, Gmat = Gmat, sigmamat = sigmamat, Err = Err))
}

# # the data generation for csSAM simulation
# ### generate the data for csSAM analysis. Only 1 binary Demo variable
# data_gen_csSAM <- function(m, n, seed = 1, outlier = T, cor = 0, balance = T, weight = F){
#   ## generate some new simulated data
#   n <- n
#   set.seed(seed)
#   CellProp <- cbind(Cell1=runif(n, .2, .3),
#                     Cell2=runif(n, .2, .3),
#                     Cell3=runif(n, .2, .3))
#   rownames(CellProp) <- paste0("Subj", 1:n)
#   Demo <- as.matrix(c(rep(0, n/2), rep(1, n/2)))
#   rownames(Demo) <- paste0("Subj", 1:n)
#
#   ## generate some true signals
#
#   m <- m                               #genes
#   bb <- c(1.5, 0, 0)                    #assoc. w. CellProp
#   sigma2.u <- 1
#   sigma2.e <- 0.6
#   sigma2.alpha <- 1.2
#   sigma2.err <- 1/16
#   #sigma2.err <- 0.5 #i.i.d. noise
#   ## mean interactions
#   aa <- rep(0, 2); aa[1] <- .75
#
#   #Bmat <- t(replicate(m, c(bb, 0, 0.5, aa)))
#   Bmat <- t(replicate(m, c(1.5, 0.5, 0.5, 0,0,0,0)))
#   #Bmat <- t(replicate(m, c(1.5, 0.5, 0.5, -0.5, -0.5, 0.75, 0.5, 0.5, 0.5, -0.5, -0.5)))
#   L <- ncol(Bmat)
#   sdvec <- (c(rep(sigma2.u, 3), rep(sigma2.alpha,4)))
#
#   if(outlier == T){
#     ### the location shift
#     idx = c(5,6)
#     str = sdvec[idx] * 3.5
#     OutSig = matrix(NA, ncol = 2, nrow = m)
#     for(i in 1:2){ #the first 250 genes are DGEs for idx 4; the 251 to 500 genes are DEGs for idx 5
#       OutSig[,i] <- c(rep(0, m * 0.05 * (i-1)), rep(str[i], m * 0.05), rep(0, m*(1 - i * 0.05)))
#     }
#     if(balance == T) OutSig <- OutSig * (2*rbinom(length(OutSig), 1, .5) - 1) # approximately 0 mean
#     else if (balance == F) OutSig <- OutSig * (2*rbinom(length(OutSig), 1, 0.9) - 1) # approximately 0 mean
#
#     Bmat[,idx] = Bmat[,idx] + OutSig
#     if(cor == 0){
#       Gmat <- matrix(rnorm(L*m),m,L) %*% diag(sdvec)
#       sigmamat = diag(sdvec^2)
#       for(i in 1:2){
#         Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]] <- Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]]
#       }
#     }
#     else{
#       cormat <- matrix(cor, ncol = 7, nrow = 7)
#       diag(cormat) = 1
#       sigmamat <- sweep(sweep(cormat, 1, sdvec, "*"), 2, sdvec, "*")
#       a = chol(sigmamat)
#       Gmat <- matrix(rnorm(L*m),m,L) %*% a
#       for(i in 1:2){
#         Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]] <- Gmat[(m * 0.05 * (i-1)+1):(m * 0.05 * i), idx[i]]
#       }
#     }
#   }
#   else{
#     Bmat = Bmat
#     if(cor == 0){
#       Gmat <- matrix(rnorm(L*m),m,L) %*% diag(sdvec)
#       sigmamat = diag(sdvec^2)
#     }
#     else{
#       cormat <- matrix(cor, ncol = 11, nrow = 11)
#       diag(cormat) = 1
#       sigmamat <- sweep(sweep(cormat, 1, sdvec, "*"), 2, sdvec, "*")
#       a = chol(sigmamat)
#       Gmat <- matrix(rnorm(L*m),m,L) %*% a
#     }
#   }
#
#   ## generate the gene expression
#   if(typeof(weight) == "logical"){
#     if(weight == F){
#       Err <- matrix(rnorm(m*n, sd=sqrt(sigma2.err)), nrow=m, ncol=n)
#     }
#   }
#   else{ # i.d but not i.i.d
#     Err <- matrix(rnorm(m*n, sd=sqrt(sigma2.err)), nrow=m, ncol=n)
#     ### the cor structure
#     weight = weight
#     e = eigen(weight)
#     w = e$vectors %*% sqrt(diag(e$values)) %*% t(e$vectors)
#     Err = Err %*% w
#   }
#   colnames(Err) <- paste0("Subj", 1:n); rownames(Err) <- paste0("Gene", 1:m)
#
#   ## Use DataPrep to get Xmat; Y is irrelevant
#   Xmat <- t(DataPrep(Err, CellProp, Demo, w = diag(rep(1, n)))$X[1:n, -1])
#
#   ## the gene expression matrix
#   GeneExp <- (Bmat + Gmat) %*% Xmat + Err
#
#   return(list(GeneExp = GeneExp, CellProp = CellProp, Demo = Demo, Bmat = Bmat, Gmat = Gmat, sigmamat = sigmamat, Err = Err))
# }


#============================================================================================#
#========================== general purposes helper functions ===============================#
#============================================================================================#
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
  if (!isSymmetric(X, tol=1e-6)) warning("X needs to be a symmetric matrix!")
  ## if X is a 1x1 matrix (number), just return sqrt(X)
  if (all(dim(X)==c(1,1))) {
    Xhalfinv <- 1/sqrt(pmax(X,0))
  } else {
    ## R's eigen has a symmetric switch
    o <- eigen(X, symmetric=TRUE); u <- o$vector
    lambda <- min.cond.num*max(o$values)
    dd <- pmax(o$values, lambda)
    Xhalfinv <- u %*% diag(1/sqrt(dd)) %*% t(u)
  }
  return(Xhalfinv)
}

cov.est <- function(bmat, var.epsilon, xx, m, coef){
  #m <- nrow(bmat)
  #vc <- diag(sqrt(coef)) %*% cov(bmat) %*% diag(sqrt(coef)) - 1/m * var.epsilon * xx
  vc <- cov(as.matrix(bmat))
  diag(vc) = diag(vc) * coef
  vc <- vc - 1/m * var.epsilon * xx
  ## We may force it to be positive semi-definite. Right now just return vc
  return(as.matrix(vc))
}

robust.cov.est <- function(bmat, var.epsilon, xx, m, robust){
  #m <- nrow(bmat)
  vc <- covRob(bmat, estim = robust)$cov - 1/m * var.epsilon * xx
  ## We may force it to be positive semi-definite. Right now just return vc
  return(as.matrix(vc))
}

#============================================================================================#
#======================= some helper functions in FastMix procedur===========================#
#============================================================================================#

## A convenience function to generate the combined covariate matrix to fit FastMix
DataPrep <- function(GeneExp, CellProp, Demo, include.demo=TRUE, w="iid"){ # here, w is the weight function
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
  ## by default, the weighting matrix is a diagonal matrix,
  ## representing the i.i.d. data
  if (w=="iid") w <- diag(nrow(X0))

  ### weighted sample
  X0 = w %*% X0
  GeneExp = GeneExp %*% w

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

## A convenience function to generate the design matrix for the response prediction

#==================================================================================================================#
#============================ new added functions for fix effect bias correction 04/15/2019 =======================#
#==================================================================================================================#

## Sfunc1 is a function that takes mu, R, and L (dimensions), and
## returns S (the conditional mean), dS (the derivative), and some
## intermediate results (T_R(mu), B_R(mu), and H_R(mu), see Eqn(4)).
Sfunc1 <- function(mu, R, L, ngrid=500){
  dx <- 2*R/ngrid; mygrid <- seq(-R, R, dx)
  phi.grid <- dnorm(mygrid-mu)
  Fchi2.grid <- pchisq(R^2-mygrid^2, L-1)
  h.grid <- Fchi2.grid*phi.grid
  Bmu <- sum(h.grid)*dx
  Tmu <- sum(h.grid*mygrid)*dx
  Hmu <- sum(h.grid*(mygrid^2))*dx
  return(c(Tmu=Tmu, Bmu=Bmu, Hmu=Hmu,
           S=Tmu/Bmu, dS=(Hmu*Bmu - Tmu^2)/Bmu^2))
}

## This is the univariate version of the inverse of S, using
## Newton-Raphson. mu1.trim: trimmed mean of mu1 (the first dimension).
Sinv1 <- function(mu1.trim, R, L, ngrid=500, tol=1e-4, max.iter=10){
  k <- 0; err <- Inf; mu1.now=mu1.trim
  while ( ( abs(err) > tol) & (k < max.iter) ) {
    rk <- Sfunc1(mu1.now,R,L,ngrid=ngrid)
    mu1.next <- mu1.now - (rk[["S"]]-mu1.trim)/rk[["dS"]]
    k <- k+1; err <- mu1.next-mu1.now
    mu1.now <- mu1.next
  }
  if ( abs(err) > tol ) warning("Maximum number of iteration reached but the error is still greater than the tolerance level.")
  return(list(mu1=mu1.now, niter=k, err=err))
}

## This is the conditional variance of X
sigmaR2 <- function(trim, L) {
  sigmaR2 <- pchisq(qchisq(1-trim, df=L), df=L+2) / (1-trim)
  return(sigmaR2)
}


## this is a wrapper for estimating the trimmed mean with bias
## correction. X is an nxL dimensional matrix (so the transpose is
## required). Sigma is the covariance matrix of X, which should be
## estimated prior to mean estimation.
RobustMeanEst <- function(X, Sigma=diag(ncol(X)), trim=0.5, tol=1e-2, max.iter=10,...){
  L <- ncol(X); S.half.inv <- rhalfinv(Sigma); S.half <- rsolve(S.half.inv)
  ## 1. standardize X
  Xtilde <- X %*% S.half.inv
  ## 2. initial estimation and trimming
  mus.now=colMeans(Xtilde)
  ## start the main loop
  k <- 0; err <- Inf;
  while ( ( abs(err) > tol) & (k < max.iter) ) {
    ## (a) the centered data
    Z <- sweep(Xtilde, 2, mus.now)
    ## (b) trimming by Euclidean length squared of Z
    ll <- sqrt(rowSums(Z^2))
    R <- as.numeric(quantile(ll, 1-trim)); idx <- which(ll <= R)
    ## (c). compute the trimmed mean, and decompose it into length/direction
    if(L == 1)  {
      Mvec.trim <- mean(Z[idx,])
    }
    else {
      Mvec.trim <- colMeans(Z[idx,])
    }
    M.trim <- sqrt(sum(Mvec.trim^2)); v <- Mvec.trim/M.trim
    ## (d) bias-correction
    o <- Sinv1(M.trim, R, L=L, ...)
    err <- o$mu1    #magnitude of the kth update
    mus.next <- mus.now + err * v
    mus.now <- mus.next
    k <- k+1
  }
  ## 4. rotate mus.now to mu.hat
  mu.est <- drop(S.half %*% mus.now)
  ## 5. compute the variance inflation coefficient
  kappa <- err/M.trim; V.inflation <- kappa^2 * sigmaR2(trim, L)/(1-trim)
  return(list(mu.est=mu.est, niter=k, err=err, R=R, V.inflation=V.inflation))
}
#================================================== end ===========================================================#
#==================================================================================================================#

#----------------------------------------------------------------------------#
# the function to calculate the chi-square type value. Similar as ols.eblup  #
#----------------------------------------------------------------------------#
hy.ols.blup.wrapper <- function(Des, Y, var.epsilon, number, random = random, vc, independent = F, trim.idx = NULL, min.cond.num=1e-6,
                                bias_term) {
  N <- length(Y)
  m <- length(unique(Des[,1]))
  n <- N/m

  vc <- as.matrix(vc)                   #works for 1x1 matrix
  if(independent == F){
    a <- eigen(vc, symmetric=TRUE)
    a$values <- pmax(a$values, var.epsilon/100)
    ## deal with length(random)==1
    if (length(a$values)==1) {
      L <- as.matrix(sqrt(a$values)); L2 <- as.matrix(a$values)
    } else {
      L <-  diag(sqrt(a$values)); L2 <-  diag(a$values)
    }
    A <- a$vectors %*% L %*% t(a$vectors)
    vc.hat <- a$vectors %*% L2 %*% t(a$vectors)
    Des.prime <- Des[,c(1+random)] %*% A # transformt the design matrix
    DZ.prime <- Des.prime
  } else {
    #the independent case
    vc = diag(vc)
    a <- eigen(vc, symmetric=TRUE)
    a$values <- pmax(a$values, var.epsilon/100)
    ## deal with length(random)==1
    if (length(a$values)==1) {
      L <- as.matrix(sqrt(a$values)); L2 <- as.matrix(a$values)
    } else {
      L <-  diag(sqrt(a$values)); L2 <-  diag(a$values)
    }
    A <- a$vectors %*% L %*% t(a$vectors)
    vc.hat <- a$vectors %*% L2 %*% t(a$vectors)
    Des.prime <- Des[,c(1+random)] %*% A # transformt the design matrix
    DZ.prime <- Des.prime
  }

  ##########Step 5: estimate \lambda
  lambda.hat <- 1/var.epsilon

  ########## estimation of the covariance matrix and the value of beta
  ZZ <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% DZ.prime[(number[i]+1):(number[i+1]),] )
  cap <-  lapply(1:m, function(i) {rsolve(diag(length(random)) + lambda.hat * ZZ[[i]])})
  XZ <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]),-1]) %*%  DZ.prime[(number[i]+1):(number[i+1]),])
  XX <-  lapply(1:m, function(i) {t(Des[(number[i]+1):(number[i+1]),-1]) %*% Des[(number[i]+1):(number[i+1]),-1]})

  if(is.null(trim.idx)){
    sigmabeta_i <-  lapply((1:m), function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
    sigmabeta <- var.epsilon * rsolve(Reduce("+", sigmabeta_i))

    ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% Y[(number[i]+1):(number[i+1])])
    XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])

    #######################################################################################
    ##### modify the trim,idx case at 6/14/2018: estimate fixed effect with trimed subjects
    #######################################################################################
    betai <- lapply((1:m), function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
    betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
  }
  else{
    sigmabeta_i <-  lapply((1:m)[trim.idx], function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
    sigmabeta <- var.epsilon * rsolve(Reduce("+", sigmabeta_i), min.cond.num=min.cond.num)
    ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% Y[(number[i]+1):(number[i+1])])
    XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])

    #######################################################################################
    ##### modify the trim,idx case at 6/14/2018: estimate fixed effect with trimed subjects
    #######################################################################################
    betai <- lapply((1:m)[trim.idx], function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
    betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
  }
  #betahat.back <- A %*% betahat

  ### 04/15/2019: new added step: bias correction
  betahat[random] = betahat[random] + bias_term

  ######EBLUP:
  gamma.hat <- lapply(1:m, function(i) lambda.hat * A %*% (ZY[[i]] - t(XZ[[i]]) %*% betahat - lambda.hat * ZZ[[i]] %*% cap[[i]] %*% ZY[[i]] + lambda.hat * ZZ[[i]] %*% cap[[i]] %*% t(XZ[[i]]) %*% betahat))
  blup <- t(do.call(cbind, gamma.hat))
  ## assign colnames to blup
  colnames(blup) <- colnames(Des)[-1][random]

  ################################################################################################
  ################## DEFINE SOME USEFUL QUANTITIES TO COMPUTE THE VARIANCE OF EBLUP ##############
  ################################################################################################
  yvar <- lapply(1:m, function(i) rsolve(ZZ[[i]] + diag(var.epsilon, length(random)), min.cond.num=min.cond.num))
  #beta.inverse <- lapply(1:m, function(i) yvar[[i]] %*% ZZ[[i]])
  #beta.inverse <- lapply(1:m, function(i) yvar[[i]] %*% XX[[i]])
  #beta.inverse2 <- lapply(1:m, function(i) yvar[[i]] %*% ZZ[[i]])
  ############# this is actually another way tp compute sigmabeta at line 88 #####################
  #varbeta <- rsolve(rsolve(A) %*% Reduce("+", beta.inverse) %*% rsolve(A))
  var.part1 <- lapply(1:m, function(i) A %*%  ZZ[[i]] %*% yvar[[i]] %*% A)
  var.part2 <- lapply(1:m, function(i) A %*% yvar[[i]] %*% t(XZ[[i]]) %*% sigmabeta %*% XZ[[i]] %*% yvar[[i]] %*% A)
  var2.part1 <- lapply(1:m, function(i) ZZ[[i]] %*% yvar[[i]])
  var2.part2 <- lapply(1:m, function(i) yvar[[i]] %*% t(XZ[[i]]) %*% sigmabeta %*% XZ[[i]] %*% yvar[[i]])
  var2.eblup <- lapply(1:m, function(i) var2.part1[[i]] - var2.part2[[i]])
  #var.part3 = -2*var.part2
  var.eblup <- lapply(1:m, function(i) var.part1[[i]] - var.part2[[i]])

  ######################## refit the B matrix  ####################################################
  ######################## the eblup values #######################################################

  ########## eta.stat is used in re.pvalue
  eta.stat <- unlist(lapply(1:m, function(i) {t(blup[i,]) %*% rsolve(var.eblup[[i]], min.cond.num=min.cond.num) %*% blup[i,]}))

  # ### recover the covariance matrix
  # eta.stat3 <- lapply(1:m, function(i) {a = eigen(rsolve(var2.eblup[[i]])); rsolve(A) %*% a$vectors %*% diag(sqrt(a$values)) %*% t(a$vectors) %*% blup[i,]})
  # eta.stat3 <- t(do.call("cbind",eta.stat3))
  #
  ## eta.stat2 <- lapply(1:m, function(i) {a = eigen(rsolve(var.eblup[[i]])); a$vectors %*% diag(sqrt(a$values)) %*% t(a$vectors) %*% blup[i,]})

  ########## eta.stat2 are used in computing re.ind.pvalue ##########
  eta.stat2 <- lapply(1:m, function(i) {rhalfinv(var.eblup[[i]], min.cond.num=min.cond.num) %*% blup[i,]})
  eta.stat2 <- t(do.call("cbind",eta.stat2))

  ########## eta.stat3 are used in Anderson-Darling test
  eta.stat3 <- lapply(1:m, function(i) {1/sqrt(diag(as.matrix(var.eblup[[i]]))) * blup[i,]})
  eta.stat3 <- t(do.call("cbind",eta.stat3))

  eta.test <- lapply(1:m, function(i) {blup[i,]/sqrt(diag(as.matrix(var.eblup[[i]])))})
  eta.test <- t(do.call("cbind",eta.test))
  ## the covariance estimation
  cov = vc.hat*var.epsilon * lambda.hat
  rownames(cov) <- colnames(cov) <- colnames(Des)[-1][random]
  ## also output var.eblup.mean for debugging purposes
  var.eblup.mean <- Reduce("+", var.eblup)/length(var.eblup)
  return(list(eta.stat = eta.stat, eta.stat2 = eta.stat2, eta.stat3 = eta.stat3,
              eta.test = eta.test, blup = blup, betahat = betahat,
              var.eblup.mean=var.eblup.mean, sigmabeta = sigmabeta,
              cov = cov, lambda.hat = lambda.hat))
}




#============================================================================================#
#======================= some helper functions for prediction score =========================#
#============================================================================================#

# see the vignettes for an concrete example

### function 1
#=========================================================================================#
#========== A function to create pseudo data sets used in the score calculation ==========#
#=================================== input ===============================================#
# GeneExp: test GeneExp data                                                              #
# CellProp: test CellProp data                                                            #
# Demo: test Demo data except for the response vector; Null if no such information        #
# train_response: the response vector from training data                                  #
# include.demo: this option should be the same as the option in FastMix. Default = T      #
# w: the weight covariance matrix among subjects                                          #
#=========================================================================================#
###  A function to generate the test data
### # here, w is the weight function. Response_idx is the idx of response variable in the Demo matrix
### Demo here is a vertor/matrix that contains demo information except for the Response vector, because this is the
### variable we want to classify.
DataPrep_test <- function(GeneExp, CellProp, Demo, train_response, include.demo=TRUE, w="iid"){
  m <- nrow(GeneExp); n <- ncol(GeneExp)
  ## assign gene names if empty
  if (is.null(rownames(GeneExp))) rownames(GeneExp) <- paste0("Gene", 1:m)
  ## turn vectors into matrix
  CellProp <- as.matrix(CellProp)
  ## assign colnames (variable names) if empty
  if (is.null(colnames(CellProp))) colnames(CellProp) <- paste0("Cell", 1:ncol(CellProp))

  ### Demo may be Null
  if(is.null(Demo)) {
    Demo0 = c()
  }
  else{
    Demo <- as.matrix(Demo)
    if(is.null(colnames(Demo))) {
      colnames(Demo) <- paste0("Var", 1:ncol(Demo))
    }

    ## standardize all clinical covariates
    Demo0 <- scale(Demo)
  }

  ### by definition,
  response_impute = sort(unique(scale(train_response)))
  #============================================================#
  #======== here we replace the Demo0 as another vector =======#
  #============================================================#

  Demo0_0 = cbind(Demo0, Res = rep(response_impute[1], n))
  Demo0_1 = cbind(Demo0, Res = rep(response_impute[2], n))
  # if(Response == 0){
  #   ### Generate the data assuming Response = 0
  #   ## actually if Demo is binary, then after scaling, the subjects with Demo = 0 should be <0.
  #   R_value = unique(Demo0[,Response_idx]) ### the unique 2 values
  #   R_value = R_value[R_value < 0] ### the value we use to impute the data
  #   Demo0[Demo0[,Response_idx] != R_value, Response_idx] = R_value ### replace the Demo0
  # }
  # else if(Response == 1){
  #   R_value = unique(Demo0[,Response_idx])
  #   R_value = R_value[R_value > 0] ### the value we use to impute the data
  #   Demo0[Demo0[,Response_idx] != R_value, Response_idx] = R_value ### replace the Demo0
  # }
  #Demo0 = Demo
  ## create interaction terms



  #============================================#
  #======= construct Data0 ====================#
  #============================================#
  Crossterms.idx <- cbind(rep(1:ncol(CellProp), ncol(Demo0_0)),
                          rep(1:ncol(Demo0_0), each=ncol(CellProp)))
  rownames(Crossterms.idx) <- paste(colnames(CellProp)[Crossterms.idx[,1]],
                                    colnames(Demo0_0)[Crossterms.idx[,2]], sep=".")
  Crossterms <- sapply(1:nrow(Crossterms.idx), function(l) {
    k <- Crossterms.idx[l,1]; p <- Crossterms.idx[l,2]
    CellProp[,k] * Demo0_0[,p]
  }); colnames(Crossterms) <- rownames(Crossterms.idx)
  ## combine all variables
  if (include.demo==TRUE){
    X0 <- cbind(CellProp, Demo0_0, Crossterms)
  } else {
    X0 <- cbind(CellProp, Crossterms)
  }

  ## by default, the weighting matrix is a diagonal matrix,
  ## representing the i.i.d. data
  if (w=="iid") w <- diag(nrow(X0))
  ### weighted sample
  X0 = w %*% X0
  GeneExp = GeneExp %*% w

  ## we ignore the svd because all Demo0 is the same.
  #ss <- svd(X0)$d
  #if (min(ss)<1e-7) stop("The design matrix is numerically singular. Please consider: (a) include less cell types in the model, or (b) use 'include.demo=FALSE' to exclude the main effects of demographic covariates from the model.")
  ## create the long table for all covariates
  X <- cbind("ID"=rep(1:m, each=n),
             do.call(rbind, replicate(m,X0,simplify=FALSE)))
  ## create the long vector of Y
  Y <- as.vector(t(GeneExp))

  Data_0 = list(X=X,Y=Y)


  #============================================#
  #======= construct Data1 ====================#
  #============================================#
  Crossterms.idx <- cbind(rep(1:ncol(CellProp), ncol(Demo0_1)),
                          rep(1:ncol(Demo0_1), each=ncol(CellProp)))
  rownames(Crossterms.idx) <- paste(colnames(CellProp)[Crossterms.idx[,1]],
                                    colnames(Demo0_1)[Crossterms.idx[,2]], sep=".")
  Crossterms <- sapply(1:nrow(Crossterms.idx), function(l) {
    k <- Crossterms.idx[l,1]; p <- Crossterms.idx[l,2]
    CellProp[,k] * Demo0_1[,p]
  }); colnames(Crossterms) <- rownames(Crossterms.idx)
  ## combine all variables
  if (include.demo==TRUE){
    X0 <- cbind(CellProp, Demo0_1, Crossterms)
  } else {
    X0 <- cbind(CellProp, Crossterms)
  }

  ### weighted sample
  X0 = w %*% X0
  GeneExp = GeneExp %*% w

  ## we ignore the svd because all Demo0 is the same.
  #ss <- svd(X0)$d
  #if (min(ss)<1e-7) stop("The design matrix is numerically singular. Please consider: (a) include less cell types in the model, or (b) use 'include.demo=FALSE' to exclude the main effects of demographic covariates from the model.")
  ## create the long table for all covariates
  X <- cbind("ID"=rep(1:m, each=n),
             do.call(rbind, replicate(m,X0,simplify=FALSE)))
  ## create the long vector of Y
  Y <- as.vector(t(GeneExp))

  Data_1 = list(X=X,Y=Y)

  return(list(Data_0 = Data_0, Data_1 = Data_1))
}



### function 3
#=========================================================================================#
#========== A function to calculated summary scores used in classification ===============#
#=================================== input ===============================================#
#=========================================================================================#
# mod1: the fitted mod using FastMix from training data                                   #
# GeneExp: test GeneExp data                                                              #
# CellProp: test CellProp data                                                            #
# Demo: test Demo data. We require it to be binary here                                   #
# idx: the index for covariates without Demo interaction                                  #
#=========================================================================================#

### m, gene number; n, subjects number;
### interaction_index, the index in design matrix for the interaction between cell and Demo
### Response_idx is the index of all temrs about the resposne variable except for the response*cell interaction terms
score_func = function(mod1, Data_0, Data_1, Response_interaction_index, Response_idx){
  ### gene-level coefficient
  gene_coef = mod1$beta.mat

  m_train = dim(mod1$beta.mat)[1]
  p_train = dim(mod1$beta.mat)[2]
  m = length(unique(Data_0$X[,1]))
  n = dim(Data_0$X)[1]/m
  p = dim(Data_0$X)[2]-1 # -1 for ID column

  if(p != p_train){
    stop("The number of covariates is different, please check!")
  }
  if(m_train != m){
    stop("The number of genes is different, please check!")
  }

  ### more than one sybjects in the test data
  if(n > 1) {
    ### because the data is stacked by ID, we also need to argument the coef_matrix
    #expand_coef = gene_coef[rep(1:nrow(gene_coef), each = n), ] # n is the number of subjects
    mu_0 = gene_coef %*% t(as.matrix(Data_0$X[1:n,-1])) ## the first col is ID, and the design matrix doesn't change over genes
    mu_1 = gene_coef %*% t(as.matrix(Data_1$X[1:n,-1])) ## the first col is ID, and the design matrix doesn't change over genes

    Y = matrix(Data_0$Y, ncol = n, nrow = m, byrow = T) ### transform the Y vector into the m by n matrix
    res = Y - mu_0   ### Y - \mu_0
    delta = mu_1 - mu_0


    #=======================================================#
    #====== step1 : generate a single summary score ========#
    #=======================================================#
    delta_scale = apply(delta, 2, function(x) sqrt(sum(x^2)))
    ### calculate the
    single_score = t(res) %*% delta
    single_score = diag(single_score) / (delta_scale)


    #=======================================================#
    #====== step2: generate multi summary scores ===========#
    #=======================================================#
    multi_score = res * delta ### m by n

    #=======================================================#
    #====== step3: generate sparse summary scores ==========#
    #=======================================================#
    ### the difference would be the calculation of delta vector

    ### we first identify DEGs in the interaction directions
    l = length(Response_interaction_index)
    sparse_coef = fix_coef = matrix(0, ncol = l+length(Response_idx), nrow = m)
    fix_coef[,1:length(Response_idx)] = mod1$beta.mat[,Response_idx] ### the coef for Demo variable
    ### fill the coef for interaction
    for(i in 1:l){
      fix_coef[,i+length(Response_idx)] = mod1$fixed.results[Response_interaction_index[i],1]
      idx = which((mod1$re.ind.pvalue)[,Response_interaction_index[i]] < 0.05)
      sparse_coef[idx,i+length(Response_idx)] = mod1$beta.mat[idx,i] - mod1$fixed.results[Response_interaction_index[i],1]
      ### only random effects
    }

    ### calculate the sparse score
    delta_sparse_fix = fix_coef %*% as.matrix(t(Data_1$X[1:n,c(Response_idx, Response_interaction_index)+1] -
                                                  Data_0$X[1:n,c(Response_idx, Response_interaction_index)+1]))
    delta_sparse_random = sparse_coef %*%  as.matrix(t(Data_1$X[1:n,c(Response_idx, Response_interaction_index)+1] -
                                                         Data_0$X[1:n,c(Response_idx, Response_interaction_index)+1]))

    delta_sparse_scale= apply(delta_sparse_fix+delta_sparse_random, 2, function(x) sqrt(sum(x^2)))
    single_sparse_score = t(res) %*% (delta_sparse_fix+delta_sparse_random)
    single_sparse_score = diag(single_sparse_score) / (delta_sparse_scale)

    #=======================================================#
    #====== step4: generate sparse multiple scores =========#
    #=======================================================#

    multi_sparse_score = res * (delta_sparse_fix+delta_sparse_random) ### m by n


    result = list(single_score = single_score, single_sparse_score = single_sparse_score, multi_score = multi_score, multi_sparse_score = multi_sparse_score)
    return(result)
  }
  else if(n == 1) {
    ### because the data is stacked by ID, we also need to argument the coef_matrix
    #expand_coef = gene_coef[rep(1:nrow(gene_coef), each = n), ] # n is the number of subjects
    mu_0 = gene_coef %*% as.matrix(Data_0$X[1:n,-1]) ## the first col is ID, and the design matrix doesn't change over genes
    mu_1 = gene_coef %*% as.matrix(Data_1$X[1:n,-1]) ## the first col is ID, and the design matrix doesn't change over genes

    Y = matrix(Data_0$Y, ncol = n, nrow = m, byrow = T) ### transform the Y vector into the m by n matrix
    res = Y - mu_0   ### Y - \mu_0
    delta = mu_1 - mu_0


    #=======================================================#
    #====== step1 : generate a single summary score ========#
    #=======================================================#
    delta_scale = apply(delta, 2, function(x) sqrt(sum(x^2)))
    ### calculate the
    single_score = t(res) %*% delta
    single_score = diag(single_score) / (delta_scale)


    #=======================================================#
    #====== step2: generate multi summary scores ===========#
    #=======================================================#
    multi_score = res * delta ### m by n

    #=======================================================#
    #====== step3: generate sparse summary scores ==========#
    #=======================================================#
    ### the difference would be the calculation of delta vector

    ### we first identify DEGs in the interaction directions
    l = length(Response_interaction_index)
    sparse_coef = fix_coef = matrix(0, ncol = l+length(Response_idx), nrow = m)
    fix_coef[,1:length(Response_idx)] = mod1$beta.mat[,Response_idx] ### the coef for Demo variable
    ### fill the coef for interaction
    for(i in 1:l){
      fix_coef[,i+length(Response_idx)] = mod1$fixed.results[Response_interaction_index[i],1]
      idx = which((mod1$re.ind.pvalue)[,Response_interaction_index[i]] < 0.05)
      sparse_coef[idx,i+length(Response_idx)] = mod1$beta.mat[idx,i] - mod1$fixed.results[Response_interaction_index[i],1]
      ### only random effects
    }

    ### calculate the sparse score
    delta_sparse_fix = fix_coef %*% as.matrix(Data_1$X[1:n,c(Response_idx, Response_interaction_index)+1] -
                                                Data_0$X[1:n,c(Response_idx, Response_interaction_index)+1])
    delta_sparse_random = sparse_coef %*%  as.matrix(Data_1$X[1:n,c(Response_idx, Response_interaction_index)+1] -
                                                       Data_0$X[1:n,c(Response_idx, Response_interaction_index)+1])

    delta_sparse_scale= apply(delta_sparse_fix+delta_sparse_random, 2, function(x) sqrt(sum(x^2)))
    single_sparse_score = t(res) %*% (delta_sparse_fix+delta_sparse_random)
    single_sparse_score = diag(single_sparse_score) / (delta_sparse_scale)


    #=======================================================#
    #====== step4: generate sparse multiple scores =========#
    #=======================================================#

    multi_sparse_score = res * (delta_sparse_fix+delta_sparse_random) ### m by n


    result = list(single_score = single_score, single_sparse_score = single_sparse_score, multi_score = multi_score, multi_sparse_score = multi_sparse_score)

    return(result)
  }

}
