## misc useful functions

.cov.est <- function(bmat, var.epsilon, xx){
  m <- nrow(bmat)
  vc <- cov(bmat) - 1/m * var.epsilon * xx
  ## We may force it to be positive semi-definite. Right now just return vc
  return(vc)
}



## A convenience function to generate the combined covariate matrix
DataPrep <- function(GeneExp, CellProp, Demo){
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
  X0 <- cbind(CellProp, Demo0, Crossterms)
  ## create the long table for all covariates
  X <- cbind("ID"=rep(1:m, each=n),
             do.call(rbind, replicate(m,X0,simplify=FALSE)))
  ## create the long vector of Y
  Y <- as.vector(t(GeneExp))
  return(list(X=X,Y=Y))
}


## A wrapper function to compute EBLUP


#----------------------------------------------------------------------------#
# the function to calculate the chi-square type value. Similar as ols.eblup  #
#----------------------------------------------------------------------------#
hy.ols.blup.wrapper <- function(Des, Y, var.epsilon, number, random = c(1,2,3,4), vc, independent = F, trim.idx = NULL) {
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
    sigmabeta_i <-  lapply((1:m)[-trim.idx], function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
    sigmabeta <- var.epsilon * solve(Reduce("+", sigmabeta_i))

    ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% Y[(number[i]+1):(number[i+1])])
    XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])

    #######################################################################################
    ##### modify the trim,idx case at 6/14/2018: estimate fixed effect with trimed sunjects
    #######################################################################################
    betai <- lapply((1:m)[-trim.idx], function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
    betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
  }
  #betahat.back <- A %*% betahat

  ######EBLUP:
  gamma.hat <- lapply(1:m, function(i) lambda.hat * A %*% (ZY[[i]] - t(XZ[[i]]) %*% betahat - lambda.hat * ZZ[[i]] %*% cap[[i]] %*%
                                                             ZY[[i]] + lambda.hat * ZZ[[i]] %*% cap[[i]] %*% t(XZ[[i]]) %*% betahat))
  blup <- t(do.call(cbind, gamma.hat))
  ## assign colnames to blup
  colnames(blup) <- colnames(Des)[-1]

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
  #var.part3 = -2*var.part2
  var.eblup <- lapply(1:m, function(i) var.part1[[i]] - var.part2[[i]])

  ######################## refit the B matrix  ####################################################
  ######################## the eblup values #######################################################
  eta.stat <- unlist(lapply(1:m, function(i) {t(blup[i,]) %*% solve(var.eblup[[i]]) %*% blup[i,]}))
  eta.test <- lapply(1:m, function(i) {blup[i,]/sqrt(diag(var.eblup[[i]]))})
  eta.test <- t(do.call("cbind",eta.test))
  ## the covariance estimation
  cov = vc.hat*var.epsilon * lambda.hat
  rownames(cov) <- colnames(cov) <- colnames(Des)[-1]
  return(list(eta.stat = eta.stat, eta.test = eta.test,
              blup = blup, betahat = betahat, sigmabeta = sigmabeta,
              cov = cov, lambda.hat = lambda.hat))
}
