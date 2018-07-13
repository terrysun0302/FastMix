#******************************************************************************#
# fast LMER function                                                           #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  Des            the design matrix ordered by time subject by subject.        #
#                 First column should be ID range from 1 to the sample size;   #
#                 the rest columns are covariates                              #
#                                                                              #
#                                                                              #
#  Y              the longitudinal outcome ordered by time subject by subject. #
#                                                                              #
#                                                                              #
#  random         the user selected covariates with random effects.            #
#                                                                              #
#  independent    if the random effects are corrected or not. T or F           #
#                                                                              #
#  trim           the trimming percentage when accounting for outliers.        #
#                                                                              #
#  Outputs                                                                     #
#                                                                              #
#  Returns a list                                                              #
#                                                                              #
#  beta.hat       the fix effect estimation                                    #
#                                                                              #
#  beta.mat       individual coefficient estimation                            #
#                                                                              #
#  Yhat           fitted response                                              #
#                                                                              #
#  sigma.beta     the covariance estimation of fixed effect                    #
#                                                                              #
#  VC             variance component estimation. The first column is the one   #
#                 for common random error. The second column is the one for    #
#                 random effects.                                              #
#                                                                              #
#  t.fixed        the t value for fixed effects                                #
#                                                                              #
#  eta            the chi sqiare type statsitics                               #
#                                                                              #
#  p.unadjust     the overall p-value for outlier detection                    #
#                                                                              #
#  p.ind.unadjust the individual p-value for outlier detection for each        #
#                 random effect.                                               #
#******************************************************************************#
library(dplyr)
########################################################################################################
################################ the function without considering outlier ##############################
########################################################################################################
ols.eblup <- function(Des, Y, random = c(1,2,3,4), independent = T, type = "F"){
  N <- length(Y) # number of total observations
  p <- dim(Des)[2]-1 #ID is not a covariate but a label
  DZ <- Des[, c(1,random+1)] #the dense form of Z matrix, not a block diagonal one
  dz <- as.data.frame(DZ) #a data.frame used for some particular functions

  m <- length(unique(Des[,1])) # the # of subjects
  ID <- sort(unique(Des[,1]))
  qprime <- m*p #the number of parameters

  ###calculate the var.epsilon
  number <- dz %>% group_by(ID) %>% summarize(num = length(ID)) %>% rbind(c(0,0), .) %>% as.matrix #define the # of obs for each subject
  number <- cumsum(number[,2]) #idx for each subject
  res <- lapply(1:m, function(i) sum(lm.fit(x = Des[(number[i]+1):(number[i+1]),-1], y = Y[(number[i]+1):(number[i+1])])$residuals^2))
  var.epsilon <- sum(unlist(res)) / (N - qprime) #sigma.e

  ####################################### calculate the OLS-based covariance matrix#################################################

  #ols <- lapply(1:m, function(i) lm.fit(Des[(number[i]+1):(number[i+1]),-1], y = Y[(number[i]+1):(number[i+1])])$coeff)
  #ols <- do.call(rbind, ols)
  #covbeta <- cov(ols)

  coef.fix <- lm.fit(Des[,-1],Y)$coeff
  ### Y - fixed effect
  new.Y <- Y - Des[,-1] %*% coef.fix
  ### estimate variance of random effects
  ols <- lapply(1:m, function(i) lm.fit(Des[(number[i]+1):(number[i+1]),(1 +random)], y = new.Y[(number[i]+1):(number[i+1])])$coeff)
  ols <- do.call(rbind, ols)
  covbeta <- cov(ols)

  xx <- lapply(1:m, function(i) solve(t(Des[(number[i]+1):(number[i+1]),(1 +random)]) %*% Des[(number[i]+1):(number[i+1]),(1 +random)]))
  xx <- Reduce("+", xx)
  vc <- covbeta - 1/m * var.epsilon * xx

  if(independent == F){
    a <- eigen(vc)
    a$values <- pmax(a$values, var.epsilon/100)
    A <- (a$vectors %*% diag(sqrt(a$values)) %*% solve(a$vectors))
    vc.hat <- (a$vectors %*% diag(a$values) %*% solve(a$vectors))

    Des.prime <- Des[,c(random+1)] %*% A # transformt the design matrix
    DZ.prime <- Des.prime
  }
  else{
    a <- pmax(diag(vc), var.epsilon/100)
    #A <- diag(sd)
    A <- diag(sqrt(a))
    vc.hat <- diag(a)

    Des.prime <- Des[,c(random+1)] %*% A # transformt the design matrix
    DZ.prime <- Des.prime
  }


  ##########Step 5: estimate \lambda
  lambda.hat <- 1/var.epsilon

  ########## estimation of the covariance matrix and the vakue of beta
  ZZ <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% DZ.prime[(number[i]+1):(number[i+1]),] )
  cap <-  lapply(1:m, function(i) {solve(diag(length(random)) + lambda.hat * ZZ[[i]])})
  XZ <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]),-1]) %*%  DZ.prime[(number[i]+1):(number[i+1]),])
  XX <-  lapply(1:m, function(i) {t(Des[(number[i]+1):(number[i+1]),-1]) %*% Des[(number[i]+1):(number[i+1]),-1]})
  XWX <-
    sigmabeta_i <-  lapply(1:m, function(i) {XX[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% t(XZ[[i]])})
  sigmabeta <- var.epsilon * solve(Reduce("+", sigmabeta_i))

  ZY <- lapply(1:m, function(i) t(DZ.prime[(number[i]+1):(number[i+1]),]) %*% Y[(number[i]+1):(number[i+1])])
  XY <-  lapply(1:m, function(i) t(Des[(number[i]+1):(number[i+1]), -1]) %*% Y[(number[i]+1):(number[i+1])])
  betai <- lapply(1:m, function(i) XY[[i]] - lambda.hat * XZ[[i]] %*% cap[[i]] %*% ZY[[i]])
  betahat <- 1/var.epsilon * sigmabeta %*% Reduce("+", betai)
  #betahat.back <- A %*% betahat

  ######EBLUP:
  gamma.hat <- lapply(1:m, function(i) lambda.hat * A %*% (ZY[[i]] - t(XZ[[i]]) %*% betahat - lambda.hat * ZZ[[i]] %*% cap[[i]] %*%
                                                             ZY[[i]] + lambda.hat * ZZ[[i]] %*% cap[[i]] %*% t(XZ[[i]]) %*% betahat))
  blup <- t(do.call(cbind, gamma.hat))
  ################################################################################################
  ################## DEFINE SOME USEFUL QUANTITIES TO COMPUTE THE VARIANCE OF EBLUP ##############
  ################################################################################################
  yvar <- lapply(1:m, function(i) solve(ZZ[[i]] + diag(var.epsilon, length(random))))
  beta.inverse <- lapply(1:m, function(i) yvar[[i]] %*% ZZ[[i]])
  ############# this is actually another way tp compute sigmabeta at line 88 #####################
  varbeta <- solve(solve(A) %*% Reduce("+", beta.inverse) %*% solve(A))

  var.part1 <- lapply(1:m, function(i) A %*%  ZZ[[i]] %*% yvar[[i]] %*% A)
  var.part2 <- lapply(1:m, function(i) A %*% beta.inverse[[i]] %*% solve(A) %*% varbeta %*%
                        solve(A) %*% t(beta.inverse[[i]]) %*% A)
  #var.part3 = -2*var.part2

  var.eblup <- lapply(1:m, function(i) var.part1[[i]] - var.part2[[i]])
  ######################### let's only save the diagonal part for p-value #########################
  #var.eblup.diag <- lapply(1:m, function(i) diag(var.eblup[[i]]))
  #var.eblup.diag <- do.call("rbind", var.eblup.diag)

  ######################## the eblup values #######################################################
  eta.stat <- unlist(lapply(1:m, function(i) {t(blup[i,]) %*% solve(var.eblup[[i]]) %*% blup[i,]}))
  eta.test <- lapply(1:m, function(i) {blup[i,]/sqrt(diag(var.eblup[[i]]))})
  eta.test <- t(do.call("cbind",eta.test))

  cut.value <- qchisq(0.95, df = length(random))
  out <- eta.stat > cut.value
  if(type == "chi"){
    p.unadjust <- 1 - pchisq(eta.stat, df = length(random))
    p.ind.unadjust <- 2*(1 - pnorm(abs(eta.test)))
  }
  else if(type == "F"){
    p.unadjust <- 1 - pf(eta.stat*(m-p)/p/(m-1), df1 = p, df2 = m - p)
    p.ind.unadjust <- 2*(1 - pt(abs(eta.test), df = m-p))
  }
  p.adjust <- 1 - pchisq(eta.stat, df = length(random))^m
  ####################### the p-values for eblup ##################################################
  #p.eblup <- lapply(1:m, function(i) {t<-blup[i,]/sqrt(var.eblup.diag[i,]); 1 - pnorm(abs(t))})
  #p.eblup <- do.call("rbind", p.eblup)
  #t.eblup <- p.eblup < 0.05

  betamat <- matrix(rep(betahat, m), nrow = m, byrow = T)
  betamat[,random] <- blup + betamat[,random]
  t.fixed <- betahat/sqrt(diag(sigmabeta))
  vc.new <- sqrt(diag(vc.hat*var.epsilon * lambda.hat))[random]
  VC <- data.frame(sqrt(var.epsilon), vc.new)
  Yhat <- lapply(1:m, function(i) Des[(number[i]+1):(number[i+1]), -1] %*% betamat[i,])
  Yhat <- do.call(rbind, Yhat)
  colnames(VC) <- c("Sigma.e", "Sigma.gamma")

  list(beta.hat=betahat, beta.mat=betamat, Yhat=Yhat, sigma.beta=sigmabeta, VC=VC, t.fixed=t.fixed,
       cov = vc.hat*var.epsilon * lambda.hat,
       #p.eblup = p.eblup, t.eblup = t.eblup, blup = blup, sd.blup = sqrt(var.eblup.diag)
       p.adjust = p.adjust, p.unadjust = p.unadjust, p.ind.unadjust=p.ind.unadjust)
}