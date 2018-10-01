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
########################################################################################################
################################ the function of proposed method #######################################
########################################################################################################

ols.eblup.trim <- function(Des, Y, random = "all", independent = T, trim = 0.2, type = "F"){
  N <- length(Y) # number of total observations
  p <- dim(Des)[2]-1 #ID is not a covariate but a label
  if(length(random) == 1 && random == "all"){
    random <- 1:p
  }
  else random <- random
  DZ <- Des[, c(1,random+1)] #the dense form of Z matrix, not a block diagonal one
  dz <- as.data.frame(DZ) #a data.frame used for some particular functions

  #----------------------------------------------------------------------------#
  # ols based step to initially calculate the covariance of random effects     #
  #----------------------------------------------------------------------------#
  m <- length(unique(Des[,1])) # the # of subjects
  ID <- sort(unique(Des[,1]))
  qprime <- m*p #the number of parameters

  ###calculate the var.epsilon
  number <- dz %>% group_by(ID) %>% summarize(num = length(ID)) %>% rbind(c(0,0), .) %>% as.matrix #define the # of obs for each subject
  number <- cumsum(number[,2]) #idx for each subject
  res <- lapply(1:m, function(i) sum(lm.fit(x = Des[(number[i]+1):(number[i+1]),-1], y = Y[(number[i]+1):(number[i+1])])$residuals^2))
  var.epsilon <- sum(unlist(res)) / (N - qprime) #sigma.e

  ####################################### calculate the OLS-based covariance matrix#################################################
  coef.fix <- lm.fit(Des[,-1],Y)$coeff
  ### Y - fixed effect
  new.Y <- Y - Des[,-1] %*% coef.fix
  ### estimate variance of random effects
  ols <- lapply(1:m, function(i) lm.fit(as.matrix(Des[(number[i]+1):(number[i+1]),(1 +random)]), y = new.Y[(number[i]+1):(number[i+1])])$coeff)
  ols <- do.call(rbind, ols)
  ## covbeta <- cov(ols)

  xx <- lapply(1:m, function(i) solve(t(Des[(number[i]+1):(number[i+1]),(1 +random)]) %*% Des[(number[i]+1):(number[i+1]),(1 +random)]))
  xx <- Reduce("+", xx)
  vc <- .cov.est(ols, var.epsilon, xx, m)

  #----------------------------------------------------------------------------#
  # trimming step based on the chi-square type statistics                      #
  #----------------------------------------------------------------------------#
  initialfit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc, independent = independent)
  ### we only use the middle subjects
  eta.stat <- initialfit$eta.stat
  trimdf <- length(random)
  if(type == "F"){
    trim.idx <- which(eta.stat %in% eta.stat[order(eta.stat)][1:(m*(1 - trim))])
    ### refit the empirical variance by trimed subjects
    vc.trim <- .cov.est(ols[trim.idx,], var.epsilon, xx)
    simf <- trimdf*(m-1)/(m-trimdf)*rf(1000000, df1 = trimdf, df2 = m - trimdf)
    lambda.quan <-  trimdf /mean((simf[order(simf)][1:(1000000*(1-trim))]))
    vc.refit <- vc.trim * lambda.quan

    if(sum(diag(vc.refit) < 0) > 0){
      testit <- function() warning("potential problem by trimming. Please try a smaller trim coefficient.", call. = FALSE)
    }

    refit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc.refit, independent = independent, trim.idx = NULL)

    p.unadjust <- 1 - pf(refit$eta.stat*(m-trimdf)/trimdf/(m-1), df1 = trimdf, df2 = m - trimdf)
    p.ind.unadjust <- 2*(1 - pt(abs(refit$eta.test), df = m-trimdf))
  }

  else if(type == "chi"){
    trim.idx <- which(eta.stat %in% eta.stat[order(eta.stat)][1:(m*(1 - trim))])
    ### refit the empirical variance by trimed subjects
    vc.trim <- .cov.est(ols[trim.idx,], var.epsilon, xx, m)
    ### I just use simulated value instead of theoretical value
    simchi <- rchisq(1000000, df = trimdf)
    lambda.quan <-  trimdf /mean((simchi[order(simchi)][1:(1000000*(1-trim))]))
    vc.refit <- vc.trim * lambda.quan

    if(sum(diag(vc.refit) < 0) > 0){
      warning("Fitting problems caused by trimming. Please try a smaller trim coefficient.", call. = FALSE)
    }

    refit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc.refit, independent = independent, trim.idx = NULL)

    p.unadjust <- 1 - pchisq(refit$eta.stat, df = trimdf)
    p.ind.unadjust <- 2*(1 - pnorm(abs(refit$eta.test)))
  }

  #----------------------------------------------------------------------------#
  # the calculation of final outputs                                           #
  #----------------------------------------------------------------------------#
  betamat <- matrix(rep(refit$betahat, m), nrow = m, byrow = T)
  colnames(betamat) <- colnames(Des)[-1]
  betamat[,random] <- refit$blup + betamat[,random]
  t.fixed <- refit$betahat/sqrt(diag(refit$sigmabeta))
  vc.new <- sqrt(pmax(diag(vc.refit),0)*var.epsilon * refit$lambda.hat)[random]
  VC <- data.frame(sqrt(var.epsilon), vc.new)
  Yhat <- lapply(1:m, function(i) Des[(number[i]+1):(number[i+1]), -1] %*% betamat[i,])
  Yhat <- do.call(rbind, Yhat)
  colnames(VC) <- c("sigma.e", "sigma.gamma")

  list(beta.hat=refit$betahat, beta.mat=betamat, Yhat=Yhat, sigma.beta=refit$sigmabeta, VC=VC,
       cov = refit$cov, t.fixed=t.fixed,
       p.unadjust = p.unadjust, eta = refit$eta.stat,p.ind.unadjust=p.ind.unadjust)
}


########## the main wrapper for deconvolution problem ##########
FastMix <- function(GeneExp, CellProp, Clinical, random="all", ...){
  gnames <- rownames(GeneExp); m <- nrow(GeneExp)
  if (is.null(gnames)) {
    rownames(GeneExp) <- gnames <- paste0("Gene", 1:m)
  }
  Data2 <- DataPrep(GeneExp, CellProp, Demo)
  # L <- ncol(Data2$X)-1
  # if (random=="all") random <- 1:L
  mod <- ols.eblup.trim(Des=Data2$X, Y=Data2$Y, random=random, ...)
  ## replace mod$Yhat by a matrix
  GeneExp.fitted <- matrix(mod$Yhat, nrow=m, byrow=TRUE)
  rownames(GeneExp.fitted) <- gnames
  colnames(GeneExp.fitted) <- colnames(GeneExp)
  mod$GeneExp.fitted <- GeneExp.fitted
  mod$Yhat <- NULL
  ## assign gene names to part of the results x
  rownames(mod$beta.mat) <- gnames
  rownames(mod$p.ind.unadjust) <- gnames
  names(mod$p.unadjust) <- gnames
  names(mod$eta) <- gnames
  return(mod)
}
