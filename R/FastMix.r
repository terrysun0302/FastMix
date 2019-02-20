### 10/22/2018 FastMix
#******************************************************************************#
# ols.eblup.trim function                                                      #
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
#  trim           the trimming percentage when accounting for outliers when    #
#                 robust = "FastMix"; default value is 0.5.                    #
#                                                                              #
#  robust         robust estimation of the covariance matrix of random effects.#
#                 "FASLE" for non-robust estimation; "mcd" for the MCD         #
#                 algorithm of Rousseeuw and Van Driessen; "weighted" for the  #
#                 Reweighted MCD; "donostah" for the Donoho-Stahel projection  #
#                 based estimator; "pairwiseQC" for the orthogonalized quadrant#
#                 correlation pairwise estimator. All these algorithms come    #
#                 from the R package "robust".                                 #
#                 "FastMix" is the proposed trimming method.                   #
#  test           the method to test DEGs. "1" is Gaussian mixture model; "2"  #
#                 is Anderson-darling normal test. The default value is "1"    #
#                                                                              #
#  trim.fix       If only consider trimmed subjects in fix effect estiamtion.  #
#                                                                              #
#  Outputs                                                                     #
#                                                                              #
#  Returns a list                                                              #
#                                                                              #
#  fixed.results   fix effects estimation and inference                        #
#                                                                              #
#  beta.mat       individual coefficient estimation                            #
#                                                                              #
#  Yhat           fitted response                                              #
#                                                                              #
#  sigma.beta     the covariance estimation of fixed effects                   #
#                                                                              #
#  VC             variance component estimation. The first column is the one   #
#                 for common random error. The second column is the one for    #
#                 random effects.                                              #
#                                                                              #
#  t.fixed        the t value for fixed effects                                #
#                                                                              #
#  eta            the chi square type statsitics                               #
#                                                                              #
#  re.pvalue      the overall p-value for outlier detection                    #
#                                                                              #
#  re.ind.pvalue the individual p-value for outlier detection for each        #
#                 random effect.                                               #
#                                                                              #
#  out_idx        the potential covariates with outliers when robust           #
#                 = "FastMix. It is NULL when robust != "FastMix"              #
#******************************************************************************#

########################################################################################################
################################ the function of proposed method #######################################
########################################################################################################

ols.eblup.trim <- function(Des, Y, random = "all", independent = F, trim = 0.5, robust = "FastMix", test = 1, trim.fix = TRUE){
  N <- length(Y) # number of total observations
  ## Exclude the first column, ID, because it is just a label
  covariates <- colnames(Des)[-1]; p <- length(covariates)
  if(length(random) == 1 && random == "all"){
    random <- 1:p
  } else {
    random <- random
  }
  DZ <- Des[, c(1,random+1)] #the dense form of Z matrix, not a block diagonal one
  dz <- as.data.frame(DZ) #a data.frame used for some particular functions

  colnames(dz)[1] <- "ID"
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

  xx <- lapply(1:m, function(i) rsolve(t(Des[(number[i]+1):(number[i+1]),(1 +random)]) %*% Des[(number[i]+1):(number[i+1]),(1 +random)]))
  XX <- Reduce("+", xx)
  coef_vector = rep(1, length(random))

  ### the case without robust estimation
  norm_idx = c()
  p_random = length(random)

  trim.idx = NULL

  if(p_random > 1){ # when the number of random effects is larger than 1
    if(robust == FALSE) {
      vc.refit <- .cov.est(ols, var.epsilon, XX, m, coef = coef_vector)
      norandom <- names(which(diag(vc.refit) <= 0))
      if(length(norandom) > 0){
        warning(paste0("Some covariates designated with random effects (",
                       paste(norandom, collapse=", "),
                       ") have zero empirical variance."))
      }
    }
    ### robust estimation by existing method from package "robust"
    else if(robust == "mcd") vc.refit <- robust.cov.est(ols, var.epsilon, XX, m, robust)
    else if(robust == "weighted") vc.refit <- robust.cov.est(ols, var.epsilon, XX, m, robust)
    else if(robust == "donostah") vc.refit <- robust.cov.est(ols, var.epsilon, XX, m, robust)
    else if(robust == "pairwiseQC") vc.refit <- robust.cov.est(ols, var.epsilon, XX, m, robust)

    ### porposed robust estimation
    else if(robust == "FastMix") {
      vc <- .cov.est(ols, var.epsilon, XX, m, coef = coef_vector)
      norandom <- names(which(diag(vc) <= 0))
      if(length(norandom) > 0){
        warning(paste0("Some covariates designated with random effects (",
                       paste(norandom, collapse=", "),
                       ") have zero empirical variance."))
      }
      #----------------------------------------------------------------------------#
      # trimming step based on the chi-square type statistics                      #
      #----------------------------------------------------------------------------#

      ### 10/11/2018 new added parts
      initialfit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc, independent = independent, trim.idx = trim.idx)
      B_cov = lapply(1:m, function(i) vc + var.epsilon * xx[[i]])
      B_cov_inv_half = lapply(1:m, function(i) {eig = eigen(rsolve(B_cov[[i]])); eig$vectors %*% sqrt(diag(eig$values)) %*% t(eig$vectors)})
      Norm_B = lapply(1:m, function(i) ols[i, ] %*% B_cov_inv_half[[i]])
      Norm_B = do.call(rbind,Norm_B)

      ### 1/20/2019  used to do test : test 1: gaussian mixed model test; test 2: anderson darling test
      if(test == 1){
        #norm_test = apply(initialfit$eta.stat2,2,function(x) shapiro.test(x)$p)
        norm_test = apply(initialfit$eta.stat3[,diag(vc) > 0],2,function(a)
          summary(Mclust(a, x = mclustBIC(a, verbose = F), verbose = F), parameters = TRUE)$G)
        norm_idx = rep(NA, p_random)
        norm_idx[diag(vc) > 0] = norm_test > 1
        norm_idx[diag(vc) <= 0] = 0
      }
      else{
        #norm_test = apply(initialfit$eta.stat3,2,function(x) shapiro.test(x)$p): this test has sample size limitation
        norm_test = apply(initialfit$eta.stat3,2,function(x) ad.test(x)$p) #: this test has sample size limitation
        norm_idx = (norm_test < 0.05) & (diag(vc) > 0)
      }
      trimdf <- sum(norm_idx)

      if(trimdf > 0){
        chi_stat = lapply(1:m, function(i) sum(Norm_B[i, norm_idx == T]^2))
        chi_stat = unlist(chi_stat)
        ## we remove extremely large chi_stat by amount 'trim' (by default=0.5)
        trim.idx <- which(chi_stat %in% chi_stat[order(chi_stat)][1:(m*(1 - trim))])
        ## truncated chi_stat follows a truncated chisq distribution;
        ## for now we use MC method to calculate lambda_alpha (Eq. 20 in the manuscript)
        simchi <- rchisq(1000000, df = trimdf)
        lambda.quan <-  trimdf /mean((simchi[order(simchi)][1:(1000000*(1-trim))]))
        if(trimdf == 1) {
          mult = var(Norm_B[trim.idx,norm_idx == T])*lambda.quan
        }
        else{
          mult = diag(cov(Norm_B[trim.idx,norm_idx == T]))*lambda.quan
        }
        coef_vector = c()
        coef_vector[norm_idx == T] = mult
        coef_vector[norm_idx == F] = 1
        vc.refit <- .cov.est(ols, var.epsilon, XX, m, coef = coef_vector)
      }
      else{
        vc.refit = vc
      }
      # if(sum(diag(vc.refit) <= 0) > 0){
      #   warning("Fitting problems caused by trimming. Please try a smaller trim coefficient.")
      # }
    }

    ### the option for fix effect
    if(trim.fix == FALSE){
      refit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc.refit, independent = independent, trim.idx = NULL)
    } else {
      refit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc.refit, independent = independent, trim.idx = trim.idx)
    }
    re.pvalue <- 1 - pchisq(refit$eta.stat, df = p_random)

    re.ind.pvalue <- 2*(1 - pnorm(abs(refit$eta.stat2)))
    # if(test == 1) re.ind.pvalue <- 2*(1 - pnorm(abs(refit$eta.stat2)))
    # else{re.ind.pvalue <- 2*(1 - pnorm(abs(refit$eta.stat3)))}
    colnames(re.ind.pvalue) <- covariates[random]
  }
  else{
    if(robust == FALSE) {
      vc.refit <- .cov.est(ols, var.epsilon, XX, m, coef = coef_vector)
      norandom <- names(which(diag(vc.refit) <= 0))
      if(length(norandom) > 0){
        warning(paste0("Some covariates designated with random effects (",
                       paste(norandom, collapse=", "),
                       ") have zero empirical variance."))
      }
    }

    ### robust estimation by existing method from package "robust"
    else if(robust == "mcd") stop("This method only works with more than 1 random effect.")
    else if(robust == "weighted") stop("This method only works with more than 1 random effect.")
    else if(robust == "donostah") stop("This method only works with more than 1 random effect.")
    else if(robust == "pairwiseQC") stop("This method only works with more than 1 random effect.")
    ### porposed robust estimation
    else if(robust == "FastMix") {
      vc <- .cov.est(ols, var.epsilon, XX, m, coef = coef_vector)
      norandom <- names(which(diag(vc) <= 0))
      if(length(norandom) > 0){
        warning(paste0("Some covariates designated with random effects (",
                       paste(norandom, collapse=", "),
                       ") have zero empirical variance."))
      }
      #----------------------------------------------------------------------------#
      # trimming step based on the chi-square type statistics                      #
      #----------------------------------------------------------------------------#

      ### 10/11/2018 new added parts
      initialfit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc, independent = independent)
      B_cov = lapply(1:m, function(i) vc + var.epsilon * xx[[i]])
      B_cov_inv_half = lapply(1:m, function(i) {eig = eigen(rsolve(B_cov[[i]])); eig$vectors %*% sqrt((eig$values)) %*% t(eig$vectors)})
      Norm_B = lapply(1:m, function(i) ols[i, ] %*% B_cov_inv_half[[i]])
      Norm_B = do.call(rbind,Norm_B)

      ### 1/20/2019  used to do test : test 1: gaussian mixed model test; test 2: anderson darling test
      if(test == 1){
        #norm_test = apply(initialfit$eta.stat2,2,function(x) shapiro.test(x)$p)
        norm_test = apply(initialfit$eta.stat3,2,function(a) summary(Mclust(a, x = mclustBIC(a, verbose = F), verbose = F), parameters = TRUE)$G)
        norm_idx = norm_test > 1
      }
      else{
        #norm_test = apply(initialfit$eta.stat3,2,function(x) shapiro.test(x)$p): this test has sample size limitation
        norm_test = apply(initialfit$eta.stat3,2,function(x) ad.test(x)$p) #: this test has sample size limitation
        norm_idx = norm_test < 0.05
      }
      trimdf <- sum(norm_idx)

      if(trimdf > 0){
        chi_stat = lapply(1:m, function(i) sum(Norm_B[i, norm_idx == T]^2))
        chi_stat = unlist(chi_stat)
        trim.idx <- which(chi_stat %in% chi_stat[order(chi_stat)][1:(m*(1 - trim))])
        simchi <- rchisq(1000000, df = trimdf)
        lambda.quan <-  trimdf /mean((simchi[order(simchi)][1:(1000000*(1-trim))]))

        mult = var(Norm_B[trim.idx,1])*lambda.quan

        coef_vector = c(mult)
        vc.refit <- .cov.est(ols, var.epsilon, XX, m, coef = coef_vector)
      }
      else{
        vc.refit = vc
      }
    }

    ### the option for fix effect
    if(trim.fix == FALSE){
      refit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc.refit, independent = independent)
    }

    else if(trim.fix == TRUE){
      refit <- hy.ols.blup.wrapper(Des, Y, var.epsilon, number, random = random, vc = vc.refit, independent = independent, trim.idx = trim.idx)
    }
    re.pvalue <- 1 - pchisq(refit$eta.stat, df = p_random)
    re.ind.pvalue <- 2*(1 - pnorm(abs(refit$eta.stat2)))
    colnames(re.ind.pvalue) <- covariates[random]
  }


  #----------------------------------------------------------------------------#
  # the calculation of final outputs                                           #
  #----------------------------------------------------------------------------#
  betamat <- matrix(rep(refit$betahat, m), nrow = m, byrow = T)
  colnames(betamat) <- colnames(Des)[-1]
  betamat[,random] <- refit$blup + betamat[,random]
  t.fixed <- drop(refit$betahat/sqrt(diag(as.matrix(refit$sigmabeta))))
  vc.new <- sqrt(pmax(diag(vc.refit),0)*var.epsilon * refit$lambda.hat)
  VC <- data.frame(sqrt(var.epsilon), vc.new)
  Yhat <- lapply(1:m, function(i) Des[(number[i]+1):(number[i+1]), -1] %*% betamat[i,])
  Yhat <- do.call(rbind, Yhat)
  colnames(VC) <- c("sigma.e", "sigma.gamma")

  ## [11/21/2018 Xing] inference for the fixed effects
  ## for now, the DF is sample size - regressors. This needs to be
  ## changed in the near future
  fixed.df <- min(m,N/m)-p-1
  fixed.p <- 2*pt(abs(drop(t.fixed)), df=fixed.df, lower.tail=FALSE)
  fixed.p.adj <- p.adjust(fixed.p, method="BH")
  fixed.results <- cbind(betahat=drop(refit$betahat), tstat=t.fixed,
                         p.value=fixed.p, p.adj=fixed.p.adj)

  list(fixed.results=fixed.results, beta.mat=betamat, Yhat=Yhat,
       sigma.beta=refit$sigmabeta, VC=VC, cov = refit$cov,
       re.pvalue = re.pvalue, eta = refit$eta.stat,
       re.ind.pvalue=re.ind.pvalue, out_idx = norm_idx)
}


########## the main wrapper for deconvolution problem ##########
FastMix <- function(GeneExp, CellProp, Demo, random="all", include.demo=TRUE, ...){
  gnames <- rownames(GeneExp); m <- nrow(GeneExp)
  if (is.null(gnames)) {
    rownames(GeneExp) <- gnames <- paste0("Gene", 1:m)
  }
  Data2 <- DataPrep(GeneExp, CellProp, Demo, include.demo=include.demo)
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
  rownames(mod$re.ind.pvalue) <- gnames
  names(mod$re.pvalue) <- gnames
  names(mod$eta) <- gnames
  return(mod)
}
