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

###  A function to generate the test data. Here, w is the weight
### function. Response_idx is the idx of response variable in the Demo
### matrix Demo here is a vertor/matrix that contains demo information
### except for the Response vector, because this is the variable we
### want to classify.
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
  if (all(w=="iid")) w <- diag(nrow(X0))
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

### m, gene number; n, subjects number; interaction_index, the index
### in design matrix for the interaction between cell and Demo
### Response_idx is the single index of the main response variable.
score_func = function(mod1, Data_0, Data_1, Response_interaction_index, Response_idx, sig.level=0.05){
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
  if (length(Response_idx) !=1) stop("You must select a single Response variable.")
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
    union_index = c()
    l = length(Response_interaction_index)
    sparse_coef = fix_coef = matrix(0, ncol = l+length(Response_idx), nrow = m)
    fix_coef[,1:length(Response_idx)] = mod1$fixed.results[Response_idx,1] ### the coef for Demo variable
    sparse_coef[,1:length(Response_idx)] = mod1$beta.mat[,Response_idx]-mod1$fixed.results[Response_idx,1]
    ### fill the coef for interaction
    for(i in 1:l){
      fix_coef[,i+length(Response_idx)] = mod1$fixed.results[Response_interaction_index[i],1]
      idx = which((mod1$re.ind.pvalue)[,Response_interaction_index[i]] < sig.level)
      union_index = union(union_index, idx)

      ## the reason for this "-" is because beta.mat is the individual
      ## coefficient (fix+ranodm)
      sparse_coef[idx,i+length(Response_idx)] = mod1$beta.mat[idx,Response_interaction_index[i]] - mod1$fixed.results[Response_interaction_index[i],1]
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

    multi_sparse_score = rbind(colSums(res * (delta_sparse_fix)),
                                        (res * delta_sparse_random)[union_index,]) ### m by n


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
    union_index = c()
    l = length(Response_interaction_index)
    sparse_coef = fix_coef = matrix(0, ncol = l+length(Response_idx), nrow = m)
    fix_coef[,1:length(Response_idx)] = mod1$fixed.results[Response_idx,1] ### the coef for Demo variable
    sparse_coef[,1:length(Response_idx)] = mod1$beta.mat[,Response_idx]-mod1$fixed.results[Response_idx,1]
    ### fill the coef for interaction
    for(i in 1:l){
      fix_coef[,i+length(Response_idx)] = mod1$fixed.results[Response_interaction_index[i],1]
      idx = which((mod1$re.ind.pvalue)[,Response_interaction_index[i]] < sig.level)
      union_index = union(union_index, idx)

      ## the eason for this "-" is because beta.mat is the individual coefficient (fix+ranodm)
      sparse_coef[idx,i+length(Response_idx)] = mod1$beta.mat[idx,Response_interaction_index[i]] -
        mod1$fixed.results[Response_interaction_index[i],1]
      ### only random effects
    }

    ### calculate the sparse score
    delta_sparse_fix = fix_coef %*% as.matrix((Data_1$X[1:n,c(Response_idx, Response_interaction_index)+1] -
                                                  Data_0$X[1:n,c(Response_idx, Response_interaction_index)+1]))
    delta_sparse_random = sparse_coef %*%  as.matrix((Data_1$X[1:n,c(Response_idx, Response_interaction_index)+1] -
                                                         Data_0$X[1:n,c(Response_idx, Response_interaction_index)+1]))

    delta_sparse_scale= apply(delta_sparse_fix+delta_sparse_random, 2, function(x) sqrt(sum(x^2)))
    single_sparse_score = t(res) %*% (delta_sparse_fix+delta_sparse_random)
    single_sparse_score = diag(single_sparse_score) / (delta_sparse_scale)

    #=======================================================#
    #====== step4: generate sparse multiple scores =========#
    #=======================================================#

    multi_sparse_score = c(colSums(res * (delta_sparse_fix)),
                               (res * delta_sparse_random)[union_index,]) ### m by n


    result = list(single_score = single_score, single_sparse_score = single_sparse_score, multi_score = multi_score, multi_sparse_score = multi_sparse_score)

    return(result)
  }

}
