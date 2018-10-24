## generate some new simulated data
n <- 50
set.seed(4567)
CellProp <- cbind(Cell1=runif(n, .2, .3),
                  Cell2=runif(n, .2, .3),
                  Cell3=runif(n, .2, .3))
rownames(CellProp) <- paste0("Subj", 1:n)
Demo <- cbind(Severity=round(runif(n, 0, 10), 1),
              Sex=c(rep(0, n/2), rep(1, n/2)))
rownames(Demo) <- paste0("Subj", 1:n)

## generate some true signals

m <- 5000                              #genes
bb <- c(1.5, 0, 0)                     #assoc. w. CellProp
sigma2.u <- 1
sigma2.e <- 0.8
sigma2.alpha <- 1.5
sigma2.err <- 1/16
#sigma2.err <- 0.5 #i.i.d. noise
## mean interactions
aa <- rep(0, 6); aa[1] <- .75

Bmat <- t(replicate(m, c(bb, 0, 0, aa))); L <- ncol(Bmat)
sdvec <-c(rep(sigma2.u, 3), rep(sigma2.e,2), rep(sigma2.alpha, 6))

Gmat <- matrix(rnorm(L*m),m,L) %*% diag(sdvec)
idx = c(1,4,7,10)
str = sdvec[idx] * 3
OutSig = matrix(NA, ncol = 4, nrow = m)
for(i in 1:4){
  OutSig[,i] <- c(rep(0, m * 0.05 * (i-1)), rep(str[i], m * 0.05), rep(0, m*(1 - i * 0.05)))
}
OutSig <- OutSig * (2*rbinom(length(OutSig), 1, .5) - 1) # approximately 0 mean
Bmat[,idx] = Bmat[,idx] + OutSig

## generate the gene expression
Err <- matrix(rnorm(m*n, sd=sqrt(sigma2.err)), nrow=m, ncol=n)
colnames(Err) <- paste0("Subj", 1:n); rownames(Err) <- paste0("Gene", 1:m)

## Use DataPrep to get Xmat; Y is irrelevant
Xmat <- t(DataPrep(Err, CellProp, Demo)$X[1:n, -1])

## the gene expression matrix
GeneExp <- (Bmat + Gmat) %*% Xmat + Err

save(L,m,n,sigma2.u,sigma2.e, sigma2.alpha, sigma2.err, aa, bb, Demo, CellProp, Xmat, GeneExp, file="dataexample.rda")


########## The actual computation ##########
library(FastMix)
rr1 <- FastMix(GeneExp, CellProp, Demo)
rr2 <- FastMix(GeneExp, CellProp, Demo, robust = "FastMix")
