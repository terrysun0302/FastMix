\name{ols.eblup.trim}
\alias{ols.eblup.trim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The main function to conduct FastMix pipeline.

}
\description{
 A new analytic pipeline, dubbed as FastMix, that combines the deconvolution step with the downstream analyses based on linear
mixed eﬀects regression (LMER) model
}
\usage{
ols.eblup.trim (Des, Y, random = "all", independent = T, trim = 0.5, robust = FALSE)
}
\arguments{
  \item{Des}{
  the design matrix ordered by gene subject by subject. First column should be identification variable, e.g., ID or subject, and the rest columns are covariates.
}
  \item{Y}{
  vectorized gene expression data.
  }
  \item{random}{
  `random` is an index vector that specifies which variable(s) requires random effects -- by default, all covariates are paired with a random effect.
  }
  \item{independent}{
  specify the correlation structure among random effects. If TRUE, random effects are assumed to be independent.
  }
  \item{trim}{
    the trimming percentage when accounting for outliers.
    Default valie is 0.5 (50\%).
  }
  \item{robust}{
  Specifies whether robust covariance estimation is implemented and which method to use:  "FALSE" for non-robust estimation; "mcd" for the MCD algorithm of    Rousseeuw and Van Driessen; "weighted" for the Reweighted MCD; "donostah" for the Donoho-Stahel projection based estimator; "pairwiseQC" for the     orthogonalized quadrant correlation pairwise estimator. All these algorithms come from the R package `robust`. "FastMix" is the proposed trimming method.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{beta.hat }{the fix effect estimation.}
\item{beta.mat }{individual coefficient estimation.}
\item{Yhat}{fitted response.}
\item{sigma.beta}{the covariance estimation of fixed effect.}
\item{VC }{variance component estimation. The first column is the one for common random error. The second column is the one for random effects.}
\item{t.fixed }{the t value for fixed effects. }
\item{eta}{the chi sqiare type statsitics used for p-value calculation.}
\item{p.unadjust}{the overall p-value for outlier detection.}
\item{p.ind.unadjust}{the individual p-value for outlier detection for each random effect.}
\item{out_idx}{he potential covariates with outliers when robust = "FastMix. It is NULL when robust != "FastMix"}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hao Sun
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## load the data example and transform the data
gnames <- rownames(GeneExp); m <- nrow(GeneExp)
if (is.null(gnames)) {
  rownames(GeneExp) <- gnames <- paste0("Gene", 1:m)
}
Data2 <- DataPrep(GeneExp, CellProp, Demo)

## fit the model
mod <- ols.eblup.trim(Des=Data2$X, Y=Data2$Y, random="all", robust = "FastMix")

}                               % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~outliers }% __ONLY ONE__ keyword per line
