\name{FastMix}
\alias{FastMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   A wrapper for FastMix pipeline.

}
\description{
A wrapper for the main function `ols.eblup.trim()` to conduct deconvolution gene expression analysis with matching cell proportions.
}
\usage{
FastMix(GeneExp, CellProp, Demo, random="all", robust = FASLE, ...)
}
\arguments{
  \item{GeneExp}{
  `GeneExp` is a m by n dimensional gene expression matrix, where m is the number of genes, and n is the number of subjects.
}
  \item{CellProp}{
  `CellProp` is a n by K dimensional matrix of cell proportions, where K is the number of cell proportions.
  }
  \item{Demo}{
  `Clinical` is a n by p dimensional matrix of clinical and demographic variables to be tested, where p is the number of covariates.
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
   \item{trim.fix}{
  Whether only consider trimmed subjects in fix effect estiamtion. The default value is FALSE.
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
\item{GeneExp.fitted }{fitted gene expression.}
\item{sigma.beta}{the covariance estimation of fixed effect.}
\item{VC }{variance component estimation. The first column is the one for common random error. The second column is the one for random effects.}
\item{t.fixed }{the t value for fixed effects. }
\item{eta}{the chi sqiare type statsitics used for p-value calculation.}
\item{p.unadjust}{the overall p-value for outlier detection.}
\item{p.ind.unadjust}{the individual p-value for outlier detection for each random effect.}
\item{out_idx}{the potential covariates with outliers when robust = "FastMix. It is NULL when robust != "FastMix"}
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
## load the data example
## fit the model
mod <- FastMix(GeneExp, CellProp, Demo, random="all", robust = FASLE)

}                               % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~outliers }% __ONLY ONE__ keyword per line
