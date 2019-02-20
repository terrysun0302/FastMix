\name{FastMix}
\alias{FastMix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   A wrapper for FastMix pipeline.

}
\description{
A wrapper for the main function \code{ols.eblup.trim()} to conduct deconvolution gene expression analysis with matching cell proportions.
}
\usage{
FastMix(GeneExp, CellProp, Demo, random="all", include.demo=TRUE, ...)
}
\arguments{
  \item{GeneExp}{
  `GeneExp` is a m by n dimensional gene expression matrix, where m is the number of genes, and n is the number of subjects.
}
  \item{CellProp}{
  `CellProp` is a n by K dimensional matrix of cell proportions, where K is the number of cell proportions.
  }
  \item{random}{
  \code{random} is an index vector that specifies which variable(s)
requires random effects -- by default, \code{random="all"}, which means that all covariates are paired with a random effect.
  }
  \item{Demo}{
  \code{Demo} is a n by p dimensional matrix of clinical and demographic
variables to be tested, where p is the number of covariates. 
  }
  \item{include.demo}{
   Whether the demographical covariates should be
included as the main effects in the model or not. Default to TRUE.
  }
  \item{...}{Additional parameters passed to \code{ols.eblup.trim()}. It
includes the following useful options
 \itemize{
  \item{independent: }{\code{independent}
  specifies the correlation structure among random effects. The default
value is TRUE, which means that all random effects are assumed to be independent.
  }
  \item{trim: }{\code{trim} is
    the trimming percentage when accounting for outliers.
    Default valie is 0.5 (50\%).
  }
   \item{test: }{
    the test method for DEGs. "1" is Gaussian mixture model, "2" is
Anderson-darling normal test. Default value is "1".
  }
  \item{robust: }{\code{robust} specifies whether robust covariance
estimation is implemented and 
which method to use:  "FALSE" for non-robust estimation; "mcd" for the
MCD algorithm of Rousseeuw and Van Driessen; "weighted" for the
Reweighted MCD; "donostah" for the Donoho-Stahel projection based
estimator; "pairwiseQC" for the orthogonalized quadrant correlation
pairwise estimator. All these algorithms come from the R package
`robust`. The default value is \code{robust="FastMix"}, which is the proposed trimming method.
  }
   \item{trim.fix: }{
  Whether only consider trimmed subjects in fix effect estiamtion. The default value is FALSE.
  }}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{fixed.results }{the estimated fix effects and their p-values. They are
overall effects shared by all genes.}
\item{beta.mat }{estimated linear coefficients for individual genes.}
\item{GeneExp.fitted }{fitted gene expressions.}
\item{sigma.beta}{the estimated covariance matrix of the fixed effects.}
\item{VC }{variance component estimation. The first column is the one for common random error. The second column is the one for random effects.}
\item{eta}{the chi-square type statsitics used for p-value calculation.}
\item{re.pvalue}{the overall p-value for detecting outliers in random effects.}
\item{re.ind.pvalue}{the individual p-value for outlier detection for each random effect.}
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
data(dataexample)

## fit the model by default parameters
mod1 <- FastMix(GeneExp, CellProp, Demo)

## some variants: only assign random effects to a subset of covariates, and uses non-robust method to estimate the covariance structure
mod2 <- FastMix(GeneExp, CellProp, Demo, random=c(1,2,10), robust = FALSE)

}                               % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~outliers }% __ONLY ONE__ keyword per line
