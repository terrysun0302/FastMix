\name{DataPrep}
\alias{DataPrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   Preparing the covariates/responses for function \code{ols.eblup.trim()}.
}
\description{
   A convenient function for generating the full covariate matrix and
preparing the response variables in a format that can be directly used
by function \code{ols.eblup.trim()}.
}
\usage{
DataPrep(GeneExp, CellProp, Demo, include.demo=TRUE)
}
\arguments{
  \item{GeneExp}{
  `GeneExp` is a m by n dimensional gene expression matrix, where m is the number of genes, and n is the number of subjects.
}
  \item{CellProp}{
  `CellProp` is a n by K dimensional matrix of cell proportions, where K is the number of cell proportions.
  }
  \item{Demo}{
  \code{Demo} is a n by p dimensional matrix of clinical and demographic
variables to be tested, where p is the number of covariates.
  }
  \item{include.demo}{Whether the demographical covariates should be
    included as the main effects in the model or not. Default to TRUE.}
  \item{w}{The weight matrix.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{X}{A combined covariate matrix with cell proportions as main
effects and all the interaction terms between cell proportions and
demographic covariates.  If \code{include.demo} is \code{TRUE}, the
demographic covariates are also included in the model as main
effects. These covariates are stacked in the row direction so that each
gene has the same covariate matrix (see our manuscript for more
details).  `X` also includes a new variable `ID` which represents each gene.
}
\item{Y}{Vectorized gene expressions to be used as the response variable
in the model.}
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
This function checks the numerical singularity of `X`, the combined
covariate matrix. If the smallest singular value is less than `1e-7`, it
stops and asks the user to reduce the complexity of the model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## load the data example
data(dataexample)

## preparing the covariate/response
des1 <- DataPrep(GeneExp, CellProp, Demo)
dim(des1$X)

## without the demographic main effects
des2 <- DataPrep(GeneExp, CellProp, Demo, include.demo=FALSE)
dim(des2$X)

}                               % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~outliers }% __ONLY ONE__ keyword per line
