\name{ols.eblup}
\alias{ols.eblup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit Linear Mixed Models by an OLS approach
}
\description{
Fit a Linear Mixed Model (LMM) to data, via an OLS-based approach.
}
\usage{
ols.eblup(Des, Y, random, independent = T, type = "F")
}
\arguments{
  \item{Des}{
the design matrix ordered by time subject by subject. First column should be ID range from 1 to the sample size;
the rest columns are covariates. See the data example for more details.
}
  \item{Y}{
  the longitudinal outcome ordered by time subject by subject.
  }
  \item{random}{
  the user selected covariates with random effects.
  }
  \item{independent}{
  if the random effects are corrected or not. T or F
  }
  \item{type}{
  the method to calculate the p-values. F or chi.
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
\item{Yhat }{fitted response.}
\item{sigma.beta}{the covariance estimation of fixed effect.}
\item{VC }{variance component estimation. The first column is the one for common random error. The second column is the one for random effects.}
\item{t.fixed }{the t value for fixed effects. }
\item{eta}{the chi sqiare type statsitics used for p-value calculation.}
\item{p.unadjust}{the overall p-value for outlier detection.}
\item{p.ind.unadjust}{the individual p-value for outlier detection for each random effect.}
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
Des <- dataexample$Des
Y <- dataexample$Y

## fit the model
mod <- ols.eblup(Des, Y, random = c(1,2,3,4), independent = T,type = "chi")

}                               % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~random effects }% __ONLY ONE__ keyword per line
