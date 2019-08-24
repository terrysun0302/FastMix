\name{dat_train}
\alias{dat_train}
\docType{data}
\title{
dat_train
}
\description{
A simulated data exmaple to show the data structure needed for the deconvolution application.
}
\usage{data(dat_train)}
\format{
  \describe{
    \item{\code{CellProp}}{`CellProp` is an n by K dimensional matrix of cell proportions, where K is the number of cell proportions.}
    \item{\code{GeneExp}}{ `GeneExp` is an m by n dimensional gene expression matrix, where m is the number of genes, and n is the number of subjects.}
    \item{\code{Demo}}{ `Clinical` is a n by p dimensional matrix of clinical and demographic variables to be tested, where p is the number of covariates.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(dat_train)
}
\keyword{datasets}
