\name{DataPrep_test}
\alias{DataPrep_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Preparing the covariates/responses for function \code{score_func()}.
}
\description{
  In reality, we do not know the true response status of the testing
  data. However, we need to create a "Response" for these subjects
  in order to do a proper standardization for the covariates.  These
  two sets of data are standardized according to two different
  "hypothetical" responses: Response==0 (Data_0), Response==1
  (Data_1).
}
\usage{
  DataPrep_test(GeneExp, CellProp, Demo, include.demo=TRUE, train_response, w)
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
  \item{train_response}{The column of reponse variable in the training data.}
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
  \item{Data_0}{The data set with pseudo respone Res = 0.
  }
  \item{Data_1}{The data set with pseudo respone Res = 1.
  }
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
  data(dat_test)
  data(dat_train)

  ## preparing the covariate/response
  test_data = DataPrep_test(dat_test$GeneExp, dat_test$CellProp, Demo=dat_test$Demo[,1], train_response = dat_train$Demo[,2], include.demo=TRUE, w)


  Data_0 = test_data$Data_0
  Data_1 = test_data$Data_1

}                               % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~outliers }% __ONLY ONE__ keyword per line