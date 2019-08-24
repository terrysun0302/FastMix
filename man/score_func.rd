\name{score_func}
\alias{score_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The main function for the score calculation.
}
\description{
  A function for discriminant analysis after the FaxtMix analysis.
}
\usage{
  score_func = function(mod1, Data_0, Data_1, Response_interaction_index = c(5,6,7), Response_idx = c(4))
}
\arguments{
  \item{mod1}{
    An object from FatxMix model using training data.
  }
  \item{Data_0}{
    See the function \code{DataPrep_test()}.
  }
  \item{Data_1}{
    See the function \code{DataPrep_test()}.
  }
  \item{Response_interaction_index}{
    The indexs of the interaction term between resposne and other variables in Data_0 and Data_1.
  }
  \item{Response_idx}{
    The indexs of the resposne response variables in Data_0 and Data_1.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  \item{single_score}{A single summary score using all genes.
  }
  \item{single_sparse_score}{A single summary score using selected genes.
  }
  \item{multi_score}{Multiple summary scores using all genes.
  }
  \item{multi_sparse_score}{Multiple summary scores using selected genes.
  }
  %% ...
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
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## load the data example and transform the data
data(dat_train)
data(dat_test)
n2 <- 200

w = diag(rep(1, n2))

gnames <- rownames(dat_test$GeneExp); m <- nrow(dat_test$GeneExp)
if (is.null(gnames)) {
  rownames(dat_test$GeneExp) <- gnames <- paste0("Gene", 1:m)
}


### this code use the data in main simulation. The difference is that Demo includes other variable than response
test_data = DataPrep_test(dat_test$GeneExp, dat_test$CellProp, Demo=dat_test$Demo[,1], train_response = dat_train$Demo[,2], include.demo=TRUE, w)


Data_0 = test_data$Data_0
Data_1 = test_data$Data_1

predicted_score = score_func(mod1, Data_0, Data_1, Response_interaction_index = c(9,10,11), Response_idx = c(5))
}                               % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~outliers }% __ONLY ONE__ keyword per line
