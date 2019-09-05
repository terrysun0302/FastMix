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
  score_func(mod1, Data_0, Data_1, Response_interaction_index,
Response_idx, sig.level=0.05)
}
\arguments{
  \item{mod1}{
    An object from FastMix model using training data.
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
  \item{sig.level}{
    By design, the sparse scores are linear combinations of genes that
    have significant interactions with the response. `sig.level` is the
    significance level that defines those significant genes. Default to 0.05.
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
  \item{multi_sparse_score}{Multiple summary scores using only genes
    with significant interactions with the response variable.
  }
  %% ...
}
%% ...

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
                       % end examples.

\keyword{models}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~outliers }% __ONLY ONE__ keyword per line
