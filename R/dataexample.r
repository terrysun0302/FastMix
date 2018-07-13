
#' A simulated data example to illustrate the data structure needed in our functions.
#'
#' A simulated data example. The sample size is 5000 and the number of time points is 50. In Des, the first column is ID,
#' the second column is the intercept and the rest columns are covariates. Y is the observed response.
#'
#' @docType data
#'
#' @usage dataexample
#'
#' @format A data frame with 250000 rows and 4 variables:
#' \describe{
#'   \item{Des}{Design matrix}
#'   \item{Y}{response}
#'   ...
#' }
#'
#'@examples
#' Des <- dataexample$Des
#' Y <- dataexample$Y


"dataexample"
