#' Bayesian Binomial Test
#'
#' This function calculates Bayes Factors for every added data point of dichotomous data.
#'
#' The first BF is calculated for 3 data points. For every subsequent data point a new BF is added.
#' The resulting BF vector indicates the change of evidence over time.
#' The function uses "proportionBF" from the BayesFactor package
#'
#' @param data The vector containing binary data.
#' @param p Probability of one result.
#' @param rscale The r of the Prior Cauchy distribution.
#' @examples
#' tbl$bf <- bfbinom(tbl@qbit)
#' @export


# Binomial Seq BF
bfbinom <- function(data, p = .5, rscale = 0.1){
  require(BayesFactor)
  bf <- c(1,1)
  cat("Calculating Sequential Bayes Factors... \n")
  pb = txtProgressBar(min = 3, max = length(data), initial = 3, style = 3)
  for (b in 3:length(data)){
    tmpbfs <- proportionBF(sum(data[1:b]), b, p = p, rscale = rscale)
    bf[b] <- exp(tmpbfs@bayesFactor$bf)
    setTxtProgressBar(pb,b)
  }
  return(bf)
}

# TODO: Bayes t-Tests
