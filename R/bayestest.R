#' Bayesian Sequential Binomial Test
#'
#' This function calculates Bayes Factors for every added data point of dichotomous data.
#'
#' The first BF is calculated for nstart data points. For every subsequent data point a new BF is added.
#' The resulting BF vector indicates the change of evidence over time.
#' The function uses "proportionBF" from the BayesFactor package
#'
#' @param data A vector containing binary data.
#' @param p Probability of one result.
#' @param prior.r The r of the Prior Cauchy distribution.
#' @param nstart How many data points should be considered before calculating the first BF (min = 2)
#' @examples
#' tbl$bf <- bfbinom(tbl@qbit)
#' @export


# Binomial Seq BF
bfbinom <- function(data, p = 0.5, prior.r = 0.1, nstart = 3){
  require(BayesFactor)
  bf <- rep(1, (nstart-1))
  cat("Calculating Sequential Bayes Factors... \n")
  pb = txtProgressBar(min = 3, max = length(data), initial = 3, style = 3)
  for (b in 3:length(data)){
    tmpbfs <- proportionBF(sum(data[1:b]), b, p = p, rscale = prior.r)
    bf[b] <- exp(tmpbfs@bayesFactor$bf)
    setTxtProgressBar(pb,b)
  }
  return(bf)
}


#' Bayesian Sequential t-Test
#'
#' This function calculates Bayes Factors for every added data point of normally distributed data, e.g. sum scores.
#'
#' The first BF is calculated for Nstart data points. For every subsequent data point a new BF is added.
#' The resulting BF vector indicates the change of evidence over time.
#' The function first calculates t-scores and p-values and subsequently uses the BFDA package to translate the data into BFs.
#' Returns a list of t-scores ($t-value), p-scores ($p-value), and BF ($BF)
#'
#' @param data A vector containing data values.
#' @param alternative Indicates the direction of the alternative hypothesis: two.sided, less, or greater
#' @param mu A number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param prior.loc Location of the cauchy distributed prior function (use 0 for an uninformed prior).
#' @param prior.r Scale of thr cauchy distributed prior function.
#' @param nstart How many data points should be considered before calculating the first BF (min = 2).
#' @examples
#' bflist <- bfttest(sumscores, alternative = "greater")
#'
#' tbl$bf <- bfttest(tbl$sums, alternative = "two.sided", mu = 50, nstart = 10, prior.loc = 0.1, prior.r = 0.05)$BF
#' @export

bfttest <- function(data, alternative = c("two.sided", "less", "greater"), mu = 0, prior.loc = 0, prior.r = 0.1, nstart = 3){
  # calculate t-scores and BFs
  bf <- t <- list()
  cat("Calculating Sequential Bayes Factors... \n")
  pb = txtProgressBar(min = nstart, max = length(data), initial = 3, style = 3)
  for (i in nstart:length(data)) {
    t[[i]] <- t.test(data[1:i], alternative = alternative, mu = mu)
    bf[[i]] <- bf10_t(t = t[[i]][[1]], n1 = i, prior.location = prior.loc, prior.scale = prior.r, prior.df = 1)
    setTxtProgressBar(pb,i)
  }
  
  tlist <- c(rep(NA,(nstart-1)),unlist(lapply(t, `[[`, 1)))
  plist <- c(rep(NA,(nstart-1)),unlist(lapply(t, `[[`, 3)))
  
  BF10 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 1)))
  BFplus0 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 2)))
  BFmin0 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 3)))
  
  if (alternative=="less"){
    bft.out <- list("t-value" = tlist, "p-value" = plist, "BF" = BFmin0)
  } else if (alternative=="greater"){
    bft.out <- list("t-value" = tlist, "p-value" = plist, "BF" = BFplus0)
  } else {
    bft.out <- list("t-value" = tlist, "p-value" = plist, "BF" = BF10)
  }
  return(bft.out)
}