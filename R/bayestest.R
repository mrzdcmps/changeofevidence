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
#' @param nullInterval Optional vector of length 2 containing lower and upper bounds of an interval hypothesis to test, in probability units.
#' @param nstart How many data points should be considered before calculating the first BF (min = 2)
#' @examples
#' tbl$bf <- bfbinom(tbl@qbit)
#' @export


bfbinom <- function(data, p = 0.5, prior.r = 0.1, nullInterval = NULL, nstart = 5){
  data <- na.omit(data)
  if(length(data) < nstart) stop("Too few observations.")
  
  require(BayesFactor)
  bf <- rep(1, (nstart-1))
  cat("N =",length(data),"\n")
  cat("Calculating Sequential Bayes Factors...\n")
  pb = txtProgressBar(min = 3, max = length(data), initial = 3, style = 3)
  for (b in nstart:length(data)){
    if(!is.null(nullInterval)){
      tmpbfs <- proportionBF(sum(data[1:b]), b, p = p, rscale = prior.r, nullInterval = nullInterval)
    } else{
      tmpbfs <- proportionBF(sum(data[1:b]), b, p = p, rscale = prior.r)
    }
    bf[b] <- exp(tmpbfs@bayesFactor$bf)[1]
    setTxtProgressBar(pb,b)
  }
  close(pb)
  
  if(is.null(nullInterval)) txt.alternative = "two.sided"
  if(0 %in% nullInterval) txt.alternative = "greater"
  if(1 %in% nullInterval) txt.alternative = "less"
  orthodoxtest <- binom.test(sum(data),length(data),p = p, alternative = txt.alternative)
  
  cat("Final Bayes Factor: ",tail(bf,n=1)," (probability of success=",orthodoxtest$estimate,"; p=",orthodoxtest$p.value,")",sep="")
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
#' @param ydata A vector containing another set of values for a paired samples test (independent samples are currently not supported). Please indicate with paired = TRUE.
#' @param alternative Indicates the direction of the alternative hypothesis: two.sided, less, or greater
#' @param mu A number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param paired Logical. Set TRUE if you want to perform a paired samples t-Test.
#' @param prior.loc Location of the cauchy distributed prior function (use 0 for an uninformed prior).
#' @param prior.r Scale of thr cauchy distributed prior function.
#' @param nstart How many data points should be considered before calculating the first BF (min = 2).
#' @examples
#' bflist <- bfttest(sumscores, alternative = "greater")
#'
#' tbl$bf <- bfttest(tbl$sums, alternative = "two.sided", mu = 50, nstart = 10, prior.loc = 0.1, prior.r = 0.05)$BF
#' @export

bfttest <- function(data, ydata = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, prior.loc = 0, prior.r = 0.1, nstart = 5){
  # calculate t-scores and BFs
  bf <- t <- list()
  data <- na.omit(data)
  cat("N =",length(data),"\n")
  cat("Calculating Sequential Bayes Factors...\n")
  pb = txtProgressBar(min = nstart, max = length(data), style = 3)
  if (!is.null(ydata)){
    ydata <- na.omit(ydata)
    if (length(ydata) != length(data)) stop("Data are not the same length!")
    if (paired == TRUE) { # Paired Samples Test
      for (i in nstart:length(data)) {
        t[[i]] <- t.test(data[1:i], ydata[1:i], alternative = alternative, paired = T)
        bf[[i]] <- bf10_t(t = t[[i]][[1]], n1 = i, independentSamples = F, prior.location = prior.loc, prior.scale = prior.r, prior.df = 1)
        setTxtProgressBar(pb,i)
      }
    } else { # Independent Samples Test
        stop("Independent Samples t-Test currently not supported!")
    }
  } else { # One Sample Test
    if(var(data[1:nstart]) == 0) stop("Cannot compute t-Test since there is no variance in the data. Please choose a larger nstart!")
      for (i in nstart:length(data)) {
        t[[i]] <- t.test(data[1:i], alternative = alternative, mu = mu)
        bf[[i]] <- bf10_t(t = t[[i]][[1]], n1 = i, prior.location = prior.loc, prior.scale = prior.r, prior.df = 1)
        setTxtProgressBar(pb,i)
      }
  }
  
  close(pb)
  
  tlist <- c(rep(NA,(nstart-1)),unlist(lapply(t, `[[`, 1)))
  plist <- c(rep(NA,(nstart-1)),unlist(lapply(t, `[[`, 3)))
  
  BF10 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 1)))
  BFplus0 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 2)))
  BFmin0 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 3)))
  
  if (alternative=="less"){
    bft.out <- list("t-value" = unname(tlist), "p-value" = unname(plist), "BF" = unname(BFmin0))
  } else if (alternative=="greater"){
    bft.out <- list("t-value" = unname(tlist), "p-value" = unname(plist), "BF" = unname(BFplus0))
  } else {
    bft.out <- list("t-value" = unname(tlist), "p-value" = unname(plist), "BF" = unname(BF10))
  }
  cat("Final Bayes Factor: ",tail(bft.out$BF,n=1)," (t=",tail(bft.out$`t-value`,n=1),"; p=",tail(bft.out$`p-value`,n=1),")",sep="")
  return(bft.out)
}


#' Bayesian Sequential Correlation Test
#'
#' This function calculates Bayes Factors for the correlation of datasets.
#'
#' The first BF is calculated for nstart data points. For every subsequent data point a new BF is added.
#' The resulting BF vector indicates the change of evidence over time.
#' The function uses "correlationBF" from the BayesFactor package
#'
#' @param x A vector containing continous data.
#' @param y A second vector containing continous data.
#' @param alternative specify the direction of the alternative hypothesis: "greater", "less", "two.sided"
#' @param prior.r Prior distribution (scaled beta)
#' @param nstart How many data points should be considered before calculating the first BF (min = 2)
#' @examples
#' bfcor(exp$sums, con$sums, nullInterval = c(-1,0))
#' @export


# Binomial Seq BF
bfcor <- function(x, y, alternative = "two.sided", prior.r = 0.1, nstart = 5){
  
  x <- na.omit(x)
  y <- na.omit(y)
  
  if(alternative == "greater") nullInterval <- c(-1,0)
  else if(alternative == "less") nullInterval <- c(0,1)
  else nullInterval <- 0
  
  require(BayesFactor)
  if (length(y) != length(x)) 
    stop("Length of y and x must be the same.")
  
  bf <- rep(1, (nstart-1))
  cat("N =",length(data),"\n")
  cat("Calculating Sequential Bayes Factors...\n")
  pb = txtProgressBar(min = nstart, max = length(x), initial = nstart, style = 3)
  for (b in nstart:length(x)){
    bf[b] <- exp(correlationBF(x[1:b], y[1:b], rscale=prior.r, nullInterval = nullInterval)@bayesFactor$bf[2])
    setTxtProgressBar(pb,b)
  }
  close(pb)
  
  orthodoxtest <- cor.test(x, y, use = "complete.obs", alternative = alternative)
  
  cat("Final Bayes Factor: ",tail(bf,n=1)," (r=",orthodoxtest$estimate,"; p=",orthodoxtest$p.value,")",sep="")
  return(bf)
}
