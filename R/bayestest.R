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
  if(all.equal(nstart, as.integer(nstart)) != TRUE) stop("nstart must be an integer!")
  if(nstart < 0) stop("nstart must be positive!")
  if(all.equal(inc, as.integer(inc)) != TRUE) stop("inc must be an integer!")
  
  require(BayesFactor)
  bf <- rep(1, (nstart-1))
  cat("N =",length(data),"\n")
  cat("Calculating Sequential Bayes Factors...\n")
  pb = txtProgressBar(min = nstart, max = length(data), initial = nstart, style = 3)
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
  if(0 %in% nullInterval) txt.alternative = "less"
  if(1 %in% nullInterval) txt.alternative = "greater"
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
#' @param x A vector containing data values.
#' @param y A vector containing another set of values for a paired samples test.
#' @param formula A formula of the form var ~ group where var is a numeric variable giving the data values and group is indicating which group the data point belongs to (must be 2 levels). Use this for an independent samples t-Test.
#' @param data Use with formula. A data frame containing the variables given in the formula. Use this for an independent samples t-Test.
#' @param alternative Indicates the direction of the alternative hypothesis: two.sided, less, or greater
#' @param mu A number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param prior.loc Location of the cauchy distributed prior function (use 0 for an uninformed prior).
#' @param prior.r Scale of thr cauchy distributed prior function.
#' @param nstart How many data points should be considered before calculating the first BF (min = 2).
#' @examples
#' bfttest(sumscores, alternative = "greater") # One-sample
#'
#' tbl$bf <- bfttest(tbl$sums, alternative = "two.sided", mu = 50, nstart = 10, prior.loc = 0.1, prior.r = 0.05)$BF
#' 
#' bfttest(formula = raw~group, data = df) # Independent samples
#' 
#' bfttest(df1$scores, df2$scores) # Paired samples
#' @export

bfttest <- function(x=NULL, y = NULL, formula = NULL, data = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, prior.loc = 0, prior.r = 0.1, nstart = 5){
  # calculate t-scores and BFs
  bf <- t <- list()
  
  if(length(alternative) > 1){
    warning("No alternative specified. Using a two-sided alternative.")
    alternative <- alternative[1]
  }
  
  if (!is.null(formula)){ #Independent Samples
    if(!is.null(x)) stop("Please use formula and data for independent and x (and y) for one-sample or paired samples tests.")
    if(is.null(data)) stop("Please specify data.")
    if(length(unique(data[,deparse(formula[[3]])])) != 2) stop("Group must have 2 levels.")
    type <- "independent" 
    samplesize <- nrow(data)
    
    cat("Independent Samples test (N = ",nrow(data),")\nCalculating Sequential Bayes Factors...\n",sep="")
    pb = txtProgressBar(min = nstart, max = nrow(data), style = 3)
    for (i in nstart:nrow(data)) {
      t[[i]] <- t.test(formula=formula, data=data[1:i,], alternative = alternative, paired=F, var.equal=TRUE)
      n1 <- table(data[1:i,deparse(formula[[3]])])[1]
      n2 <- table(data[1:i,deparse(formula[[3]])])[2]
      bf[[i]] <- bf10_t(t = t[[i]][[1]], n1 = n1, n2 = n2, independentSamples = T, prior.location = prior.loc, prior.scale = prior.r, prior.df = 1)
      setTxtProgressBar(pb,i)
    }
  } else {
    if(is.null(x)) stop("Please use formula and data for independent and x (and y) for one-sample or paired samples tests.")
    x <- na.omit(x)
    pb = txtProgressBar(min = nstart, max = length(x), style = 3)
    
    if (!is.null(y)){ # Paired Samples Test
      y <- na.omit(y)
      if (length(y) != length(x)) stop("Data are not the same length!")
      type <- "paired"
      samplesize <- length(x)
      
      cat("Paired Samples test (N = ",length(x),")\nCalculating Sequential Bayes Factors...\n",sep="")
      for (i in nstart:length(x)) {
        t[[i]] <- t.test(x[1:i], y[1:i], alternative = alternative, paired = T, var.equal=TRUE)
        bf[[i]] <- bf10_t(t = t[[i]][[1]], n1 = i, independentSamples = F, prior.location = prior.loc, prior.scale = prior.r, prior.df = 1)
        setTxtProgressBar(pb,i)
      }
      
    } else { # One Sample Test
      if(var(x[1:nstart]) == 0) stop("Cannot compute t-Test since there is no variance in the data. Please choose a larger nstart!")
      type <- "one-sample"
      samplesize <- length(x)
      
      cat("One Sample test (N = ",length(x),")\nCalculating Sequential Bayes Factors...\n",sep="")
      for (i in nstart:length(x)) {
        t[[i]] <- t.test(x[1:i], alternative = alternative, mu = mu)
        bf[[i]] <- bf10_t(t = t[[i]][[1]], n1 = i, prior.location = prior.loc, prior.scale = prior.r, prior.df = 1)
        setTxtProgressBar(pb,i)
      }
    }
  }
  
  close(pb)
  
  tlist <- c(rep(NA,(nstart-1)),unlist(lapply(t, `[[`, 1)))
  plist <- c(rep(NA,(nstart-1)),unlist(lapply(t, `[[`, 3)))
  
  BF10 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 1)))
  BFplus0 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 2)))
  BFmin0 <- c(rep(1,(nstart-1)),unlist(lapply(bf, `[[`, 3)))
  
  if (alternative=="less"){
    bft.out <- list("t-value" = unname(tlist), "p-value" = unname(plist), "BF" = unname(BFmin0), "test type" = type, "prior" = list("Cauchy", "prior location" = prior.loc, "prior scale" = prior.r), "sample size" = samplesize, "alternative" = alternative)
  } else if (alternative=="greater"){
    bft.out <- list("t-value" = unname(tlist), "p-value" = unname(plist), "BF" = unname(BFplus0), "test type" = type, "prior" = list("Cauchy", "prior location" = prior.loc, "prior scale" = prior.r), "sample size" = samplesize, "alternative" = alternative)
  } else {
    bft.out <- list("t-value" = unname(tlist), "p-value" = unname(plist), "BF" = unname(BF10), "test type" = type, "prior" = list("Cauchy", "prior location" = prior.loc, "prior scale" = prior.r), "sample size" = samplesize, "alternative" = alternative)
  }
  cat("Final Bayes Factor: ",tail(bft.out$BF,n=1)," (t=",tail(bft.out$`t-value`,n=1),"; p=",tail(bft.out$`p-value`,n=1),")\n",sep="")
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
