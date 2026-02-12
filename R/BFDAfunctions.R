# ==============================================================================
# These are the functions for t-tests with informed priors
# ==============================================================================
# see also https://arxiv.org/abs/1704.02479 for the formulae

#@import hypergeo

# helper functions for the computation of the Bayes factor with informed priors

A <- function(t, n, nu, mu.delta, g) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = mu.delta^2*t^2/
                             (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
}

B <- function(t, n, nu, mu.delta, g) {
  
  out <- mu.delta*t/sqrt(1/2*(1/n + g)*((1 + n*g)*nu + t^2)) *
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2)) *
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = mu.delta^2*t^2/
                               (2*(1/n + g)*((1 + n*g)*nu + t^2))))
  
  return(out)
  
}


C <- function(delta, t, n, nu) {
  
  Re(hypergeo::genhypergeo(U = (nu + 1)/2, L = 1/2,
                           z = n*t^2*delta^2/(2*(nu + t^2))))
  
}

D <- function(delta, t, n, nu) {
  
  out <- t*delta*sqrt(2*n/(nu + t^2))*
    exp(lgamma((nu + 2)/2) - lgamma((nu + 1)/2))*
    Re(hypergeo::genhypergeo(U = (nu + 2)/2, L = 3/2,
                             z = n*t^2*delta^2/(2*(nu + t^2))))
  
  return(out)
  
}

term_normalprior <- function(t, n, nu, mu.delta, g) {
  
  (1 + n*g)^(-1/2) * exp(-mu.delta^2/(2*(1/n + g))) *
    (1 + t^2/(nu*(1 + n*g)))^(-(nu + 1)/2) *
    (A(t, n, nu, mu.delta, g) + B(t, n, nu, mu.delta, g))
  
}

integrand <- function(g, t, n, nu, mu.delta, r, kappa) {
  
  tmp <- term_normalprior(t = t, n = n, nu = nu, mu.delta = mu.delta, g = g)
  pg_log <- kappa/2*(2*log(r) + log(kappa/2)) - lgamma(kappa/2) -
    (kappa/2 + 1)*log(g) - r^2*kappa/(2*g)
  pg <- exp(pg_log)
  out <- tmp*pg
  
  return(out)
  
}

dtss <- function(delta, mu.delta, r, kappa, log = FALSE) {
  
  out <- - log(r) + lgamma((kappa + 1)/2) - .5*(log(pi) + log(kappa)) -
    lgamma(kappa/2) - (kappa + 1)/2 * log(1 + ((delta - mu.delta)/r)^2/kappa)
  
  if ( ! log)
    out <- exp(out)
  
  return(out)
  
}

posterior_t_tmp <- function(delta, t, n1, n2 = NULL, independentSamples = FALSE,
                            prior.location, prior.scale, prior.df,
                            rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, n1*n2/(n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.location
  r <- prior.scale
  kappa <- prior.df
  
  numerator <- exp(-neff/2*delta^2)*(1 + t^2/nu)^(-(nu + 1)/2)*
    (C(delta, t, neff, nu) + D(delta, t, neff, nu))*
    dtss(delta, mu.delta, r, kappa)
  
  denominator <- integrate(integrand, lower = 0, upper = Inf,
                           t = t, n = neff, nu = nu, mu.delta = mu.delta,
                           r = r, kappa = kappa, rel.tol = rel.tol)$value
  
  out <- numerator/denominator
  
  if ( is.na(out))
    out <- 0
  
  return(out)
  
}

posterior_t <- Vectorize(posterior_t_tmp, "delta")

cdf_t <- function(x, t, n1, n2 = NULL, independentSamples = FALSE,
                  prior.location, prior.scale, prior.df) {
  
  area <- integrate(posterior_t, lower = -Inf, upper = x, t = t, n1 = n1, n2 = n2,
                    independentSamples = independentSamples,
                    prior.location = prior.location, prior.scale = prior.scale,
                    prior.df = prior.df)$value
  
  if(area > 1)
    area <- 1
  
  return(area)
  
}


# Function to compute the Bayes factor with t distribution as prior

bf10_t <- function(t, n1, n2 = NULL, independentSamples = FALSE, prior.location,
                   prior.scale, prior.df, rel.tol = .Machine$double.eps^0.25) {
  
  neff <- ifelse(independentSamples, n1*n2/(n1 + n2), n1)
  nu <- ifelse(independentSamples, n1 + n2 - 2, n1 - 1)
  
  mu.delta <- prior.location
  r <- prior.scale
  kappa <- prior.df
  numerator <- integrate(integrand, lower = 0, upper = Inf,
                         t = t, n = neff, nu = nu, mu.delta = mu.delta,
                         r = r, kappa = kappa,
                         rel.tol = rel.tol)$value
  denominator <- (1 + t^2/nu)^(-(nu + 1)/2)
  
  BF10 <- numerator/denominator
  priorAreaSmaller0 <- integrate(dtss, lower = -Inf, upper = 0,
                                 mu.delta = prior.location, r = prior.scale,
                                 kappa = prior.df)$value
  postAreaSmaller0 <- cdf_t(x = 0, t = t, n1 = n1, n2 = n2,
                            independentSamples = independentSamples,
                            prior.location = prior.location,
                            prior.scale = prior.scale, prior.df = prior.df)
  BFmin1 <- postAreaSmaller0/priorAreaSmaller0
  BFplus1 <- (1 - postAreaSmaller0)/(1 - priorAreaSmaller0)
  BFmin0 <- BFmin1 * BF10
  BFplus0 <- BFplus1 * BF10
  
  return(list(BF10 = BF10, BFplus0 = BFplus0, BFmin0 = BFmin0))
  
}


# ==============================================================================
# These are functions for nonparametric t-tests
# ==============================================================================
# van Doorn: Bayesian Rank-Based Hypothesis Testing for the Rank Sum Test, the Signed Rank Test, and Spearman's rho
# see https://osf.io/gny35

rankSumGibbsSampler <- function(xVals, yVals, nSamples = 1e3, cauchyPriorParameter = 1/sqrt(2),  progBar = TRUE, 
                                nBurnin = 1, nGibbsIterations = 10, nChains = 10) {
  
  if (progBar) {
    myBar <- txtProgressBar(min = 1, max = nSamples*nChains, initial = 1, char = "*",style=3,width=50)
  }
  
  n1 <- length(xVals)
  n2 <- length(yVals)
  
  allRanks <- rank(c(xVals,yVals))
  xRanks <- allRanks[1:n1]
  yRanks <- allRanks[(n1+1):(n1+n2)]
  
  deltaSamples <- numeric(nSamples)
  deltaSamplesMatrix <- matrix(ncol = nChains, nrow = nSamples-nBurnin)
  totalIterCount <- 0
  
  for(thisChain in 1:nChains) {
    
    currentVals <- sort(rnorm((n1+n2)))[allRanks] # initial values
    
    oldDeltaProp <- 0
    
    for (j in 1:nSamples) {
      
      for (i in sample(1:(n1+n2))) {
        
        currentRank <- allRanks[i]
        
        currentBounds <- upperLowerTruncation(ranks=allRanks, values=currentVals, currentRank=currentRank)
        if (i <= n1) {
          oldDeltaProp <- -0.5*oldDeltaProp
        } else {
          oldDeltaProp <- 0.5*oldDeltaProp
        }
        
        currentVals[i] <- truncNormSample(currentBounds[["under"]], currentBounds[["upper"]], mu=oldDeltaProp, sd=1)
        
      }
      
      xVals <- currentVals[1:n1]
      yVals <- currentVals[(n1+1):(n1+n2)]
      
      gibbsResult <- sampleGibbsTwoSampleWilcoxon(x = xVals, y = yVals, nIter = nGibbsIterations,
                                                  rscale = cauchyPriorParameter)
      
      deltaSamples[j] <- oldDeltaProp <- gibbsResult
      if (progBar) setTxtProgressBar(myBar, j + ( (thisChain-1) * nSamples)) 
      
    }
    
    if (nBurnin > 0) {
      deltaSamples <- -deltaSamples[-(1:nBurnin)]
    } else {
      deltaSamples <- -deltaSamples
    }
    
    deltaSamplesMatrix[, thisChain] <- deltaSamples
    
  }
  
  betweenChainVar <- (nSamples / (nChains - 1)) * sum((apply(deltaSamplesMatrix, 2, mean)  - mean(deltaSamplesMatrix))^2)
  withinChainVar <- (1/ nChains) * sum(apply(deltaSamplesMatrix, 2, var))
  
  fullVar <- ((nSamples - 1) / nSamples) * withinChainVar + (betweenChainVar / nSamples)
  rHat <- sqrt(fullVar/withinChainVar)
  
  return(list(deltaSamples = as.vector(deltaSamplesMatrix), rHat = rHat))
}

sampleGibbsTwoSampleWilcoxon <- function(x, y, nIter = 10, rscale = 1/sqrt(2)) {
  meanx <- mean(x)
  meany <- mean(y)
  n1 <- length(x)
  n2 <- length(y)
  sigmaSq <- 1 # Arbitrary number for sigma
  g <- 1
  for(i in 1:nIter){
    #sample mu
    varMu <- (4 * g * sigmaSq) / ( 4 + g * (n1 + n2) )
    meanMu <- (2 * g * (n2 * meany - n1 * meanx)) / ((g * (n1 + n2) + 4))
    mu <- rnorm(1, meanMu, sqrt(varMu))
    # sample g
    betaG <- (mu^2 + sigmaSq * rscale^2) / (2*sigmaSq)
    g <- 1/rgamma(1, 1, betaG)
    # convert to delta
    delta <- mu / sqrt(sigmaSq)
  }
  return(delta)
}


truncNormSample <- function(lBound = -Inf, uBound = Inf, mu = 0, sd = 1) {
  
  lBoundUni <- pnorm(lBound, mean = mu, sd = sd)
  uBoundUni <- pnorm(uBound, mean = mu, sd = sd)
  mySample <- qnorm(runif(1, lBoundUni, uBoundUni), mean = mu, sd = sd)
  
  return(mySample)
}

upperLowerTruncation <- function(ranks, values, currentRank) {
  
  if (currentRank == min(ranks)) {
    under <- -Inf
  } else {
    under <- max(values[ranks < currentRank])
  }
  
  if (currentRank == max(ranks)) {
    upper <- Inf
  } else {
    upper <- min(values[ranks > currentRank])
  }
  
  return(list(under=under, upper=upper))
}

signRankGibbsSampler <- function(xVals, yVals = NULL, nSamples = 1e3, cauchyPriorParameter = 1/sqrt(2), testValue = 0, 
                                 progBar = TRUE, nBurnin = 1, nGibbsIterations = 10, nChains = 10) {
  
  if (progBar) {
    myBar <- txtProgressBar(min = 1, max = nSamples*nChains, initial = 1, char = "*",style=3,width=50)
  }
  
  n <- length(xVals)
  
  if (!is.null(yVals)) { 
    differenceScores <- xVals - yVals
  } else {
    differenceScores <- xVals - testValue
  }
  
  differenceSigns <- (sign(differenceScores))
  absDifferenceRanked <- rank(abs(differenceScores))
  prodSignAbsRank <- differenceSigns * absDifferenceRanked
  
  diffSamples <- numeric(n)
  
  
  deltaSamples <- numeric(nSamples)
  deltaSamplesMatrix <- matrix(ncol = nChains, nrow = nSamples-nBurnin)
  oldDeltaProp <- 0
  
  for(thisChain in 1:nChains) {
    
    initDiffSamples <- sort(abs(rnorm(n)))[absDifferenceRanked]
    sampledDiffsAbs <- abs(initDiffSamples)
    
    for (j in 1:nSamples) {
      
      for (i in sample(1:n)) {
        
        currentRank <- absDifferenceRanked[i]
        
        currentBounds <- upperLowerTruncation(ranks=absDifferenceRanked, values=sampledDiffsAbs, currentRank=currentRank)
        if (is.infinite(currentBounds[["under"]])) {currentBounds[["under"]] <- 0}
        
        sampledDiffsAbs[i] <- truncNormSample(currentBounds[["under"]], currentBounds[["upper"]], mu = abs(oldDeltaProp), sd=1)
      }
      diffSamples <- sampledDiffsAbs * differenceSigns
      
      if (any(differenceSigns == 0)) {
        nullSamples <- sampledDiffsAbs[differenceSigns == 0] * sample(c(-1,1), size = sum(differenceSigns == 0), replace = TRUE)
        diffSamples[which(differenceSigns == 0)] <- nullSamples
      }
      
      gibbsOutput <- sampleGibbsOneSampleWilcoxon(diffScores = diffSamples, nIter = nGibbsIterations, rscale = cauchyPriorParameter)
      
      deltaSamples[j] <- oldDeltaProp <- gibbsOutput
      if (progBar) setTxtProgressBar(myBar, j + ( (thisChain-1) * nSamples))
      
    }
    
    if (nBurnin > 0) {
      deltaSamples <- deltaSamples[-(1:nBurnin)]
    } else {
      deltaSamples <- deltaSamples
    }
    deltaSamplesMatrix[, thisChain] <- deltaSamples
  }
  
  betweenChainVar <- (nSamples / (nChains - 1)) * sum((apply(deltaSamplesMatrix, 2, mean)  - mean(deltaSamplesMatrix))^2)
  withinChainVar <- (1/ nChains) * sum(apply(deltaSamplesMatrix, 2, var))
  
  fullVar <- ((nSamples - 1) / nSamples) * withinChainVar + (betweenChainVar / nSamples)
  rHat <- sqrt(fullVar/withinChainVar)
  
  return(list(deltaSamples = as.vector(deltaSamplesMatrix), rHat = rHat))
}

sampleGibbsOneSampleWilcoxon <- function(diffScores, nIter = 10, rscale = 1/sqrt(2)){
  ybar <- mean(diffScores)
  n <- length(diffScores)
  sigmaSq <- 1
  mu <- ybar
  g <- ybar^2 / sigmaSq + 1
  
  for(i in 1:nIter){   
    #sample mu
    varMu  <- sigmaSq / (n + (1 / g))
    meanMu <- (n * ybar) / (n + (1 / g))
    mu <- rnorm(1, meanMu, sqrt(varMu) )
    
    # sample g
    scaleg <- (mu^2 + sigmaSq * rscale^2) / (2*sigmaSq)
    g = 1 / rgamma(1, 1, scaleg )
    
    delta <- mu / sqrt(sigmaSq)
  }
  return(delta)
}


# this function computes  BF10 for both Wilcoxon tests and Spearman's rho, as specified in whichTest. Recommended values
# for Wilcoxon are 1/sqrt(2) and 1 for Spearman. These should be the same as specified in the Gibbs sampler function call.
# The oneSided argument can be FALSE (for two-sided tests), "right" for positive one-sided tests, and "left" for negative
# one-sided tests.
computeBayesFactorOneZero <- function(posteriorSamples, priorParameter = 1, oneSided = FALSE, whichTest = "Wilcoxon") {
  
  postDens <- logspline::logspline(posteriorSamples)
  densZeroPoint <- logspline::dlogspline(0, postDens)
  
  corFactorPosterior <- logspline::plogspline(0, postDens)
  if (oneSided == "right")
    corFactorPosterior <- 1 - corFactorPosterior
  
  if (whichTest == "Wilcoxon") {
    # priorParameter should be the Cauchy scale parameter
    priorDensZeroPoint <- dcauchy(0, scale = priorParameter)
    corFactorPrior <-  pcauchy(0, scale = priorParameter, lower.tail = (oneSided != "right" ))
  } else if (whichTest == "Spearman") {
    # priorParameter should be kappa
    priorDensZeroPoint <- dbeta(0.5, 1/priorParameter, 1/priorParameter) / 2
    corFactorPrior <-  pbeta(0.5, 1/priorParameter, 1/priorParameter, lower.tail = (oneSided != "right" ))
  }
  
  if (isFALSE(oneSided)) {
    bf10 <- priorDensZeroPoint / densZeroPoint
  } else {
    bf10 <- (priorDensZeroPoint / corFactorPrior) / (densZeroPoint / corFactorPosterior)
  }
  
  return(bf10)
}
