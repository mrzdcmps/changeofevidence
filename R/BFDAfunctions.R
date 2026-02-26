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


#' Extract Cauchy Prior from Posterior of a Bayesian t-Test
#'
#' Fits a Cauchy distribution to the posterior distribution of a parametric
#' Bayesian t-test. The fitted Cauchy location and scale can be used as informed
#' prior parameters for an exact replication study.
#'
#' @param x A \code{seqbf} object created with \code{bfttest()} (parametric only)
#' @param n.points Number of grid points for fitting (default: 1000)
#' @param delta.range Optional numeric vector of length 2 specifying the range
#'   of delta values for the fitting grid. If NULL (default), automatically
#'   determined from the posterior mode and width.
#' @return An object of class \code{"cauchyFit"} containing:
#'   \itemize{
#'     \item \code{location}: Fitted Cauchy location parameter
#'     \item \code{scale}: Fitted Cauchy scale parameter
#'     \item \code{r_squared}: R-squared goodness-of-fit measure
#'     \item \code{max_abs_deviation}: Maximum absolute deviation between
#'       posterior and fitted Cauchy
#'     \item \code{prior}: List with original prior location and scale
#'     \item \code{test_type}: Type of t-test (one-sample, paired, independent)
#'     \item \code{sample_size}: Sample size(s) from the original test
#'   }
#' @details
#' The function extracts the final t-statistic and sample sizes from the
#' \code{seqbf} object, evaluates the analytical posterior distribution using
#' \code{posterior_t()}, and fits a Cauchy distribution by minimizing the
#' integrated squared error via Nelder-Mead optimization. The posterior mode
#' and half-width at half-maximum are used as starting values.
#'
#' Only parametric t-tests are supported. Non-parametric tests do not have an
#' analytically available posterior, and binomial/correlation tests use different
#' prior families.
#' @examples
#' \dontrun{
#' # Original study
#' original <- bfttest(rnorm(30, 0.5), prior.loc = 0, prior.r = 0.1)
#' fit <- extractPrior(original)
#' print(fit)
#'
#' # Use fitted parameters as informed prior for replication
#' replication <- bfttest(new_data, prior.loc = fit$location, prior.r = fit$scale)
#' }
#' @importFrom stats optim optimize uniroot dcauchy
#' @export
extractPrior <- function(x, n.points = 1000, delta.range = NULL) {

  # Validate input

  if (!inherits(x, "seqbf")) {
    stop("x must be a seqbf object created with bfttest()")
  }

  test_type <- x$`test type`

  if (test_type %in% c("binomial")) {
    stop("extractPrior() is not available for binomial tests. ",
         "Binomial tests use a Logistic prior, not a Cauchy prior.")
  }
  if (test_type %in% c("correlation")) {
    stop("extractPrior() is not available for correlation tests. ",
         "Correlation tests use a Beta prior, not a Cauchy prior.")
  }
  if (!test_type %in% c("one-sample", "paired", "independent")) {
    stop("Unsupported test type: ", test_type)
  }

  # Handle backward compatibility: old seqbf objects don't have 'parametric' field
  # If missing, assume parametric (old versions only supported parametric tests)
  is_parametric <- if (is.null(x$parametric)) TRUE else x$parametric

  if (!isTRUE(is_parametric)) {
    stop("extractPrior() requires a parametric t-test. ",
         "Non-parametric tests do not have an analytically available posterior.")
  }

  # Extract test parameters
  t_val <- tail(na.omit(x[[1]]), n = 1)
  prior_loc <- x$prior$location
  prior_scale <- x$prior$scale

  # Validate t-value
  if (length(t_val) == 0 || !is.finite(t_val)) {
    stop("Could not extract a valid t-statistic from the seqbf object. ",
         "t-value: ", if(length(t_val) == 0) "missing" else t_val)
  }

  # Check for extreme t-values that may cause numerical issues
  if (abs(t_val) > 100) {
    warning("Extremely large t-statistic (|t| = ", round(abs(t_val), 2), "). ",
            "This may cause numerical precision issues in posterior calculation.")
  }

  # Check for very large sample sizes that may cause numerical issues
  sample_size_check <- if (test_type == "independent") {
    min(x$`sample size`[1], x$`sample size`[2])
  } else {
    x$`sample size`[1]
  }

  if (sample_size_check > 1000) {
    warning("Very large sample size (n = ", sample_size_check, "). ",
            "With narrow priors, this may cause numerical underflow in posterior calculations. ",
            "Consider using a wider prior (e.g., prior.r = 1) if you encounter errors.")
  }

  # Validate prior parameters
  if (!is.finite(prior_loc)) {
    stop("Invalid prior location: ", prior_loc, ". ",
         "Prior location must be a finite number.")
  }
  if (!is.finite(prior_scale) || prior_scale <= 0) {
    stop("Invalid prior scale: ", prior_scale, ". ",
         "Prior scale must be a positive finite number.")
  }

  if (test_type == "independent") {
    n1 <- x$`sample size`[1]
    n2 <- x$`sample size`[2]
    indep <- TRUE
  } else {
    n1 <- if (length(x$`sample size`) == 1) x$`sample size` else x$`sample size`[1]
    n2 <- NULL
    indep <- FALSE
  }

  # Wrapper around posterior_t
  post_density <- function(delta) {
    posterior_t(delta, t = t_val, n1 = n1, n2 = n2,
                independentSamples = indep,
                prior.location = prior_loc,
                prior.scale = prior_scale,
                prior.df = 1)
  }

  # Estimate reasonable search interval based on t-value and sample size
  # For large samples, posterior concentrates near delta ≈ t/sqrt(n)
  neff <- if (indep) n1 * n2 / (n1 + n2) else n1
  delta_estimate <- t_val / sqrt(neff)

  # Test if posterior can be evaluated at the expected location
  test_vals <- post_density(c(0, delta_estimate, prior_loc))
  if (all(!is.finite(test_vals)) || all(test_vals == 0)) {
    stop("Posterior density function returns only zeros or non-finite values. ",
         "This indicates severe numerical underflow, likely due to:\n",
         "  - Very large sample size (n = ", neff, ")\n",
         "  - Very narrow prior (scale = ", prior_scale, ")\n",
         "Suggestions:\n",
         "  1. Use a wider prior: re-run bfttest() with prior.r >= 0.5\n",
         "  2. Use fewer observations for extractPrior()\n",
         "  3. This function may not be suitable for very large datasets")
  }

  # Use a wide interval centered on the estimate
  search_interval <- c(delta_estimate - 20, delta_estimate + 20)

  # Find posterior mode
  mode_result <- optimize(post_density, interval = search_interval, maximum = TRUE)
  post_mode <- mode_result$maximum
  half_max <- mode_result$objective / 2

  # Validate posterior mode
  if (!is.finite(post_mode)) {
    stop("Could not find a finite posterior mode. ",
         "This may indicate numerical issues with the posterior distribution.")
  }

  # Check if mode is at boundary (within 1% of interval endpoints)
  interval_width <- diff(search_interval)
  if (abs(post_mode - search_interval[1]) < 0.01 * interval_width ||
      abs(post_mode - search_interval[2]) < 0.01 * interval_width) {
    warning("Posterior mode found at search interval boundary (mode = ",
            round(post_mode, 4), ", interval = [",
            round(search_interval[1], 2), ", ", round(search_interval[2], 2), "]). ",
            "Results may be unreliable.")
  }

  # Find HWHM (half-width at half-maximum) for scale starting value
  hwhm <- tryCatch({
    root_result <- uniroot(function(d) post_density(d) - half_max,
                           interval = c(post_mode, post_mode + 10),
                           extendInt = "upX")
    abs(root_result$root - post_mode)
  }, error = function(e) {
    # Fallback: use prior scale as starting value
    prior_scale
  })

  # Validate HWHM: ensure it's positive and finite
  if (!is.finite(hwhm) || hwhm <= 0) {
    warning("Could not determine a valid HWHM. Using prior scale as fallback.")
    hwhm <- prior_scale
  }

  # Safety check: ensure hwhm has a minimum value to prevent numerical issues
  hwhm <- max(hwhm, 1e-6)

  # Determine delta range if not specified
  if (is.null(delta.range)) {
    delta.range <- c(post_mode - 5 * hwhm, post_mode + 5 * hwhm)
  }

  # Create evaluation grid
  delta_grid <- seq(delta.range[1], delta.range[2], length.out = n.points)
  post_vals <- post_density(delta_grid)

  # Validate posterior density values
  if (all(!is.finite(post_vals)) || all(post_vals == 0)) {
    stop("Could not evaluate posterior density. ",
         "All values are either non-finite or zero. ",
         "This may indicate numerical issues with the prior/posterior calculation.\n",
         "Diagnostics:\n",
         "  t-value: ", round(t_val, 4), "\n",
         "  Sample size: ", if(indep) paste0(n1, ", ", n2) else n1, "\n",
         "  Prior: location=", prior_loc, ", scale=", prior_scale, "\n",
         "  Posterior mode: ", round(post_mode, 4), "\n",
         "  Delta range: [", round(delta.range[1], 4), ", ", round(delta.range[2], 4), "]\n",
         "Try adjusting delta.range manually or check if your data produces extreme statistics.")
  }

  if (any(!is.finite(post_vals))) {
    warning("Some posterior density values are non-finite. ",
            "These will be set to zero for fitting.")
    post_vals[!is.finite(post_vals)] <- 0
  }

  # Fit Cauchy via minimizing integrated squared error
  # Optimize on log-scale for scale to enforce positivity
  objective <- function(par) {
    loc <- par[1]
    sc <- exp(par[2])  # log-scale
    cauchy_vals <- dcauchy(delta_grid, location = loc, scale = sc)
    sum((post_vals - cauchy_vals)^2)
  }

  fit <- optim(par = c(post_mode, log(hwhm)),
               fn = objective,
               method = "Nelder-Mead")

  fitted_loc <- fit$par[1]
  fitted_scale <- exp(fit$par[2])

  # Compute diagnostics
  cauchy_fitted <- dcauchy(delta_grid, location = fitted_loc, scale = fitted_scale)

  ss_res <- sum((post_vals - cauchy_fitted)^2)
  ss_tot <- sum((post_vals - mean(post_vals))^2)
  r_sq <- 1 - ss_res / ss_tot

  max_dev <- max(abs(post_vals - cauchy_fitted))

  result <- list(
    location = fitted_loc,
    scale = fitted_scale,
    r_squared = r_sq,
    max_abs_deviation = max_dev,
    prior = list(location = prior_loc, scale = prior_scale),
    test_type = test_type,
    sample_size = x$`sample size`
  )
  class(result) <- "cauchyFit"
  return(result)
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
