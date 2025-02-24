#' Bayesian Sequential Binomial Test
#'
#' Calculates Bayes Factors sequentially for dichotomous data, showing evidence evolution over time.
#' Uses 'proportionBF' from the BayesFactor package to compute Bayes Factors starting from a minimum
#' number of observations and updating with each new data point.
#'
#' @param data A vector containing binary data (0s and 1s)
#' @param p Probability of success under the null hypothesis (default: 0.5)
#' @param prior.r Scale parameter for the prior distribution (default: 0.1)
#' @param alternative Direction of alternative hypothesis: "two.sided", "greater", or "less"
#' @param nstart Minimum number of data points before first BF calculation (>= 2)
#' @param exact Logical. If TRUE, calculates BF for all points; if FALSE, uses steps for efficiency
#' @param nullInterval Optional vector of length 2 for interval hypothesis testing (deprecated)
#' @return A list of class "seqbf" containing:
#'   \itemize{
#'     \item probability of success: estimated probability
#'     \item p-value: from classical binomial test
#'     \item BF: vector of sequential Bayes Factors
#'     \item test type: "binomial"
#'     \item prior: list with prior distribution details
#'     \item sample size: total number of observations
#'     \item alternative: chosen alternative hypothesis
#'   }
#' @examples
#' data <- rbinom(100, 1, 0.7)
#' result <- bfbinom(data)
#' @importFrom BayesFactor proportionBF
#' @importFrom stats binom.test na.omit
#' @export

bfbinom <- function(data, p = 0.5, prior.r = 0.1, 
                    alternative = c("two.sided", "greater", "less"), 
                    nstart = 5, exact = TRUE, nullInterval = NULL) {
  
  # Input validation
  if (!is.numeric(data) || !all(data %in% c(0, 1, NA))) {
    stop("Data must be binary (0s and 1s)")
  }
  
  data <- na.omit(data)
  n_data <- length(data)
  
  # Validate parameters
  if (n_data < nstart) {
    stop(sprintf("Too few observations (%d). Need at least %d.", n_data, nstart))
  }
  if (!is.finite(prior.r) || prior.r <= 0) {
    stop("prior.r must be positive and finite")
  }
  if (!is.numeric(p) || p <= 0 || p >= 1) {
    stop("p must be between 0 and 1")
  }
  if (!(nstart == floor(nstart)) || nstart < 2) {
    stop("nstart must be an integer >= 2")
  }
  
  # Handle alternative hypothesis
  alternative <- match.arg(alternative)
  if (!is.null(nullInterval)) {
    warning("nullInterval is deprecated. Using 'alternative' parameter instead.")
  }
  
  # Set nullInterval based on alternative
  nullInterval <- switch(alternative,
                         "greater" = c(p, 1),
                         "less" = c(0, p),
                         "two.sided" = NULL
  )
  
  # Calculate steps for BF computation
  steps <- if (exact) {
    seq(nstart, n_data, 1)
  } else {
    .seqlast(nstart, n_data, .nstep(n_data)) #stepwise
  }
  
  # Initialize BF vector
  bf <- c(rep(1, nstart - 1), numeric(length(steps)))
  
  # Progress reporting
  message(sprintf("N = %d", n_data))
  message("Calculating Sequential Bayes Factors...")
  
  pb <- txtProgressBar(min = nstart, max = n_data, initial = nstart, style = 3)
  
  # Calculate BFs
  for (i in seq_along(steps)) {
    n <- steps[i]
    successes <- sum(data[1:n])
    
    bf_result <- tryCatch({
      tmpbfs <- proportionBF(successes, n, p = p, rscale = prior.r, 
                             nullInterval = nullInterval)
      exp(tmpbfs@bayesFactor$bf)[1]
    }, error = function(e) {
      warning(sprintf("Error at step %d: %s", n, e$message))
      NA
    })
    
    bf[nstart - 1 + i] <- bf_result
    setTxtProgressBar(pb, n)
  }
  
  close(pb)
  
  # Compute classical test
  orthodox_test <- binom.test(sum(data), n_data, p = p, alternative = alternative)
  
  # Prepare output
  bf_out <- list(
    "probability of success" = orthodox_test$estimate,
    "p-value" = orthodox_test$p.value,
    "BF" = bf,
    "test type" = "binomial",
    "prior" = list(
      "distribution" = "Logistic",
      "location" = 0,
      "scale" = prior.r
    ),
    "sample size" = n_data,
    "alternative" = alternative
  )
  
  # Report final results
  final_bf <- tail(bf, n = 1)
  if (final_bf > 1) {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f (probability of success = %.3f; p = %.3f)",
      final_bf, orthodox_test$estimate, orthodox_test$p.value
    ))
  } else {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f; BF01=%.3f (probability of success = %.3f; p = %.3f)",
      final_bf, 1/final_bf, orthodox_test$estimate, orthodox_test$p.value
    ))
  }
  
  class(bf_out) <- "seqbf"
  return(bf_out)
}


#' Bayesian Sequential t-Test
#'
#' Calculates sequential Bayes Factors for t-tests (one-sample, paired, or independent samples),
#' showing how evidence evolves as data accumulates. Uses the BFDA package for Bayes Factor
#' calculations.
#'
#' @param x A numeric vector of observations or a formula for independent samples test
#' @param y Optional numeric vector for paired samples test
#' @param formula Formula of form response ~ group where group has exactly 2 levels
#' @param data Optional data frame containing the variables in formula
#' @param alternative Direction of alternative hypothesis: "two.sided", "greater", or "less"
#' @param mu Null hypothesis value (default = 0)
#' @param prior.loc Location parameter for Cauchy prior (default = 0)
#' @param prior.r Scale parameter for Cauchy prior (default = 0.1)
#' @param nstart Minimum observations before first BF calculation ("auto" or >= 2)
#' @param exact Logical. If TRUE, calculates BF for all points
#' @return A list of class "seqbf" containing:
#'   \itemize{
#'     \item t-value: Sequential t-statistics
#'     \item p-value: Sequential p-values
#'     \item BF: Sequential Bayes Factors
#'     \item test type: "one-sample", "paired", or "independent"
#'     \item prior: List with prior distribution details
#'     \item sample size: Number of observations
#'     \item alternative: Chosen alternative hypothesis
#'   }
#' @examples
#' # One-sample test
#' x <- rnorm(30, 0.5, 1)
#' result1 <- bfttest(x, alternative = "greater")
#'
#' # Independent samples test
#' group <- rep(c("A", "B"), each = 15)
#' values <- c(rnorm(15), rnorm(15, 0.5))
#' df <- data.frame(values = values, group = group)
#' result2 <- bfttest(values ~ group, data = df)
#'
#' # Paired samples test
#' pre <- rnorm(20)
#' post <- pre + rnorm(20, 0.5)
#' result3 <- bfttest(pre, post)
#' @importFrom stats t.test var complete.cases na.omit
#' @export

bfttest <- function(x = NULL, y = NULL, formula = NULL, data = NULL, 
                    alternative = c("two.sided", "less", "greater"), 
                    mu = 0, prior.loc = 0, prior.r = 0.1, 
                    nstart = "auto", exact = TRUE) {
  
  # Match alternative argument
  alternative <- match.arg(alternative)
  
  # Validate nstart
  if (nstart != "auto") {
    if (!isTRUE(nstart == floor(nstart))) {
      stop("nstart must be an integer")
    }
    if (nstart < 2) {
      stop("nstart must be >= 2")
    }
  }
  
  # Handle formula input
  if (inherits(x, "formula")) {
    formula <- x
    x <- NULL
  }
  
  # Determine test type and validate inputs
  if (!is.null(formula)) {
    # Independent samples test
    if (!is.null(x)) {
      stop("Use either formula or x/y input, not both")
    }
    if (is.null(data)) {
      stop("Data frame required with formula input")
    }
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame")
    }
    
    # Extract variables from formula
    vars <- all.vars(formula)
    if (length(vars) != 2) {
      stop("Formula must contain exactly two variables")
    }
    
    # Clean data
    data <- data[complete.cases(data[vars]), ]
    if (nrow(data) == 0) {
      stop("No complete cases in data")
    }
    
    response_var <- vars[1]
    group_var <- vars[2]
    
    # Validate group variable
    groups <- unique(data[[group_var]])
    if (length(groups) != 2) {
      stop("Group variable must have exactly 2 levels")
    }
    
    test_type <- "independent"
    sample_size <- table(data[[group_var]])
    total_sample_size <- sum(sample_size)
    
    # Determine starting point
    if (nstart == "auto") {
      nstart <- .determine_min_n_independent(data, group_var, vars[1], alternative, prior.loc, prior.r)
    }
    
  } else if (!is.null(y)) {
    # Paired samples test
    if (length(x) != length(y)) {
      stop("x and y must have same length for paired test")
    }
    complete_cases <- complete.cases(x, y)
    x <- x[complete_cases]
    y <- y[complete_cases]
    
    test_type <- "paired"
    sample_size <- length(x)
    total_sample_size <- sample_size
    
    if (nstart == "auto") {
      nstart <- .determine_min_n_paired(x, y, alternative, prior.loc, prior.r)
    }
    
  } else if (!is.null(x)) {
    # One sample test
    x <- na.omit(x)
    if (length(x) == 0) {
      stop("No valid observations in x")
    }
    
    test_type <- "one-sample"
    sample_size <- length(x)
    total_sample_size <- sample_size
    
    if (nstart == "auto") {
      nstart <- .determine_min_n_one_sample(x, alternative, prior.loc, prior.r)
    }
    
  } else {
    stop("No valid input provided")
  }
  
  # Calculate sequential steps
  steps <- if (exact) {
    seq(nstart, total_sample_size, 1)
  } else {
    .seqlast(nstart, total_sample_size, .nstep(total_sample_size))  # stepwise
  }
  
  # Initialize vectors for results
  t_values <- rep(NA, total_sample_size)
  p_values <- rep(NA, total_sample_size)
  bf_values <- rep(NA, total_sample_size)
  
  # Pre-fill BF values for 1:(nstart-1) with 1
  if (nstart > 1) {
    bf_values[1:(nstart - 1)] <- 1
  }
  
  # Progress reporting
  message(sprintf("%s test (N = %d%s)",
                  .capitalize(test_type),
                  total_sample_size,
                  ifelse(test_type == "independent",
                         sprintf(" [%d + %d]", sample_size[1], sample_size[2]),
                         "")))
  message("Calculating Sequential Bayes Factors...")
  
  # Calculate sequential BFs
  pb <- txtProgressBar(min = nstart, max = total_sample_size, initial = nstart, style = 3)
  
  for (i in seq_along(steps)) {
    n <- steps[i]
    
    result <- tryCatch({
      if (test_type == "independent") {
        subset_data <- data[1:n, ]
        t_result <- t.test(formula, data = subset_data, 
                           alternative = alternative, var.equal = TRUE)
        n1 <- table(subset_data[[group_var]])[1]
        n2 <- table(subset_data[[group_var]])[2]
        bf_result <- bf10_t(t = t_result$statistic,
                            n1 = n1, n2 = n2,
                            independentSamples = TRUE,
                            prior.location = prior.loc,
                            prior.scale = prior.r,
                            prior.df = 1)
      } else if (test_type == "paired") {
        t_result <- t.test(x[1:n], y[1:n],
                           alternative = alternative,
                           paired = TRUE)
        bf_result <- bf10_t(t = t_result$statistic,
                            n1 = n,
                            independentSamples = FALSE,
                            prior.location = prior.loc,
                            prior.scale = prior.r,
                            prior.df = 1)
      } else {
        t_result <- t.test(x[1:n],
                           alternative = alternative,
                           mu = mu)
        bf_result <- bf10_t(t = t_result$statistic,
                            n1 = n,
                            prior.location = prior.loc,
                            prior.scale = prior.r,
                            prior.df = 1)
      }
      list(t = unname(t_result$statistic),
           p = t_result$p.value,
           bf = .get_directional_bf(bf_result, alternative))
    }, error = function(e) {
      warning(sprintf("Error at step %d: %s", n, e$message))
      list(t = NA, p = NA, bf = NA)
    })
    
    t_values[n] <- result$t
    p_values[n] <- result$p
    bf_values[n] <- result$bf
    setTxtProgressBar(pb, n)
  }
  
  close(pb)
  
  # Prepare output
  bf_out <- list(
    "t-value" = t_values,
    "p-value" = p_values,
    "BF" = bf_values,
    "test type" = test_type,
    "prior" = list(
      "distribution" = "Cauchy",
      "location" = prior.loc,
      "scale" = prior.r
    ),
    "sample size" = sample_size,
    "alternative" = alternative
  )
  
  # Get final results
  final_bf <- tail(na.omit(bf_out$BF), n = 1)
  final_t <- tail(na.omit(bf_out$`t-value`), n = 1)
  final_p <- tail(na.omit(bf_out$`p-value`), n = 1)
  
  # Report final results
  if (is.na(final_bf)) {
    message("Final Bayes Factor could not be calculated")
  } else if (final_bf > 1) {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f (t = %.3f; p = %.3f)",
      final_bf, final_t, final_p
    ))
  } else if (final_bf == 1) {
    message(sprintf(
      "Final Bayes Factor: BF10 = 1 (no evidence for either hypothesis; t = %.3f; p = %.3f)",
      final_t, final_p
    ))
  } else {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f; BF01 = %.3f (t = %.3f; p = %.3f)",
      final_bf, 1/final_bf, final_t, final_p
    ))
  }
  
  class(bf_out) <- "seqbf"
  return(bf_out)
}


#' Bayesian Sequential Correlation Test
#'
#' Calculates sequential Bayes Factors for correlation between two continuous variables,
#' showing how evidence evolves as data accumulates. Uses 'correlationBF' from the
#' BayesFactor package.
#'
#' @param x A vector containing continuous data
#' @param y A vector containing continuous data
#' @param alternative Direction of alternative hypothesis: "two.sided", "greater", or "less"
#' @param prior.r Scale parameter for the prior distribution (default: 0.1)
#' @param nstart Minimum number of data points before first BF calculation (>= 2)
#' @param exact Logical. If TRUE, calculates BF for all points; if FALSE, uses steps for efficiency
#' @return A list of class "seqbf" containing:
#'   \itemize{
#'     \item r: Pearson correlation coefficient
#'     \item p-value: from classical correlation test
#'     \item BF: vector of sequential Bayes Factors
#'     \item test type: "correlation"
#'     \item prior: list with prior distribution details
#'     \item sample size: total number of observations
#'     \item alternative: chosen alternative hypothesis
#'   }
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100, 0, 0.5)  # Correlated data
#' result <- bfcor(x, y, alternative = "greater")
#' @importFrom BayesFactor correlationBF
#' @importFrom stats cor.test na.omit
#' @export

bfcor <- function(x, y, alternative = c("two.sided", "greater", "less"), 
                  prior.r = 0.1, nstart = 5, exact = TRUE) {
  
  # Match alternative argument
  alternative <- match.arg(alternative)
  
  # Handle missing values
  complete_cases <- complete.cases(x, y)
  x <- x[complete_cases]
  y <- y[complete_cases]
  
  # Input validation
  if (length(x) == 0 || length(y) == 0) {
    stop("Data has no valid observations after removing missing values")
  }
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Both x and y must be numeric vectors")
  }
  if (!all(is.finite(x)) || !all(is.finite(y))) {
    stop("Data must be finite")
  }
  if (length(y) != length(x)) {
    stop("Length of y and x must be the same")
  }
  if (!isTRUE(nstart == floor(nstart)) || nstart < 2) {
    stop("nstart must be an integer >= 2")
  }
  if (!is.finite(prior.r) || prior.r <= 0) {
    stop("prior.r must be positive and finite")
  }
  
  # Set nullInterval based on alternative
  nullInterval <- switch(alternative,
                         "greater" = c(-1, 0),   # Null: correlation <= 0
                         "less" = c(0, 1),       # Null: correlation >= 0
                         "two.sided" = 0
  )
  
  n_data <- length(x)
  
  if (n_data < nstart) {
    stop(sprintf("Too few observations (%d). Need at least %d.", n_data, nstart))
  }
  
  # Calculate steps for BF computation
  steps <- if (exact) {
    seq(nstart, n_data, 1)
  } else {
    # Create reasonable steps for larger datasets
    step_size <- max(1, floor(n_data / 100))
    seq(nstart, n_data, step_size)
  }
  
  # Initialize BF vector
  bf <- c(rep(1, nstart - 1), numeric(length(steps)))
  
  # Progress reporting
  message(sprintf("N = %d", n_data))
  message("Calculating Sequential Bayes Factors...")
  
  pb <- txtProgressBar(min = nstart, max = n_data, initial = nstart, style = 3)
  
  # Calculate BFs
  for (i in seq_along(steps)) {
    n <- steps[i]
    bf_result <- tryCatch({
      tmpbfs <- correlationBF(x[1:n], y[1:n], rscale = prior.r, 
                              nullInterval = nullInterval)
      exp(tmpbfs@bayesFactor$bf[2])
    }, error = function(e) {
      warning(sprintf("Error at step %d: %s", n, e$message))
      NA
    })
    
    bf[nstart - 1 + i] <- bf_result
    setTxtProgressBar(pb, n)
  }
  
  close(pb)
  
  # Compute classical test
  orthodox_test <- cor.test(x, y, alternative = alternative)
  
  # Prepare output
  bf_out <- list(
    "r" = orthodox_test$estimate,
    "p-value" = orthodox_test$p.value,
    "BF" = bf,
    "test type" = "correlation",
    "prior" = list(
      "distribution" = "Beta",
      "location" = 0,
      "scale" = prior.r
    ),
    "sample size" = n_data,
    "alternative" = alternative
  )
  
  # Get final BF
  final_bf <- tail(bf, n = 1)
  
  # Report final results with improved BF interpretation
  if (is.na(final_bf)) {
    message("Final Bayes Factor could not be calculated")
  } else if (final_bf > 1) {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f (r = %.3f; p = %.3f)",
      final_bf, orthodox_test$estimate, orthodox_test$p.value
    ))
  } else if (final_bf == 1) {
    message(sprintf(
      "Final Bayes Factor: BF10 = 1 (no evidence for either hypothesis; r = %.3f; p = %.3f)",
      orthodox_test$estimate, orthodox_test$p.value
    ))
  } else {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f; BF01 = %.3f (r = %.3f; p = %.3f)",
      final_bf, 1/final_bf, orthodox_test$estimate, orthodox_test$p.value
    ))
  }
  
  class(bf_out) <- "seqbf"
  return(bf_out)
}


#' Robustness Analysis for Bayesian Sequential tests
#'
#' This function compares Bayes Factors of different priors. Returns a list containing BFMatrix, a data frame comprising the Bayes Factors of all possible combinations of the prior parameters (prior.location and prior.width)
#'
#'
#' @param x A seqbf object created with bfttest(), bfbinom(), or bfcor()
#' @param informed Logical. If you don't want to specify prior locations and scales manually, informed = FALSE will test multiple prior widths and location 0, while informed = TRUE will test different locations and widths.
#' @param prior.loc Range of locations of the cauchy distributed prior function.
#' @param prior.r Range of scales of the cauchy distributed prior function.
#' @examples
#' bfrInformed <- bfRobustness(seqbf) 
#' bfrUninformed <- bfRobustness(seqbf, informed = FALSE) 
#' print(bfrInformed)
#' plot(bfrInformed)
#' @export

bfRobustness <- function(x = NULL, informed = TRUE, prior.loc = NULL, prior.r = NULL){
  
  if(inherits(x,"seqbf") == FALSE) stop("Please provide a seqbf object created with bfttest() or bfbinom() or bfcor()")
  
  # specify prior parameters if not set by user
  
  if (is.null(prior.loc)) {
    if (informed == TRUE) prior.loc <- seq(0,1,0.05)
    else prior.loc <- 0
  }
  if (is.null(prior.r)) {
    if (informed == TRUE) prior.r <- seq(0.05, 1, 0.05)
    else prior.r <- seq(0.01, 1.41, 0.01)
  }
  
  if(x$alternative == "two.sided") aa <- 1
  else if(x$alternative == "greater") aa <- 2
  else if(x$alternative == "less") aa <- 3
  
  t <- tail(x$`t-value`, n=1)
  grid <- expand.grid(prior.loc, prior.r)
  
  if(x$`test type` == "independent"){
    
    n1 <- x$`sample size`[1]
    n2 <- x$`sample size`[2]
    grid$bf <- pbapply::pbapply(grid,1,function(g) bf10_t(t = t, n1 = n1, n2 = n2, independentSamples = T, prior.location = g[1], prior.scale = g[2], prior.df = 1)[[aa]])
    
  } else if(x$`test type` == "paired"){
    
    grid$bf <- pbapply::pbapply(grid,1,function(g) bf10_t(t = t, n1 = x$`sample size`, independentSamples = F, prior.location = g[1], prior.scale = g[2], prior.df = 1)[[aa]])
    
  } else if(x$`test type` == "one-sample"){
    
    grid$bf <- pbapply::pbapply(grid,1,function(g) bf10_t(t = t, n1 = x$`sample size`, prior.location = g[1], prior.scale = g[2], prior.df = 1)[[aa]])
    
  }
  
  
  names(grid)[1] <- "prior.loc"
  names(grid)[2] <- "prior.r"
  
  bft.out <- list("BFMatrix" = grid, "test type" = x$`test type`, "prior" = list(x$prior[[1]], "prior location" = x$prior[[2]], "prior scale" = x$prior[[3]]), "sample size" = x$`sample size`, "alternative" = x$alternative)
  
  cat("Highest BF = ", round(max(grid$bf),2), " with prior: Cauchy(",grid$prior.loc[grid$bf==max(grid$bf)],", ",grid$prior.r[grid$bf==max(grid$bf)],")", sep="")
  
  class(bft.out) <- "bfRobustness"
  return(bft.out)
  
}


#' @export
#' @method print seqbf
print.seqbf <- function(x, ...) {
    
    cat(paste0("
  Sequential Bayesian Testing
  --------------------------------	
  Test type: ", x$`test type`, "
  Sample size: ", paste(x$`sample size`, collapse=", "), "
  Final Bayes Factor: BF10=", round(tail(x$BF, n=1), 3), "; BF01=", round(1/tail(x$BF, n=1), 3), "
  Parameter prior: ", x$prior[[1]],"(",x$prior[[2]],", ",x$prior[[3]],")", "
  Directionality of H1 analysis prior: ", x$alternative, "
  Orthodox Test: ",names(x)[1],"=",round(tail(x[[1]],n=1), 3),"; p=",round(tail(x$`p-value`,n=1), 3), "
              \n"
    ))
}


#' @export
#' @method print bfRobustness
print.bfRobustness <- function(x, ...) {
  
  grid <- x$BFMatrix
  
  cat(paste0("
  Prior Robustness Analysis
  --------------------------------	
  Test type: ", x$`test type`, "
  Sample size: ", paste(x$`sample size`, collapse=", "), "
  Tested priors: 
  -- distribution: ", x$prior[1], "
  -- location: ", paste(unique(grid$prior.loc), collapse=","), "
  -- scale: ", paste(unique(grid$prior.r), collapse=","), "
  Highest Bayes Factor: ", round(max(grid$bf),3), " with prior: Cauchy(",grid$prior.loc[grid$bf==max(grid$bf)],", ",grid$prior.r[grid$bf==max(grid$bf)],")
  Lowest Bayes Factor: ", round(min(grid$bf),3), " with prior: Cauchy(",grid$prior.loc[grid$bf==min(grid$bf)],", ",grid$prior.r[grid$bf==min(grid$bf)],")
  Median Bayes Factor: ", round(median(grid$bf),3), " 
              \n"
  ))
}


# Helper functions
.seqlast <- function(from, to, by){
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

.nstep <- function(size){
  step <- floor((size-1)/100)+1
  return(step)
}

# Try running bfttest with specific n
.try_bft_calculation <- function(n, ...) {
  nstart <- n - 1
  tryCatch({
    
    sink(tempfile())  # Redirect output to a temp file
    result <- suppressMessages(suppressWarnings(
      bfttest(..., nstart = nstart, exact = TRUE)
    ))
    sink()  # Restore normal output
    
    
    # Check if we got valid BF values
    #!is.na(tail(result$BF, 1))
    !is.na(result$BF[length(result$BF) - 1])
  }, error = function(e) {
    #warning(sprintf("Error at step %d: %s", nstart, e$message))
    FALSE
  })
}

# Determine minimum n for independent samples test
.determine_min_n_independent <- function(data, group_var, response_var, 
                                         alternative = "two.sided", prior.loc = 0, prior.r = 0.1) {
  n <- 3
  nstart <- n - 1
  while (n <= nrow(data)) {
    subset <- data[1:n, ]
    subsetstart <- data[1:nstart, ]
    # Check basic conditions first
    if (length(unique(subsetstart[[group_var]])) == 2 &&
        var(subsetstart[[response_var]]) > 0) {
      # Create formula for the test
      formula <- as.formula(paste(response_var, "~", group_var))
      # Try running bfttest
      if (.try_bft_calculation(n, 
                               formula = formula, 
                               data = subset,
                               alternative = alternative,
                               prior.loc = prior.loc,
                               prior.r = prior.r)) {
        return(nstart)
      }
    }
    n <- n + 1
    nstart <- n - 1
  }
  stop("Could not find valid starting point for Bayes Factor calculation")
}

# Determine minimum n for paired samples test
.determine_min_n_paired <- function(x, y, alternative = "two.sided", prior.loc = 0, prior.r = 0.1) {
  n <- 3
  nstart <- n - 1
  while (n <= length(x)) {
    if (var(x[1:nstart] - y[1:nstart]) > 0) {
      # Try running bfttest
      if (.try_bft_calculation(n,
                               x = x[1:n],
                               y = y[1:n],
                               alternative = alternative,
                               prior.loc = prior.loc,
                               prior.r = prior.r)) {
        return(nstart)
      }
    }
    n <- n + 1
    nstart <- n - 1
  }
  stop("Could not find valid starting point for Bayes Factor calculation")
}

# Determine minimum n for one sample test
.determine_min_n_one_sample <- function(x, alternative = "two.sided", prior.loc = 0, prior.r = 0.1) {
  n <- 3
  nstart <- n - 1
  while (n <= length(x)) {
    if (var(x[1:nstart]) > 0) {
      # Try running bfttest
      if (.try_bft_calculation(n,
                               x = x[1:n],
                               alternative = alternative,
                               prior.loc = prior.loc,
                               prior.r = prior.r)) {
        return(nstart)
      }
    }
    n <- n + 1
    nstart <- n - 1
  }
  stop("Could not find valid starting point for Bayes Factor calculation")
}

# Get directional BF based on alternative
.get_directional_bf <- function(bf_result, alternative) {
  switch(alternative,
         "less" = bf_result$BFmin0,
         "greater" = bf_result$BFplus0,
         "two.sided" = bf_result$BF10)
}

# Capitalize first letter
.capitalize <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
