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
#'     \item data: list with element \code{x} (the cleaned binary data vector)
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
    .seqlast(nstart, n_data, .nstep(n_data))
  }
  
  # Initialize BF vector for full sample size
  bf <- rep(NA, n_data)
  
  # Pre-fill BF values for 1:(nstart-1) with 1
  if (nstart > 1) {
    bf[1:(nstart - 1)] <- 1
  }
  
  # Progress reporting
  message(sprintf("Binomial test (N = %d)", n_data))
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
    
    bf[n] <- bf_result  # Assign to position n, not sequential position
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
      "location" = log(p / (1 - p)),  # logit(p)
      "scale" = prior.r
    ),
    "sample size" = n_data,
    "alternative" = alternative,
    "null_p" = p,
    "data_type" = "binary",
    "data" = list(x = data)
  )
  
  # Report final results
  final_bf <- tail(na.omit(bf), n = 1)  # Use na.omit here
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


#' Bayesian Sequential t-Test or Non-Parametric Test
#'
#' Calculates sequential Bayes Factors for parametric t-tests or non-parametric 
#' rank-based tests (Wilcoxon/Mann-Whitney), showing how evidence evolves as data accumulates.
#' For parametric tests, uses BFDA package. For non-parametric tests, uses van Doorn et al.'s
#' latent normal approach with Cauchy priors on effect size.
#'
#' @param x A numeric vector of observations or a formula for independent samples test
#' @param y Optional numeric vector for paired samples test
#' @param formula Formula of form response ~ group where group has exactly 2 levels
#' @param data Optional data frame containing the variables in formula
#' @param alternative Direction of alternative hypothesis: "two.sided", "greater", or "less"
#' @param mu Null hypothesis value (default = 0, only for one-sample parametric test)
#' @param parametric Logical. If TRUE (default), uses t-test. If FALSE, uses Wilcoxon/Mann-Whitney
#' @param prior.loc Location parameter for Cauchy prior (default = 0)
#' @param prior.r Scale parameter for Cauchy prior (default = 0.1). Use sqrt(2)/2 for default JZS prior.
#' @param nstart Minimum observations before first BF calculation ("auto" or >= 2)
#' @param nsamples Number of MCMC samples for non-parametric tests (default = adaptive sampling based on N)
#' @param exact Logical. If TRUE, calculates BF for all points
#' @param parallel Logical. If TRUE (default), uses parallel processing across available
#'   CPU cores (detectCores() - 1), providing substantial speedup especially for
#'   non-parametric tests with large samples. Set to FALSE for reproducibility across
#'   sessions or when debugging errors.
#' @param data_type Type of input data: "unknown" (default, auto-detected), "summed_bits"
#'   (sum of binary bits), "continuous", or "integer".
#' @param n_bits Number of bits summed per observation (only relevant for data_type = "summed_bits").
#' @param bit_probability Probability of a 1-bit (only relevant for data_type = "summed_bits",
#'   default = 0.5).
#' @return A list of class "seqbf" containing:
#'   \itemize{
#'     \item t-value or W-value: Sequential test statistics
#'     \item p-value: Sequential p-values
#'     \item BF: Sequential Bayes Factors
#'     \item delta: Sequential delta estimates (posterior medians for effect size)
#'     \item delta.lower: Lower bounds of 95% credible intervals for delta
#'     \item delta.upper: Upper bounds of 95% credible intervals for delta
#'     \item test type: "one-sample", "paired", or "independent"
#'     \item parametric: Logical indicating test type
#'     \item prior: List with prior distribution details
#'     \item sample size: Number of observations
#'     \item alternative: Chosen alternative hypothesis
#'     \item data: list with \code{x} (one-sample/paired) or \code{x} and \code{y} (paired)
#'       or \code{data} and \code{formula} (independent samples)
#'   }
#' @note Adaptive sampling uses 
#' - 250 samples for intermediate steps and 1000 samples for final BF for N <= 500
#' - 150 samples for intermediate steps and 750 samples for final BF for 500 < N <= 2000
#' - 100 samples for intermediate steps and 500 samples for final BF for N > 2000
#' @note Parallel processing can significantly speed up computation, especially for non-parametric 
#' tests with large samples. Set parallel = FALSE for reproducible results or debugging.
#' @examples
#' \dontrun{
#' # Parametric one-sample test
#' x <- rnorm(30, 0.5, 1)
#' result1 <- bfttest(x, alternative = "greater")
#'
#' # Non-parametric independent samples test with parallel processing (default, fast)
#' group <- rep(c("A", "B"), each = 15)
#' values <- c(rnorm(15), rnorm(15, 0.5))
#' df <- data.frame(values = values, group = group)
#' result2 <- bfttest(values ~ group, data = df, parametric = FALSE)
#'
#' # For reproducibility or debugging, disable parallel processing
#' result2b <- bfttest(values ~ group, data = df, parametric = FALSE, parallel = FALSE)
#'
#' # Non-parametric paired samples test
#' pre <- rnorm(20)
#' post <- pre + rnorm(20, 0.5)
#' result3 <- bfttest(pre, post, parametric = FALSE)
#' }
#' @importFrom stats t.test wilcox.test var complete.cases na.omit quantile
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom pbapply pblapply
#' @export

bfttest <- function(x = NULL, y = NULL, formula = NULL, data = NULL,
                    alternative = c("two.sided", "less", "greater"),
                    mu = 0, parametric = TRUE,
                    prior.loc = 0,
                    prior.r = 0.1,
                    nstart = "auto",
                    nsamples = "auto",
                    exact = TRUE,
                    parallel = TRUE,
                    # NEW PARAMETERS:
                    data_type = c("unknown", "summed_bits", "continuous", "integer"),
                    n_bits = NULL,
                    bit_probability = NULL) {
  
  # Match alternative argument
  alternative <- match.arg(alternative)

  # Match data_type argument
  data_type <- match.arg(data_type)

  # Set defaults based on data_type
  if (!is.null(n_bits) && data_type == "unknown") {
    data_type <- "summed_bits"
  }

  # Set default bit_probability
  if (is.null(bit_probability)) {
    bit_probability <- 0.5
  }

  # Validate parallel argument
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("parallel must be TRUE or FALSE")
  }
  
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
      if (parametric) {
        nstart <- .determine_min_n_independent(data, group_var, vars[1], alternative, prior.loc, prior.r)
      } else {
        # For non-parametric tests, find the first n where both groups have sufficient data
        min_per_group <- 5  # Minimum observations per group for stable MCMC
        
        for (n in seq_len(total_sample_size)) {
          group_counts <- table(data[[group_var]][1:n])
          if (length(group_counts) == 2 && all(group_counts >= min_per_group)) {
            nstart <- n
            break
          }
        }
        
        if (is.character(nstart) && nstart == "auto") {
          stop("Could not find valid starting point where both groups have sufficient data")
        }
      }
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
      if (parametric) {
        nstart <- .determine_min_n_paired(x, y, alternative, prior.loc, prior.r)
      } else {
        nstart <- 10  # Minimum for stable MCMC
      }
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
      if (parametric) {
        nstart <- .determine_min_n_one_sample(x, alternative, prior.loc, prior.r)
      } else {
        nstart <- 10  # Minimum for stable MCMC
      }
    }
    
  } else {
    stop("No valid input provided")
  }

  # Auto-detect data_type from the cleaned data if still unknown
  if (data_type == "unknown") {
    vals <- if (test_type == "independent") {
      data[[response_var]]
    } else if (test_type == "paired") {
      c(x, y)
    } else {
      x
    }

    data_type <- if (all(vals == floor(vals), na.rm = TRUE)) {
      if (min(vals, na.rm = TRUE) >= 0) "summed_bits" else "integer"
    } else {
      "continuous"
    }

    message(sprintf(
      "data_type auto-detected as '%s'. Set data_type= explicitly in bfttest() to suppress this message.",
      data_type
    ))
  }

  # Validate nstart doesn't exceed sample size
  if (nstart > total_sample_size) {
    stop(sprintf("nstart (%d) cannot exceed total sample size (%d)", 
                 nstart, total_sample_size))
  }
  
  # Calculate sequential steps
  steps <- if (exact) {
    seq(nstart, total_sample_size, 1)
  } else {
    .seqlast(nstart, total_sample_size, .nstep(total_sample_size))
  }
  
  # Determine MCMC sampling strategy upfront (for non-parametric tests)
  use_adaptive_sampling <- FALSE
  intermediate_samples <- NULL
  final_samples <- NULL
  
  if (!parametric) {
    if (nsamples == "auto") {
      use_adaptive_sampling <- TRUE
      # Get intermediate and final sample counts
      intermediate_samples <- .adaptive_nsamples(total_sample_size, is_final = FALSE)
      final_samples <- .adaptive_nsamples(total_sample_size, is_final = TRUE)
    } else {
      # Validate nsamples
      if (!isTRUE(nsamples == floor(nsamples)) || nsamples <= 0) {
        stop("nsamples must be a positive integer or 'auto'")
      }
      intermediate_samples <- nsamples
      final_samples <- nsamples
    }
  }
  
  # Initialize vectors for results
  stat_values <- rep(NA, total_sample_size)
  p_values <- rep(NA, total_sample_size)
  bf_values <- rep(NA, total_sample_size)
  delta_values <- rep(NA, total_sample_size)
  delta_lower <- rep(NA, total_sample_size)
  delta_upper <- rep(NA, total_sample_size)
  
  # Pre-fill BF values for 1:(nstart-1) with 1
  if (nstart > 1) {
    bf_values[1:(nstart - 1)] <- 1
  }
  
  # Progress reporting
  test_name <- if (parametric) {
    sprintf("%s t-test", .capitalize(test_type))
  } else {
    if (test_type == "independent") {
      "Mann-Whitney U test"
    } else {
      sprintf("%s Wilcoxon signed-rank test", .capitalize(test_type))
    }
  }
  
  message(sprintf("%s (N = %d%s)",
                  test_name,
                  total_sample_size,
                  ifelse(test_type == "independent",
                         sprintf(" [%d + %d]", sample_size[1], sample_size[2]),
                         "")))
  
  # Display MCMC sampling strategy for non-parametric tests
  if (!parametric) {
    if (use_adaptive_sampling) {
      message(sprintf("MCMC samples: %d (intermediate) \u2192 %d (final)",
                      intermediate_samples, final_samples))
    } else {
      message(sprintf("MCMC samples: %d", final_samples))
    }
  }
  
  # Display parallel processing info
  if (parallel) {
    cores <- parallel::detectCores() - 1
    message(sprintf("Using parallel processing with %d cores", cores))
  }
  
  message("Calculating Sequential Bayes Factors...")
  
  # Define worker function for parallel/sequential execution
  calculate_step <- function(i) {
    n <- steps[i]
    is_final <- (i == length(steps))
    
    # Determine current sample count for this iteration
    if (!parametric) {
      current_nsamples <- if (is_final) final_samples else intermediate_samples
    }
    
    result <- tryCatch({
      if (parametric) {
        # PARAMETRIC TESTS
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
          # Approximate delta from t-statistic
          delta_est <- unname(t_result$statistic) * sqrt(1/n1 + 1/n2)
          
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
          # Approximate delta from t-statistic
          delta_est <- unname(t_result$statistic) / sqrt(n)
          
        } else {  # one-sample
          t_result <- t.test(x[1:n],
                             alternative = alternative,
                             mu = mu)
          bf_result <- bf10_t(t = t_result$statistic,
                              n1 = n,
                              prior.location = prior.loc,
                              prior.scale = prior.r,
                              prior.df = 1)
          # Approximate delta from t-statistic
          delta_est <- unname(t_result$statistic) / sqrt(n)
        }
        
        list(
          n = n,
          stat = unname(t_result$statistic),
          p = t_result$p.value,
          bf = .get_directional_bf(bf_result, alternative),
          delta = delta_est,
          delta_lower = NA,
          delta_upper = NA
        )
        
      } else {
        # NON-PARAMETRIC TESTS
        if (test_type == "independent") {
          # Mann-Whitney U test
          subset_data <- data[1:n, ]
          group1_data <- subset_data[[response_var]][subset_data[[group_var]] == groups[1]]
          group2_data <- subset_data[[response_var]][subset_data[[group_var]] == groups[2]]
          
          # Frequentist test for p-value
          w_result <- wilcox.test(group1_data, group2_data,
                                  alternative = alternative,
                                  exact = FALSE)
          
          # Bayesian rank sum test
          mcmc_result <- rankSumGibbsSampler(
            xVals = group1_data, 
            yVals = group2_data,
            nSamples = current_nsamples,
            cauchyPriorParameter = prior.r,
            progBar = FALSE
          )
          
          # Compute BF
          bf_result <- .get_bf_nonparametric(
            mcmc_result$deltaSamples, 
            alternative,
            prior.r
          )
          
          # Summarize delta
          delta_quant <- quantile(mcmc_result$deltaSamples, 
                                  probs = c(0.025, 0.5, 0.975))
          
          list(
            n = n,
            stat = unname(w_result$statistic),
            p = w_result$p.value,
            bf = bf_result,
            delta = delta_quant[2],
            delta_lower = delta_quant[1],
            delta_upper = delta_quant[3]
          )
          
        } else {
          # Wilcoxon signed-rank test (paired or one-sample)
          if (test_type == "paired") {
            diff_data <- x[1:n] - y[1:n]
            w_result <- wilcox.test(x[1:n], y[1:n],
                                    alternative = alternative,
                                    paired = TRUE,
                                    exact = FALSE)
          } else {  # one-sample
            diff_data <- x[1:n] - mu
            w_result <- wilcox.test(x[1:n],
                                    alternative = alternative,
                                    mu = mu,
                                    exact = FALSE)
          }
          
          # Bayesian signed rank test
          mcmc_result <- signRankGibbsSampler(
            xVals = diff_data,
            yVals = NULL,
            nSamples = current_nsamples,
            cauchyPriorParameter = prior.r,
            progBar = FALSE
          )
          
          # Compute BF
          bf_result <- .get_bf_nonparametric(
            mcmc_result$deltaSamples, 
            alternative,
            prior.r
          )
          
          # Summarize delta
          delta_quant <- quantile(mcmc_result$deltaSamples, 
                                  probs = c(0.025, 0.5, 0.975))
          
          list(
            n = n,
            stat = unname(w_result$statistic),
            p = w_result$p.value,
            bf = bf_result,
            delta = delta_quant[2],
            delta_lower = delta_quant[1],
            delta_upper = delta_quant[3]
          )
        }
      }
    }, error = function(e) {
      warning(sprintf("Error at step %d: %s", n, e$message))
      list(n = n, stat = NA, p = NA, bf = NA, delta = NA, 
           delta_lower = NA, delta_upper = NA)
    })
    
    return(result)
  }
  
  # Execute calculations (parallel or sequential)
  if (parallel && length(steps) > 1) {
    # Parallel execution
    cl <- parallel::makeCluster(cores)
    
    # Export necessary objects and functions
    parallel::clusterEvalQ(cl, {
      library(changeofevidence)
    })
    
    parallel::clusterExport(cl, 
                            c("steps", "test_type", "parametric", 
                              "alternative", "prior.loc", "prior.r", "mu",
                              "intermediate_samples", "final_samples"),
                            envir = environment())
    
    # Export test-specific data
    if (test_type == "independent") {
      parallel::clusterExport(cl, c("data", "formula", "response_var", "group_var", "groups"), 
                              envir = environment())
    } else {
      parallel::clusterExport(cl, c("x", "y"), envir = environment())
    }
    
    # Run parallel calculations with progress bar
    results_list <- pbapply::pblapply(seq_along(steps), calculate_step, cl = cl)
    
    parallel::stopCluster(cl)
    
  } else {
    # Sequential execution with progress bar
    results_list <- pbapply::pblapply(seq_along(steps), calculate_step)
  }
  
  # Organize results into vectors
  for (result in results_list) {
    n <- result$n
    stat_values[n] <- result$stat
    p_values[n] <- result$p
    bf_values[n] <- result$bf
    delta_values[n] <- result$delta
    delta_lower[n] <- result$delta_lower
    delta_upper[n] <- result$delta_upper
  }
  
  # Prepare output
  stat_name <- if (parametric) "t-value" else "W-value"
  
  # Store tested data consistently
  tested_data <- if (test_type == "independent") {
    list(data = data, formula = formula)
  } else if (test_type == "paired") {
    list(x = x, y = y)
  } else {
    list(x = x)
  }

  bf_out <- list(
    stat_name = stat_values,
    "p-value" = p_values,
    "BF" = bf_values,
    "delta" = delta_values,
    "delta.lower" = delta_lower,
    "delta.upper" = delta_upper,
    "test type" = test_type,
    "parametric" = parametric,
    "prior" = list(
      "distribution" = "Cauchy",
      "location" = prior.loc,
      "scale" = prior.r
    ),
    "sample size" = sample_size,
    "alternative" = alternative,
    "null_mu" = mu,
    "data_type" = data_type,
    "n_bits" = n_bits,
    "bit_probability" = bit_probability,
    "data" = tested_data
  )
  names(bf_out)[1] <- stat_name
  
  # Get final results
  final_bf <- tail(na.omit(bf_out$BF), n = 1)
  final_stat <- tail(na.omit(bf_out[[stat_name]]), n = 1)
  final_p <- tail(na.omit(bf_out$`p-value`), n = 1)
  final_delta <- tail(na.omit(bf_out$delta), n = 1)
  
  # Report final results
  if (is.na(final_bf)) {
    message("Final Bayes Factor could not be calculated")
  } else if (final_bf > 1) {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f (%s = %.3f; p = %.3f; \u03b4 = %.3f)",
      final_bf, stat_name, final_stat, final_p, final_delta
    ))
  } else if (final_bf == 1) {
    message(sprintf(
      "Final Bayes Factor: BF10 = 1 (no evidence; %s = %.3f; p = %.3f; \u03b4 = %.3f)",
      stat_name, final_stat, final_p, final_delta
    ))
  } else {
    message(sprintf(
      "Final Bayes Factor: BF10 = %.3f; BF01 = %.3f (%s = %.3f; p = %.3f; \u03b4 = %.3f)",
      final_bf, 1/final_bf, stat_name, final_stat, final_p, final_delta
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
#'     \item data: list with elements \code{x} and \code{y} (cleaned data vectors)
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
    .seqlast(nstart, n_data, .nstep(n_data))
  }
  
  # Initialize BF vector for full sample size
  bf <- rep(NA, n_data)
  
  # Pre-fill BF values for 1:(nstart-1) with 1
  if (nstart > 1) {
    bf[1:(nstart - 1)] <- 1
  }
  
  # Progress reporting
  message(sprintf("Correlation test (N = %d)", n_data))
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
    
    bf[n] <- bf_result  # Assign to position n, not sequential position
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
      "alpha" = 1/prior.r,
      "beta" = 1/prior.r
    ),
    "sample size" = n_data,
    "alternative" = alternative,
    "null_rho" = 0,
    "data_type" = "correlation",
    "data" = list(x = x, y = y)
  )
  
  # Get final BF
  final_bf <- tail(na.omit(bf), n = 1)  # Use na.omit here
  
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
#' \dontrun{
#' seqbf <- bfttest(rnorm(30, 0.5, 1), alternative = "two.sided", parallel = FALSE)
#'
#' bfrInformed <- bfRobustness(seqbf)
#' print(bfrInformed)
#'
#' bfrUninformed <- bfRobustness(seqbf, informed = FALSE)
#' plot(bfrInformed)
#' }
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
  
  # Handle binomial tests
  if (!is.null(x$`test type`) && x$`test type` == "binomial") {
    final_bf <- tail(na.omit(x$BF), n = 1)
    final_prob <- tail(na.omit(x$`probability of success`), n = 1)
    final_p <- tail(na.omit(x$`p-value`), n = 1)
    
    cat(sprintf("
  Sequential Bayesian Testing
  --------------------------------
  Test: Binomial proportion test
  Sample size: %d
  Final Bayes Factor: BF10 = %.3f; BF01 = %.3f
  Prior: %s(%g, %g)
  Alternative hypothesis: %s
  Probability of success: %.3f; p = %.3f
  \n",
                x$`sample size`,
                final_bf,
                1 / final_bf,
                x$prior$distribution,
                x$prior$location,
                x$prior$scale,
                x$alternative,
                final_prob,
                final_p
    ))
    
    return(invisible(x))
  }
  
  # Handle correlation tests
  if (!is.null(x$`test type`) && x$`test type` == "correlation") {
    final_bf <- tail(na.omit(x$BF), n = 1)
    final_r <- tail(na.omit(x$r), n = 1)
    final_p <- tail(na.omit(x$`p-value`), n = 1)
    
    cat(sprintf("
  Sequential Bayesian Testing
  --------------------------------
  Test: Pearson correlation test
  Sample size: %d
  Final Bayes Factor: BF10 = %.3f; BF01 = %.3f
  Prior: %s(%g, %g)
  Alternative hypothesis: %s
  Correlation: r = %.3f; p = %.3f
  \n",
                x$`sample size`,
                final_bf,
                1 / final_bf,
                x$prior$distribution,
                x$prior$location,
                x$prior$scale,
                x$alternative,
                final_r,
                final_p
    ))
    
    return(invisible(x))
  }
  
  # Handle t-tests (parametric and non-parametric)
  # Check if parametric field exists
  if (is.null(x$parametric)) {
    stop("Cannot determine test type: 'parametric' field is missing")
  }
  
  # Determine test name
  if (x$parametric) {
    test_name <- paste(.capitalize(x$`test type`), "t-test")
  } else {
    test_name <- if (x$`test type` == "independent") {
      "Mann-Whitney U test"
    } else {
      paste(.capitalize(x$`test type`), "Wilcoxon signed-rank test")
    }
  }
  
  # Get final values
  final_bf <- tail(na.omit(x$BF), n = 1)
  final_stat <- tail(na.omit(x[[1]]), n = 1)
  final_p <- tail(na.omit(x$`p-value`), n = 1)
  
  # Build delta string only if delta exists
  delta_string <- ""
  if (!is.null(x$delta) && !is.null(x$delta.lower) && !is.null(x$delta.upper)) {
    final_delta <- tail(na.omit(x$delta), n = 1)
    final_delta_lower <- tail(na.omit(x$delta.lower), n = 1)
    final_delta_upper <- tail(na.omit(x$delta.upper), n = 1)
    
    if (!is.na(final_delta)) {
      if (!x$parametric) {
        delta_string <- sprintf(
          "  Effect size: \u03b4 = %.3f, 95%% CI [%.3f, %.3f]",
          final_delta, final_delta_lower, final_delta_upper
        )
      } else {
        delta_string <- sprintf(
          "  Effect size: \u03b4 \u2248 %.3f (approximation from t-statistic)",
          final_delta
        )
      }
    }
  }
  
  # Print output
  cat(sprintf("
  Sequential Bayesian Testing
  --------------------------------
  Test: %s%s
  Sample size: %s
  Final Bayes Factor: BF10 = %.3f; BF01 = %.3f
  Prior: %s(%g, %g)
  Alternative hypothesis: %s
  %s: %.3f; p = %.3f%s
  \n",
              test_name,
              ifelse(x$parametric, " (parametric)", " (non-parametric)"),
              paste(x$`sample size`, collapse = ", "),
              final_bf,
              1 / final_bf,
              x$prior$distribution,
              x$prior$location,
              x$prior$scale,
              x$alternative,
              names(x)[1],
              final_stat,
              final_p,
              ifelse(delta_string != "", paste0("\n", delta_string), "")
  ))
  
  invisible(x)
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

# Adaptive number of samples for non-parametric tests
.adaptive_nsamples <- function(n, is_final = FALSE) {
  
  if (n < 500) {
    base_samples <- 1000
  } else if (n < 2000) {
    base_samples <- 750
  } else {
    base_samples <- 500
  }
  
  # Final calculation gets full precision
  if (is_final) {
    return(base_samples)
  }
  
  # Intermediate calculations: use reduced samples (but not too low)
  if (n < 500) {
    return(max(250, base_samples * 0.25))      # 250 samples
  } else if (n < 2000) {
    return(max(150, base_samples * 0.20))      # 150 samples
  } else {
    return(max(100, base_samples * 0.20))      # 100 samples (not 75)
  }
}

# Try running bfttest with specific n
.try_bft_calculation <- function(n, ...) {
  nstart <- n - 1
  tryCatch({
    sink(tempfile())  # Redirect output to a temp file
    # Explicitly call the PACKAGE version to avoid recursion
    result <- suppressMessages(suppressWarnings(
      changeofevidence::bfttest(..., nstart = nstart, exact = TRUE)
    ))
    sink()  # Restore normal output
    
    # Check if we got valid BF values
    !is.na(result$BF[length(result$BF) - 1])
  }, error = function(e) {
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
.determine_min_n_one_sample <- function(x, alternative = "two.sided", prior.loc = 0, prior.r = 0.1, mu = 0) {
  n <- 3
  nstart <- n - 1
  while (n <= length(x)) {
    if (var(x[1:nstart]) > 0) {
      # Try running bfttest - NOW PASSING mu!
      if (.try_bft_calculation(n,
                               x = x[1:n],
                               alternative = alternative,
                               prior.loc = prior.loc,
                               prior.r = prior.r,
                               mu = mu)) {
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

# Helper function for computing BF from posterior samples (non-parametric)
.get_bf_nonparametric <- function(delta_samples, alternative, prior_scale) {
  
  if (alternative == "two.sided") {
    # Two-sided test using Savage-Dickey ratio
    bf10 <- computeBayesFactorOneZero(
      posteriorSamples = delta_samples,
      priorParameter = prior_scale,
      oneSided = FALSE,
      whichTest = "Wilcoxon"
    )
  } else if (alternative == "greater") {
    # One-sided test: H1: delta > 0
    bf10 <- computeBayesFactorOneZero(
      posteriorSamples = delta_samples,
      priorParameter = prior_scale,
      oneSided = "right",
      whichTest = "Wilcoxon"
    )
  } else {  # "less"
    # One-sided test: H1: delta < 0
    bf10 <- computeBayesFactorOneZero(
      posteriorSamples = delta_samples,
      priorParameter = prior_scale,
      oneSided = "left",
      whichTest = "Wilcoxon"
    )
  }
  
  return(bf10)
}

# Capitalize first letter
.capitalize <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
