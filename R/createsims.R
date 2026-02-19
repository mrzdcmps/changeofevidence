# Helper Functions for Simulation Parameter Inference and Validation

# Infer n_bits from mu and p
.infer_n_bits <- function(mu, p = 0.5) {
  if (mu <= 0 && p == 0.5) {
    stop("Cannot infer n_bits from mu=0. Please specify n_bits explicitly or set data_type='continuous'")
  }

  n_bits_raw <- abs(mu / p)
  n_bits_rounded <- round(n_bits_raw)

  if (n_bits_rounded < 2) {
    stop(sprintf("Inferred n_bits = %d is too small (mu=%.2f, p=%.2f). Consider data_type='continuous'",
                 n_bits_rounded, mu, p))
  }

  if (n_bits_rounded > 1000) {
    warning(sprintf("Inferred n_bits = %d is very large (mu=%.2f, p=%.2f). Consider specifying manually.",
                    n_bits_rounded, mu, p))
  }

  if (abs(n_bits_raw - n_bits_rounded) > 0.1) {
    warning(sprintf(
      "mu/p = %.2f is not an integer. Rounding to n_bits=%d. Expected mean will be %.2f",
      n_bits_raw, n_bits_rounded, n_bits_rounded * p
    ))
  }

  return(n_bits_rounded)
}

# Validate t-test simulation parameters
.validate_t_params <- function(N, mu, n_bits, p, data_type) {
  if (N < 3) {
    stop("N must be at least 3 for meaningful simulations")
  }

  if (data_type == "summed_bits") {
    expected_mean <- n_bits * p

    if (abs(mu - expected_mean) > 0.01) {
      warning(sprintf(
        paste("mu (%.2f) does not match expected mean from n_bits*p (%.2f).",
              "Simulations will have expected mean %.2f"),
        mu, expected_mean, expected_mean
      ))
    }

    if (n_bits < 1) {
      stop("n_bits must be at least 1 for summed_bits data type")
    }
  }

  if (p <= 0 || p >= 1) {
    stop("p (bit_probability) must be between 0 and 1")
  }

  invisible(TRUE)
}

# Extract nstart from seqbf object
.extract_nstart <- function(seqbf_obj) {
  bf_vec <- seqbf_obj$BF
  first_calculated <- which(!is.na(bf_vec) & bf_vec != 1)[1]
  if (is.na(first_calculated)) return(5)  # Default fallback
  return(first_calculated)
}

# Generate binary data from Quantis device
.generate_binary_quantis <- function(n, min = 0, max = 1) {
  system(sprintf('EasyQuantis -u 0 -n %d --min %d --max %d -i "tmp_quantis.txt"',
                 n, min, max), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (!file.exists("tmp_quantis.txt")) {
    stop("Could not create random data from QRNG. Is the device connected?")
  }
  data <- read.table("tmp_quantis.txt", header = FALSE)$V1
  unlink("tmp_quantis.txt")
  return(data)
}

# Generate float data from Quantis device
.generate_float_quantis <- function(n, min = 0, max = 1) {
  system(sprintf('EasyQuantis -u 0 -n %d --min %.6f --max %.6f -f "tmp_quantis.txt"',
                 n, min, max), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (!file.exists("tmp_quantis.txt")) {
    stop("Could not create random data from QRNG. Is the device connected?")
  }
  data <- read.table("tmp_quantis.txt", header = FALSE)$V1
  unlink("tmp_quantis.txt")
  return(data)
}

# NULL coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Find the simplest rational approximation p_num/p_den for a float p.
# Tries denominators 1..max_den and returns the first where |round(p*d)/d - p| < tol.
# Works reliably for any p that was itself specified as a simple fraction (0.5, 0.6, 0.3, etc.).
.float_to_rational <- function(p, max_den = 10000, tol = .Machine$double.eps^0.5) {
  for (d in 1:max_den) {
    n <- round(p * d)
    if (abs(n / d - p) < tol) {
      g <- function(a, b) if (b == 0L) a else Recall(b, a %% b)  # GCD via Euclidean algorithm
      gc <- g(as.integer(n), as.integer(d))
      return(list(num = as.integer(n / gc), den = as.integer(d / gc)))
    }
  }
  stop(sprintf(
    "Cannot represent p=%.15g as a fraction with denominator <= %d. Supply p_num and p_den explicitly.",
    p, max_den
  ))
}

# Helper to find quantum data file automatically
.find_quantum_data <- function() {
  # Priority order:
  # 1. Previously downloaded full dataset (from options)
  # 2. Package extdata (full dataset)
  # 3. Sample data (package extdata)
  # 4. NULL (will trigger fallback)

  # Check saved path from previous download
  saved_path <- getOption("changeofevidence.quantum_data_path")
  if (!is.null(saved_path) && file.exists(saved_path)) {
    return(list(path = saved_path, type = "full", size = length(readRDS(saved_path))))
  }

  # Check package extdata for full dataset
  full_path <- system.file("extdata", "quantis_bits.rds", package = "changeofevidence")
  if (file.exists(full_path) && file.size(full_path) > 1e6) {  # > 1MB means full dataset
    return(list(path = full_path, type = "full", size = length(readRDS(full_path))))
  }

  # Check for sample data
  sample_path <- system.file("extdata", "quantis_bits_sample.rds", package = "changeofevidence")
  if (file.exists(sample_path)) {
    return(list(path = sample_path, type = "sample", size = length(readRDS(sample_path))))
  }

  # No data found
  return(NULL)
}

# Helper to check if quantum data is sufficient and handle fallback
.check_quantum_data <- function(required_bits, data_info, allow_fallback = TRUE) {
  if (is.null(data_info)) {
    if (allow_fallback) {
      message(sprintf("No quantum data found (%s bits required).",
                      format(required_bits, big.mark = ",")))
      return(.prompt_download_or_pseudo())
    } else {
      stop("No quantum data found. Run download_quantum_data() or use method='pseudo'")
    }
  }

  if (data_info$size < required_bits) {
    if (data_info$type == "sample") {
      if (allow_fallback) {
        message(sprintf(
          "Sample quantum data insufficient: %s bits needed, %s available.",
          format(required_bits, big.mark = ","),
          format(data_info$size, big.mark = ",")
        ))
        return(.prompt_download_or_pseudo())
      } else {
        stop(sprintf(
          "Sample data insufficient (%s bits needed, %s available).\nDownload full dataset: download_quantum_data()",
          format(required_bits, big.mark = ","),
          format(data_info$size, big.mark = ",")
        ))
      }
    } else {
      stop(sprintf(
        "Insufficient data: need %s bits, have %s bits.\nReduce n.sims or N.",
        format(required_bits, big.mark = ","),
        format(data_info$size, big.mark = ",")
      ))
    }
  }

  # Sufficient data
  return(list(method = "files", path = data_info$path))
}

# Prompt user to download full dataset or fall back to pseudo
.prompt_download_or_pseudo <- function() {
  if (!interactive()) {
    stop(paste(
      "Insufficient quantum data for the requested simulations.",
      "In a non-interactive session, either run download_quantum_data() beforehand",
      "or use method='pseudo'.",
      sep = " "
    ))
  }

  choice <- utils::menu(
    choices = c(
      "Download full quantum dataset (~125 MB) via download_quantum_data()",
      "Use pseudo-random generation instead"
    ),
    title = "Quantum data insufficient. How would you like to proceed?"
  )

  if (choice == 1) {
    download_quantum_data()
    data_info <- .find_quantum_data()
    if (is.null(data_info)) {
      stop("Download appears to have failed. Please try again or use method='pseudo'.")
    }
    return(list(method = "files", path = data_info$path))
  } else {
    message("Using pseudo-random generation.")
    return(list(method = "pseudo", path = NULL))
  }
}

#' Generate Simulations for Change of Evidence Analysis
#'
#' Convenience wrapper that automatically generates simulations based on a seqbf object,
#' or provides legacy interface (deprecated).
#'
#' @param x A seqbf object from bfttest(), bfbinom(), or bfcor(), OR trials (deprecated)
#' @param n.sims Number of simulations (default: 1000)
#' @param N Override sample size (extracted from seqbf if not provided)
#' @param ... Additional parameters passed to type-specific functions
#' @param mu Override null hypothesis mean (for t-tests)
#' @param n_bits Override n_bits (for t-tests with summed bits)
#' @param p Override bit probability (for binomial or t-tests)
#' @param rho Override null correlation (for correlation tests)
#' @param data_type Override data type
#' @param int_min,int_max Override integer range
#' @param mean.scores DEPRECATED - use n_bits instead
#' @param trials DEPRECATED - use N instead
#' @param method Random generation method: "pseudo", "files", or "quantis"
#' @param filespath Path to quantum random data file
#' @param parallel Use parallel processing (default: TRUE)
#' @param nstart Minimum observations before calculating first BF
#' @param alternative Alternative hypothesis
#' @param prior.loc Prior location parameter
#' @param prior.r Prior scale parameter
#' @param use.files DEPRECATED - use method="files"
#' @param use.quantis DEPRECATED - use method="quantis"
#'
#' @return Dataframe with simulation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Automatic simulation from test result
#' result <- bfttest(rnorm(50, 0.5), mu = 0)
#' sims <- simcreate(result, n.sims = 1000)
#'
#' # Manual binomial simulation
#' sims <- simcreate_bin(N = 100, p = 0.5, n.sims = 1000)
#' }
simcreate <- function(x = NULL,
                     n.sims = 1000,
                     N = NULL,
                     ...,
                     # Override parameters
                     mu = NULL,
                     n_bits = NULL,
                     p = NULL,
                     rho = NULL,
                     data_type = NULL,
                     int_min = NULL,
                     int_max = NULL,
                     # Legacy parameters (deprecated)
                     trials = NULL,
                     mean.scores = NULL,
                     # Common parameters
                     method = NULL,
                     filespath = NULL,
                     parallel = TRUE,
                     nstart = NULL,
                     exact = TRUE,
                     alternative = NULL,
                     prior.loc = NULL,
                     prior.r = NULL,
                     use.files = NULL,
                     use.quantis = NULL) {

  # Handle legacy use.files/use.quantis
  if (!is.null(use.files) && use.files == TRUE) {
    message("The argument `use.files` is deprecated. Please use `method = 'files'` instead.")
    method <- "files"
  }
  if (!is.null(use.quantis) && use.quantis == TRUE) {
    message("The argument `use.quantis` is deprecated. Please use `method = 'quantis'` instead.")
    method <- "quantis"
  }

  # Handle new interface: x is seqbf object
  if (inherits(x, "seqbf")) {
    test_type <- x$`test type`

    # Warn if exact levels may not match
    if (!exact && length(na.omit(x$BF)) > 100) {
      warning(paste(
        "The seqbf object has more than 100 calculated BF values (exact=TRUE data),",
        "but simulations will use exact=FALSE (~100 steps).",
        "This mismatch may bias CoE statistics. Consider setting exact=TRUE for simulations",
        "or recomputing the seqbf object with exact=FALSE."
      ))
    }

    # Extract sample size
    N_to_use <- N %||% {
      sample_size <- x$`sample size`
      if (length(sample_size) > 1) sum(sample_size) else sample_size  # Handle independent samples
    }

    # Extract nstart
    nstart_to_use <- nstart %||% .extract_nstart(x)

    # Dispatch based on test type
    if (test_type == "binomial") {
      return(simcreate_bin(
        N = N_to_use,
        n.sims = n.sims,
        p = p %||% x$null_p %||% 0.5,
        method = method %||% "pseudo",
        filespath = filespath,
        parallel = parallel,
        nstart = nstart_to_use,
        exact = exact,
        prior.r = prior.r %||% x$prior$scale,
        alternative = alternative %||% x$alternative
      ))

    } else if (test_type %in% c("one-sample", "paired", "independent")) {
      # T-test
      mu_to_use <- mu %||% x$null_mu %||% 0
      data_type_to_use <- data_type %||% x$data_type %||% "unknown"
      n_bits_to_use <- n_bits %||% x$n_bits
      p_to_use <- p %||% x$bit_probability %||% 0.5

      return(simcreate_t(
        N = N_to_use,
        n.sims = n.sims,
        mu = mu_to_use,
        data_type = data_type_to_use,
        n_bits = n_bits_to_use,
        p = p_to_use,
        int_min = int_min,
        int_max = int_max,
        method = method %||% "pseudo",
        filespath = filespath,
        parallel = parallel,
        nstart = nstart_to_use,
        exact = exact,
        alternative = alternative %||% x$alternative,
        prior.loc = prior.loc %||% x$prior$location,
        prior.r = prior.r %||% x$prior$scale
      ))

    } else if (test_type == "correlation") {
      # For correlation, prior.r is stored as 1/alpha (Beta distribution)
      prior_r_cor <- prior.r %||% {
        if (!is.null(x$prior$alpha)) 1/x$prior$alpha else 0.353
      }

      return(simcreate_cor(
        N = N_to_use,
        n.sims = n.sims,
        rho = rho %||% x$null_rho %||% 0,
        method = method %||% "pseudo",
        parallel = parallel,
        nstart = nstart_to_use,
        exact = exact,
        prior.r = prior_r_cor
      ))

    } else {
      stop(sprintf("Unknown test type: %s", test_type))
    }
  }

  # Handle legacy interface (DEPRECATED)
  if (!is.null(trials) || (!is.null(x) && is.numeric(x) && !inherits(x, "seqbf"))) {
    .Deprecated(
      new = "simcreate_bin() / simcreate_t() / simcreate_cor() or simcreate(seqbf_object)",
      old = "simcreate(trials, mean.scores, ...)",
      package = "changeofevidence"
    )

    # Map old parameters
    if (is.numeric(x) && is.null(trials)) {
      trials <- x
    }

    if (is.null(mean.scores)) {
      # Binomial case
      warning("Using deprecated interface. Please use simcreate_bin() instead.")
      return(simcreate_bin(
        N = trials,
        n.sims = n.sims,
        p = p %||% 0.5,
        method = method %||% "pseudo",
        filespath = filespath,
        parallel = parallel,
        nstart = nstart %||% 5,
        prior.r = prior.r %||% 0.1,
        alternative = alternative %||% "two.sided"
      ))
    } else {
      # T-test case with summed bits
      # Need to calculate N from trials and mean.scores
      n_bits_legacy <- mean.scores
      p_legacy <- p %||% 0.5
      N_legacy <- trials %/% n_bits_legacy

      warning(sprintf(
        paste("Using deprecated interface. Calculated N=%d from trials=%d and mean.scores=%d.",
              "Please use simcreate_t() with N parameter instead."),
        N_legacy, trials, n_bits_legacy
      ))

      return(simcreate_t(
        N = N_legacy,
        n.sims = n.sims,
        mu = mu %||% mean.scores,
        data_type = "summed_bits",
        n_bits = n_bits_legacy,
        p = p_legacy,
        method = method %||% "pseudo",
        filespath = filespath,
        parallel = parallel,
        nstart = nstart %||% 5,
        alternative = alternative %||% "two.sided",
        prior.loc = prior.loc %||% 0,
        prior.r = prior.r %||% 0.353
      ))
    }
  }

  # If we get here, missing required arguments
  stop(paste(
    "simcreate() requires a seqbf object as first argument.",
    "For manual simulation generation, use simcreate_bin(), simcreate_t(), or simcreate_cor().",
    "See ?simcreate for examples."
  ))
}

# Helper function to suppress all output (messages, warnings, print, cat, etc.)
.quiet <- function(x) {
  # Suppress all messages, warnings, and output
  suppressMessages(suppressWarnings({
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }))
}


#' Generate Binomial Test Simulations
#'
#' Generates Monte Carlo simulations for binomial tests with specified sample size
#' and null hypothesis probability.
#'
#' @param N Sample size (number of binary observations)
#' @param n.sims Number of simulations (default: 1000)
#' @param p Null hypothesis probability (default: 0.5)
#' @param method Random generation method: "pseudo", "files", or "quantis"
#' @param filespath Path to quantum random data file (for method="files")
#' @param parallel Use parallel processing (default: TRUE)
#' @param nstart Minimum observations before calculating first BF (default: 5)
#' @param prior.r Prior scale parameter (default: 0.1)
#' @param alternative Alternative hypothesis: "two.sided", "less", or "greater"
#'
#' @return Dataframe with columns: simid, index, raw, rw, density.rw, bf, density.bf
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate binomial simulations
#' sims <- simcreate_bin(N = 100, p = 0.5, n.sims = 100)
#' }
simcreate_bin <- function(N,
                         n.sims = 1000,
                         p = 0.5,
                         method = c("pseudo", "files", "quantis"),
                         filespath = NULL,
                         parallel = TRUE,
                         nstart = 5,
                         exact = TRUE,
                         prior.r = 0.1,
                         alternative = c("two.sided", "less", "greater")) {

  method <- match.arg(method)
  alternative <- match.arg(alternative)

  if (N < 3) stop("N must be at least 3")
  if (p <= 0 || p >= 1) stop("p must be between 0 and 1")

  # Infer exact rational p_num/p_den from the float p, then set up rejection sampling.
  rat <- .float_to_rational(p)
  p_num <- rat$num
  p_den <- rat$den
  valid_range <- (65536L %/% p_den) * p_den   # largest multiple of p_den <= 65536
  chunks_per_sim <- ceiling(N * 1.01)          # 1% buffer absorbs rare rejections

  # Auto-detect quantum data and handle fallback
  quantum_data <- NULL
  quantum_filepath <- NULL
  if (method == "files") {
    if (is.null(filespath)) {
      # Automatic detection with smart fallback
      required_bits <- chunks_per_sim * n.sims * 16L
      data_info <- .find_quantum_data()
      fallback_result <- .check_quantum_data(required_bits, data_info, allow_fallback = TRUE)

      method <- fallback_result$method
      if (method == "files") {
        quantum_filepath <- fallback_result$path
        quantum_data <- if (!parallel) readRDS(fallback_result$path) else NULL
      }

    } else {
      # User specified path - strict checking (no fallback)
      if (!file.exists(filespath)) {
        stop(sprintf("File not found: %s", filespath))
      }
      quantum_filepath <- filespath
      quantum_data <- if (!parallel) readRDS(filespath) else NULL  # Only load for non-parallel
      required_bits <- chunks_per_sim * n.sims * 16L
      if (!parallel) {
        if (length(quantum_data) < required_bits) {
          stop(sprintf(
            "Insufficient data in %s: need %s bits, have %s bits.",
            basename(filespath),
            format(required_bits, big.mark = ","),
            format(length(quantum_data), big.mark = ",")
          ))
        }
      } else {
        # For parallel, check size without loading all data
        actual_size <- length(readRDS(filespath))
        if (actual_size < required_bits) {
          stop(sprintf(
            "Insufficient data in %s: need %s bits, have %s bits.",
            basename(filespath),
            format(required_bits, big.mark = ","),
            format(actual_size, big.mark = ",")
          ))
        }
      }
    }
  }

  # For parallel + files: pre-load only the required bits (avoids per-worker full-file reads)
  quantum_bits <- NULL
  if (method == "files" && parallel && !is.null(quantum_filepath)) {
    required_bits <- chunks_per_sim * n.sims * 16L
    all_bits <- readRDS(quantum_filepath)
    quantum_bits <- all_bits[1:required_bits]
    rm(all_bits)
  }

  # For quantis: pre-generate all integers from device in one sequential call,
  # then parallel workers access their slice by index — no per-worker device access needed.
  quantum_int <- NULL
  if (method == "quantis") {
    total_ints <- chunks_per_sim * n.sims
    message(sprintf("Generating %d integers from Quantis device...", total_ints))
    quantum_int <- .generate_binary_quantis(total_ints, min = 0L, max = 65535L)
  }

  # Define single simulation function
  generate_simulation <- function(i) {
    # Generate binary data with exact probability p_num/p_den via rejection sampling.
    # For method="files"/"quantis": draw 16-bit integers uniform in [0, 65535], reject
    # those >= valid_range, then output 1 iff accepted_int %% p_den < p_num.
    if (method == "pseudo") {
      sim_data <- rbinom(N, 1, p)
    } else if (method == "files") {
      start_idx <- (i - 1) * chunks_per_sim * 16L + 1
      end_idx <- i * chunks_per_sim * 16L
      if (!is.null(quantum_bits)) {
        raw_quantum <- quantum_bits[start_idx:end_idx]
      } else {
        raw_quantum <- quantum_data[start_idx:end_idx]
      }
      bit_groups <- matrix(raw_quantum, nrow = chunks_per_sim, ncol = 16L, byrow = TRUE)
      integers <- as.numeric(bit_groups %*% 2^(0:15))
      accepted <- integers[integers < valid_range]
      sim_data <- as.integer(accepted[1:N] %% p_den < p_num)
    } else if (method == "quantis") {
      start_idx <- (i - 1) * chunks_per_sim + 1
      end_idx <- i * chunks_per_sim
      raw_ints <- quantum_int[start_idx:end_idx]
      accepted <- raw_ints[raw_ints < valid_range]
      sim_data <- as.integer(accepted[1:N] %% p_den < p_num)
    }

    # Create dataframe
    sim <- data.frame(
      simid = i,
      index = 1:N,
      V1 = sim_data
    )

    # Calculate random walk (convert 0 to -1 for cumsum)
    qbitmin1 <- ifelse(sim$V1 == 0, -1, 1)
    sim$rw <- cumsum(qbitmin1)

    # Calculate BFs
    sim$bf <- .quiet(changeofevidence::bfbinom(
      sim$V1,
      p = p,
      prior.r = prior.r,
      nstart = nstart,
      exact = exact,
      alternative = alternative
    )$BF)

    return(sim)
  }

  # Run simulations
  if (parallel) {
    cores <- parallel::detectCores() - 1
    message(sprintf("Generating %d binomial simulations (N=%d) using %d cores...", n.sims, N, cores))
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("generate_simulation", "quantum_bits", "quantum_data", "quantum_int", ".quiet",
                                   "N", "p", "p_num", "p_den", "valid_range", "chunks_per_sim",
                                   "method", "nstart", "exact", "prior.r", "alternative"),
                           envir = environment())
    parallel::clusterEvalQ(cl, library(changeofevidence))
    sims_list <- pbapply::pblapply(1:n.sims, generate_simulation, cl = cl)
    parallel::stopCluster(cl)
  } else {
    message(sprintf("Generating %d binomial simulations (N=%d)...", n.sims, N))
    sims_list <- pbapply::pblapply(1:n.sims, generate_simulation)
  }

  # Combine and process
  sims_df <- do.call(rbind, sims_list)
  sims_df$raw <- sims_df$V1
  sims_df$V1 <- NULL

  # FFT processing (suppress all messages and output)
  sims_df <- .quiet(simredo(sims_df, N, rw = TRUE))

  return(sims_df)
}
#' Generate T-Test Simulations
#'
#' Generates Monte Carlo simulations for t-tests with support for different data types:
#' summed bits, continuous normal data, or random integers.
#'
#' @param N Sample size (number of observations after summing/grouping)
#' @param n.sims Number of simulations (default: 1000)
#' @param mu Null hypothesis mean (default: 0 for continuous/integer, n_bits*p for summed_bits)
#' @param data_type Type of data: "summed_bits", "continuous", "integer", or "unknown"
#' @param n_bits Number of bits to sum per observation (for summed_bits)
#' @param p Probability for each bit (default: 0.5)
#' @param int_min Minimum integer value (for integer data_type)
#' @param int_max Maximum integer value (for integer data_type)
#' @param method Random generation method: "pseudo", "files", or "quantis"
#' @param filespath Path to quantum random data file (for method="files")
#' @param parallel Use parallel processing (default: TRUE)
#' @param nstart Minimum observations before calculating first BF (default: 5)
#' @param alternative Alternative hypothesis: "two.sided", "less", or "greater"
#' @param prior.loc Prior location parameter (default: 0)
#' @param prior.r Prior scale parameter (default: 0.353)
#'
#' @return Dataframe with columns: simid, index, raw, rw, density.rw, bf, density.bf
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate t-test simulations with summed bits
#' sims <- simcreate_t(N = 100, mu = 10, data_type = "summed_bits",
#'                    n_bits = 20, p = 0.5, n.sims = 100)
#'
#' # Generate t-test simulations with continuous data
#' sims <- simcreate_t(N = 100, mu = 0, data_type = "continuous", n.sims = 100)
#' }
simcreate_t <- function(N,
                       n.sims = 1000,
                       mu = NULL,
                       data_type = c("summed_bits", "continuous", "integer", "unknown"),
                       n_bits = NULL,
                       p = 0.5,
                       int_min = NULL,
                       int_max = NULL,
                       method = c("pseudo", "files", "quantis"),
                       filespath = NULL,
                       parallel = TRUE,
                       nstart = 5,
                       exact = TRUE,
                       alternative = c("two.sided", "less", "greater"),
                       prior.loc = 0,
                       prior.r = 0.353) {

  method <- match.arg(method)
  data_type <- match.arg(data_type)

  alternative <- match.arg(alternative)

  # Infer data_type if unknown
  if (data_type == "unknown") {
    if (!is.null(n_bits)) {
      data_type <- "summed_bits"
      message("data_type set to 'summed_bits' based on n_bits parameter")
    } else if (!is.null(int_min) && !is.null(int_max)) {
      data_type <- "integer"
      message("data_type set to 'integer' based on int_min/int_max parameters")
    } else {
      data_type <- "continuous"
      warning(paste("data_type='unknown' and no additional parameters specified.",
                   "Defaulting to 'continuous' (normal distribution).",
                   "Specify n_bits for summed_bits or int_min/int_max for integers."))
    }
  }

  # Handle summed_bits: infer n_bits and/or mu if needed
  if (data_type == "summed_bits") {
    if (is.null(n_bits)) {
      n_bits <- .infer_n_bits(mu %||% 0, p)
      warning(sprintf("n_bits inferred as %d from mu=%.2f and p=%.2f", n_bits, mu %||% 0, p))
    }
    if (is.null(mu)) {
      mu <- n_bits * p
      message(sprintf("mu not specified: defaulting to n_bits * p = %g for summed_bits data", mu))
    }
    .validate_t_params(N, mu, n_bits, p, data_type)
    trials <- N * n_bits  # Total raw bits needed
  } else {
    if (is.null(mu)) mu <- 0
    trials <- N  # For continuous/integer, trials = N
  }

  # Handle integer data_type
  if (data_type == "integer") {
    if (is.null(int_min) || is.null(int_max)) {
      stop("data_type='integer' requires int_min and int_max parameters")
    }
    if (int_min >= int_max) {
      stop("int_min must be less than int_max")
    }
  }

  # File method only supported for binary/summed_bits
  if (method == "files" && data_type != "summed_bits") {
    stop(sprintf("method='files' only supported for data_type='summed_bits'. You specified '%s'", data_type))
  }

  # For p != 0.5, use 2-byte integers (uniform in [0, 65535]) with rejection sampling to
  # achieve exact Bernoulli(p) where p = mu/n_bits (a rational number with denominator n_bits).
  # Rejection rate = (65536 mod n_bits) / 65536 <= (n_bits-1)/65536  (< 0.002 for n_bits <= 100).
  BITS_PER_OUTPUT <- if (method == "files" && data_type == "summed_bits" && !isTRUE(all.equal(p, 0.5))) 16L else 1L
  rat         <- if (BITS_PER_OUTPUT == 16L) .float_to_rational(p) else NULL
  p_num       <- if (BITS_PER_OUTPUT == 16L) rat$num else NULL
  p_den       <- if (BITS_PER_OUTPUT == 16L) rat$den else NULL
  valid_range <- if (BITS_PER_OUTPUT == 16L) (65536L %/% p_den) * p_den else NULL     # largest multiple of p_den <= 65536
  # Draw slightly more chunks than trials to absorb rare rejections (1% buffer >> expected rejections)
  chunks_per_sim <- if (BITS_PER_OUTPUT == 16L) ceiling(trials * 1.01) else trials

  # Auto-detect quantum data and handle fallback
  quantum_data <- NULL
  quantum_filepath <- NULL
  if (method == "files") {
    if (is.null(filespath)) {
      # Automatic detection with smart fallback
      required_bits <- chunks_per_sim * n.sims * BITS_PER_OUTPUT
      data_info <- .find_quantum_data()
      fallback_result <- .check_quantum_data(required_bits, data_info, allow_fallback = TRUE)

      method <- fallback_result$method
      if (method == "files") {
        quantum_filepath <- fallback_result$path
        quantum_data <- if (!parallel) readRDS(fallback_result$path) else NULL
      }

    } else {
      # User specified path - strict checking (no fallback)
      if (!file.exists(filespath)) {
        stop(sprintf("File not found: %s", filespath))
      }
      quantum_filepath <- filespath
      quantum_data <- if (!parallel) readRDS(filespath) else NULL  # Only load for non-parallel
      required_bits <- chunks_per_sim * n.sims * BITS_PER_OUTPUT
      if (!parallel) {
        if (length(quantum_data) < required_bits) {
          stop(sprintf(
            "Insufficient data in %s: need %s bits, have %s bits.",
            basename(filespath),
            format(required_bits, big.mark = ","),
            format(length(quantum_data), big.mark = ",")
          ))
        }
      } else {
        actual_size <- length(readRDS(filespath))
        if (actual_size < required_bits) {
          stop(sprintf(
            "Insufficient data in %s: need %s bits, have %s bits.",
            basename(filespath),
            format(required_bits, big.mark = ","),
            format(actual_size, big.mark = ",")
          ))
        }
      }
    }
  }

  # For parallel + files: pre-load only the required bits (avoids per-worker full-file reads)
  quantum_bits <- NULL
  if (method == "files" && parallel && !is.null(quantum_filepath)) {
    required_bits <- chunks_per_sim * n.sims * BITS_PER_OUTPUT
    all_bits <- readRDS(quantum_filepath)
    quantum_bits <- all_bits[1:required_bits]
    rm(all_bits)
  }

  # For quantis: pre-generate all data from device in one sequential call,
  # then parallel workers access their slice by index — no per-worker device access needed.
  quantum_int <- NULL
  if (method == "quantis") {
    if (data_type == "summed_bits") {
      total_ints <- chunks_per_sim * n.sims
      max_val <- if (BITS_PER_OUTPUT == 16L) 65535L else 1L
      message(sprintf("Generating %d integers from Quantis device...", total_ints))
      quantum_int <- .generate_binary_quantis(total_ints, min = 0L, max = max_val)
    } else if (data_type == "continuous") {
      total_floats <- N * n.sims
      message(sprintf("Generating %d floats from Quantis device...", total_floats))
      quantum_int <- .generate_float_quantis(total_floats, min = 0, max = 1)
    } else if (data_type == "integer") {
      total_ints <- N * n.sims
      message(sprintf("Generating %d integers from Quantis device...", total_ints))
      quantum_int <- .generate_binary_quantis(total_ints, min = int_min, max = int_max)
    }
  }

  # Define single simulation function
  generate_simulation <- function(i) {
    # Generate data based on type
    if (data_type == "summed_bits") {
      # Generate binary data
      if (method == "pseudo") {
        raw_bits <- rbinom(trials, 1, p)
      } else if (method == "files") {
        start_idx <- (i - 1) * chunks_per_sim * BITS_PER_OUTPUT + 1
        end_idx <- i * chunks_per_sim * BITS_PER_OUTPUT
        if (!is.null(quantum_bits)) {
          raw_quantum <- quantum_bits[start_idx:end_idx]  # Pre-sliced subset (parallel)
        } else {
          raw_quantum <- quantum_data[start_idx:end_idx]  # Full data (non-parallel)
        }
        if (BITS_PER_OUTPUT == 1L) {
          raw_bits <- raw_quantum
        } else {
          # Convert 16-bit chunks to integers uniform in [0, 65535],
          # then apply rejection sampling for exact Bernoulli(p_num/p_den).
          n_chunks <- length(raw_quantum) %/% 16L
          bit_groups <- matrix(raw_quantum[1:(n_chunks * 16L)], nrow = n_chunks, ncol = 16L, byrow = TRUE)
          integers <- as.numeric(bit_groups %*% 2^(0:15))
          # Reject integers >= valid_range so the remainder is exactly uniform in [0, p_den)
          accepted <- integers[integers < valid_range]
          raw_bits <- as.integer(accepted[1:trials] %% p_den < p_num)
        }
      } else if (method == "quantis") {
        start_idx <- (i - 1) * chunks_per_sim + 1
        end_idx <- i * chunks_per_sim
        raw_ints <- quantum_int[start_idx:end_idx]
        if (BITS_PER_OUTPUT == 1L) {
          raw_bits <- raw_ints  # pre-loaded 0/1 bits (p = 0.5)
        } else {
          accepted <- raw_ints[raw_ints < valid_range]
          raw_bits <- as.integer(accepted[1:trials] %% p_den < p_num)
        }
      }

      # Sum in groups
      groups <- rep(1:N, each = n_bits)
      sim_data <- as.numeric(tapply(raw_bits, groups, sum))

    } else if (data_type == "continuous") {
      # Generate continuous normal data
      if (method == "pseudo") {
        sim_data <- rnorm(N, mean = mu, sd = 1)
      } else if (method == "quantis") {
        uniform_data <- quantum_int[(i - 1) * N + 1:(i * N)]
        sim_data <- qnorm(uniform_data, mean = mu, sd = 1)
      }

    } else if (data_type == "integer") {
      # Generate random integers
      if (method == "pseudo") {
        sim_data <- sample(int_min:int_max, N, replace = TRUE)
      } else if (method == "quantis") {
        sim_data <- quantum_int[(i - 1) * N + 1:(i * N)]
      }
    }

    # Create dataframe
    sim <- data.frame(
      simid = i,
      index = 1:N,
      raw = sim_data
    )

    # Calculate random walk (cumulative sum centered at mu)
    sim$rw <- cumsum(sim$raw - mu)

    # Calculate BFs
    sim$bf <- .quiet(changeofevidence::bfttest(
      sim$raw,
      mu = mu,
      prior.loc = prior.loc,
      prior.r = prior.r,
      alternative = alternative,
      nstart = nstart,
      exact = exact,
      parallel = FALSE
    )$BF)

    return(sim)
  }

  # Run simulations
  if (parallel) {
    cores <- parallel::detectCores() - 1
    message(sprintf("Generating %d t-test simulations (N=%d, data_type=%s) using %d cores...",
                   n.sims, N, data_type, cores))
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("generate_simulation", "quantum_bits", "quantum_data", "quantum_int", ".quiet",
                                   "N", "mu", "data_type", "n_bits", "p", "int_min", "int_max",
                                   "trials", "method", "nstart", "exact", "prior.loc", "prior.r", "alternative",
                                   "BITS_PER_OUTPUT", "chunks_per_sim",
                                   "p_num", "p_den", "valid_range"),
                           envir = environment())
    parallel::clusterEvalQ(cl, library(changeofevidence))
    sims_list <- pbapply::pblapply(1:n.sims, generate_simulation, cl = cl)
    parallel::stopCluster(cl)
  } else {
    message(sprintf("Generating %d t-test simulations (N=%d, data_type=%s)...",
                   n.sims, N, data_type))
    sims_list <- pbapply::pblapply(1:n.sims, generate_simulation)
  }

  # Combine and process
  sims_df <- do.call(rbind, sims_list)

  # FFT processing (suppress all messages and output)
  sims_df <- .quiet(simredo(sims_df, N, rw = TRUE))

  return(sims_df)
}


#' Generate Correlation Test Simulations
#'
#' Generates Monte Carlo simulations for correlation tests with specified sample size
#' and null hypothesis correlation.
#'
#' @param N Sample size (number of paired observations)
#' @param n.sims Number of simulations (default: 1000)
#' @param rho Null hypothesis correlation (default: 0)
#' @param method Random generation method: "pseudo" or "quantis" (files not supported)
#' @param parallel Use parallel processing (default: TRUE)
#' @param nstart Minimum observations before calculating first BF (default: 5)
#' @param prior.r Prior scale parameter (default: 0.353)
#'
#' @return Dataframe with columns: simid, index, raw (X values), raw2 (Y values), rw (cumulative r), bf, density.bf
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate correlation simulations
#' sims <- simcreate_cor(N = 100, rho = 0, n.sims = 100)
#' }
simcreate_cor <- function(N,
                         n.sims = 1000,
                         rho = 0,
                         method = c("pseudo", "quantis"),
                         parallel = TRUE,
                         nstart = 5,
                         exact = TRUE,
                         prior.r = 0.353) {

  method <- match.arg(method)

  if (method == "quantis" && parallel) {
    message("method='quantis' requires hardware device access and cannot be parallelized. Setting parallel=FALSE.")
    parallel <- FALSE
  }

  if (N < 3) stop("N must be at least 3")
  if (rho < -1 || rho > 1) stop("rho must be between -1 and 1")

  # Define single simulation function
  generate_simulation <- function(i) {
    # Generate bivariate data with correlation rho
    if (method == "pseudo") {
      # Standard method: bivariate normal
      x <- rnorm(N)
      y <- rho * x + sqrt(1 - rho^2) * rnorm(N)

    } else if (method == "quantis") {
      # Generate X and Y, then induce correlation
      x_uniform <- .generate_float_quantis(N, min = 0, max = 1)
      y_uniform <- .generate_float_quantis(N, min = 0, max = 1)

      x <- qnorm(x_uniform)
      y_indep <- qnorm(y_uniform)
      y <- rho * x + sqrt(1 - rho^2) * y_indep
    }

    # Create dataframe
    sim <- data.frame(
      simid = i,
      index = 1:N,
      raw = x,
      raw2 = y
    )

    # Calculate sequential correlations for "random walk"
    seq_cor <- sapply(nstart:N, function(n) {
      cor(x[1:n], y[1:n])
    })
    sim$rw <- c(rep(NA, nstart - 1), seq_cor)

    # Calculate BFs
    sim$bf <- .quiet(changeofevidence::bfcor(
      x, y,
      prior.r = prior.r,
      nstart = nstart,
      exact = exact
    )$BF)

    return(sim)
  }

  # Run simulations
  if (parallel) {
    cores <- parallel::detectCores() - 1
    message(sprintf("Generating %d correlation simulations (N=%d) using %d cores...", n.sims, N, cores))
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("generate_simulation", ".quiet", "N", "rho",
                                   "method", "nstart", "exact", "prior.r", ".generate_float_quantis"),
                           envir = environment())
    parallel::clusterEvalQ(cl, library(changeofevidence))
    sims_list <- pbapply::pblapply(1:n.sims, generate_simulation, cl = cl)
    parallel::stopCluster(cl)
  } else {
    message(sprintf("Generating %d correlation simulations (N=%d)...", n.sims, N))
    sims_list <- pbapply::pblapply(1:n.sims, generate_simulation)
  }

  # Combine and process
  sims_df <- do.call(rbind, sims_list)

  # FFT processing (suppress all messages and output)
  sims_df <- .quiet(simredo(sims_df, N, rw = TRUE))

  return(sims_df)
}


#' Download Pregenerated Quantum Random Data
#'
#' Downloads pregenerated quantum random bits from GitHub releases for use with method="files".
#' The file contains 1,000 simulations worth of quantum random bits (1 million bits each).
#' After download, the path is automatically saved and used by simulation functions.
#'
#' @param path Directory to save the file (default: package extdata directory)
#' @param force Re-download even if file exists (default: FALSE)
#'
#' @return Path to downloaded file (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' # Download once (file location is remembered automatically)
#' download_quantum_data()
#'
#' # Then use in simulations (finds file automatically)
#' result <- bfbinom(rbinom(100, 1, 0.6))
#' sims <- simcreate(result, n.sims = 1000, method = "files")
#' }
download_quantum_data <- function(path = NULL, force = FALSE) {
  # Use package extdata directory by default
  if (is.null(path)) {
    path <- system.file("extdata", package = "changeofevidence")
    if (path == "" || !dir.exists(path)) {
      # Package not installed or extdata doesn't exist, use user data directory
      if (.Platform$OS.type == "windows") {
        path <- file.path(Sys.getenv("APPDATA"), "changeofevidence")
      } else {
        path <- file.path(Sys.getenv("HOME"), ".changeofevidence")
      }
      message("Using user data directory: ", path)
    }
  }

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  dest_file <- file.path(path, "quantis_bits.rds")

  if (file.exists(dest_file) && !force) {
    message("Quantum data file already exists at: ", dest_file)
    message("Use force=TRUE to re-download.")
    # Save path for future use
    options(changeofevidence.quantum_data_path = dest_file)
    return(invisible(dest_file))
  }

  # URL to GitHub release
  # TODO: Update with actual release tag after creating release
  url <- "https://github.com/mrzdcmps/changeofevidence/releases/download/data-v1.0/quantum_bits_1000sims.rds"

  message("Downloading pregenerated quantum random data...")
  message("(1,000 simulations \u00d7 1,000,000 bits = ~125 MB)")
  message("This may take a few minutes...")

  tryCatch({
    download.file(url, dest_file, mode = "wb", quiet = FALSE)
    message("\nDownload complete!")
    message("File saved to: ", dest_file)

    # Save path in options for automatic detection
    options(changeofevidence.quantum_data_path = dest_file)
    message("Path saved for future use.")

  }, error = function(e) {
    stop(sprintf(
      "Download failed: %s\n\nPlease check:\n  1. Internet connection\n  2. GitHub release exists\n  3. Or generate manually with data-raw/generate_quantum_data.R",
      e$message
    ))
  })

  invisible(dest_file)
}
