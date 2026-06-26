# Change of Evidence Analysis for Random Walks
#
# This file provides coe_rw() and its helper functions (maxrw, energyrw),
# which mirror the coe() family but operate on random walk (rw) trajectories
# rather than Bayes Factor trajectories. FFT comparison uses the 'density.rw'
# column already stored in simulation dataframes from simcreate().


#' Maximum Random Walk Analysis
#'
#' Evaluates the likelihood of the highest absolute excursion reached by a random walk.
#'
#' The peak absolute value of the observed random walk is compared to those of all
#' simulations. An unusual peak can be assumed if fewer than 5% of simulations show
#' a larger absolute excursion at any point.
#'
#' @param rw A numeric vector containing the random walk trajectory.
#' @param sims.df A dataframe containing simulations with columns "simid" and "rw".
#' @return A list containing the maximum absolute RW value, its position, the signed
#'   value at that position, and the proportion of simulations with a higher peak.
#' @examples
#' \dontrun{
#' r.maxrw <- maxrw(cumsum(rbinom(100, 1, 0.5) * 2 - 1), sims)
#' }
#' @export
maxrw <- function(rw, sims.df = NULL) {
  if (is.null(sims.df)) stop("Please provide a simulation dataframe via sims.df.")
  if (!is.numeric(rw)) stop("rw must be a numeric vector!")
  if (!"rw" %in% names(sims.df)) stop("sims.df must contain a 'rw' column.")

  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)

  if (length(rw) != nrow(sims.df[sims.df$simid == simids[1], ]))
    stop("rw is not the same length as simulations!")

  sim.maxrw <- tapply(sims.df$rw, sims.df$simid, function(x) max(abs(x), na.rm = TRUE))
  max_abs_rw <- max(abs(rw), na.rm = TRUE)
  max_abs_rw_n <- which.max(abs(rw))
  sims_with_higher <- (sum(sim.maxrw >= max_abs_rw, na.rm = TRUE) / u.nsims) * 100

  out <- list(
    MaxRW = max_abs_rw,
    MaxRW_N = max_abs_rw_n,
    MaxRW_signed = rw[max_abs_rw_n],
    Sims_with_higher_RW = sims_with_higher,
    Num_Sims = u.nsims,
    data = rw
  )
  class(out) <- "coe_rw_result"
  return(out)
}


#' Random Walk Energy Analysis
#'
#' Calculates the energy (signed area under the random walk curve) and evaluates
#' its likelihood relative to simulations.
#'
#' Energy is calculated as the trapezoidal integral of the absolute random walk trajectory.
#' Using absolute values means both positive and negative excursions contribute equally,
#' and energy is always non-negative. The null energy for a random walk is 0.
#'
#' @param rw A numeric vector containing the random walk trajectory.
#' @param sims.df A dataframe containing simulations with columns "simid" and "rw".
#' @return A list containing the energy, simulation energy statistics, and the
#'   proportion of simulations with higher energy.
#' @examples
#' \dontrun{
#' r.energyrw <- energyrw(cumsum(rnorm(100)), sims)
#' }
#' @export
energyrw <- function(rw, sims.df = NULL) {
  if (is.null(sims.df)) stop("Please provide a simulation dataframe via sims.df.")
  if (!is.numeric(rw)) stop("rw must be a numeric vector!")
  if (!"rw" %in% names(sims.df)) stop("sims.df must contain a 'rw' column.")

  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)

  if (length(rw) != nrow(sims.df[sims.df$simid == simids[1], ]))
    stop("rw is not the same length as simulations!")

  .trapz_rw <- function(x) {
    pracma::trapz(as.numeric(seq_along(x)), abs(x))
  }

  sim.energy <- tapply(sims.df$rw, sims.df$simid, .trapz_rw)
  real_energy <- .trapz_rw(rw)
  sims_with_higher <- (sum(sim.energy >= real_energy, na.rm = TRUE) / u.nsims) * 100

  out <- list(
    Energy = real_energy,
    Simenergy_M = mean(sim.energy, na.rm = TRUE),
    Simenergy_SD = sd(sim.energy, na.rm = TRUE),
    Sims_with_higher_energy = sims_with_higher,
    Num_Sims = u.nsims,
    data = rw
  )
  class(out) <- "coe_rw_result"
  return(out)
}


#' Change of Evidence Analysis for Random Walks
#'
#' Performs a comprehensive analysis of a random walk trajectory by combining
#' Maximum RW, Energy RW, and FFT analyses, analogous to \code{coe()} for Bayes Factors.
#'
#' The function compares the observed random walk against the \code{rw} and
#' \code{density.rw} columns stored in simulation dataframes from \code{simcreate()}.
#' It automatically adjusts simulation length if simulations are longer than the data.
#'
#' @param rw A numeric vector containing the random walk trajectory.
#' @param sims.df A dataframe from \code{simcreate()} containing at minimum
#'   columns "simid", "rw", and "density.rw".
#' @return A list of class "coe_rw" containing results from maxrw, energyrw,
#'   and ffttest analyses, plus the harmonic mean p-value.
#' @examples
#' \dontrun{
#' result <- bfttest(x, mu = 0)
#' sims <- simcreate(result, n.sims = 1000)
#' rw <- cumsum(x - 0)
#' coe_rw_result <- coe_rw(rw, sims)
#' }
#' @export
coe_rw <- function(rw, sims.df) {
  if (!is.numeric(rw)) stop("rw must be a numeric vector.")
  if (!is.data.frame(sims.df)) stop("sims.df must be a dataframe.")
  if (!"simid" %in% names(sims.df)) stop("sims.df must contain a 'simid' column.")
  if (!"rw" %in% names(sims.df)) stop("sims.df must contain a 'rw' column.")
  if (!"density.rw" %in% names(sims.df)) stop("sims.df must contain a 'density.rw' column.")

  data_length <- length(rw)
  sim_ids <- unique(sims.df$simid)
  if (length(sim_ids) == 0) stop("No simulations found in sims.df.")

  sim_lengths <- sapply(sim_ids, function(id) sum(sims.df$simid == id))
  if (!all(sim_lengths == sim_lengths[1])) stop("Not all simulations have the same length.")
  sim_length <- sim_lengths[1]

  if (sim_length < data_length) {
    stop(paste0(
      "Simulation length (", sim_length, ") is shorter than rw length (",
      data_length, "). Cannot proceed."
    ))
  }

  if (sim_length > data_length) {
    message(paste0(
      "Simulation length (", sim_length, ") is longer than rw length (",
      data_length, "). Adjusting simulations using simredo."
    ))
    sims.df <- simredo(sims.df, data_length, rw = TRUE)
  }

  tryCatch({
    maxrw_result <- maxrw(rw, sims.df)
    energyrw_result <- energyrw(rw, sims.df)

    fftdata_rw <- fftcreate(rw)
    ffttest_result <- ffttest(fftdata_rw, sims.df, sims.df.col = "density.rw")

    p_maxrw  <- maxrw_result$Sims_with_higher_RW / 100
    p_energy <- energyrw_result$Sims_with_higher_energy / 100
    p_fft    <- ffttest_result$Sims_with_higher_Amplitude / 100

    p_values <- c(p_maxrw, p_energy, p_fft)
    if (any(is.na(p_values))) {
      warning("One or more p-values are NA. Check your simulation data.")
      harmonic_p <- NA
    } else if (any(p_values == 0)) {
      harmonic_p <- 0
    } else {
      harmonic_p <- 3 / sum(1 / p_values)
    }

    results <- list(
      Data_Length = data_length,
      Num_Sims = length(sim_ids),
      maxrw = maxrw_result,
      energyrw = energyrw_result,
      ffttest = ffttest_result,
      harmonic_p = harmonic_p
    )

    class(results) <- c("coe_rw", "list")
    return(results)

  }, error = function(e) {
    stop(paste("Error during analysis:", e$message))
  })
}


#' @export
#' @method print coe_rw_result
print.coe_rw_result <- function(x, ..., header = TRUE) {
  if (header) {
    cat("*** Change of Evidence (Random Walk) Results ***\n\n")
    if (!is.null(x$Num_Sims)) cat("Number of Simulations:", x$Num_Sims, "\n")
    cat("-----------------------------------\n")
  }

  if (!is.null(x$MaxRW)) {
    cat("Max |RW|:", round(x$MaxRW, 3), "(signed:", round(x$MaxRW_signed, 3), ") at N =", x$MaxRW_N, "\n")
    cat("Sims with \u2265 this |RW|:", x$Sims_with_higher_RW, "%\n")
  }

  if (!is.null(x$Energy)) {
    cat("Energy:", round(x$Energy, 3), "\n")
    cat("Simulated Energy: M =", round(x$Simenergy_M, 3), ", SD =", round(x$Simenergy_SD, 3), "\n")
    cat("Sims with \u2265 this Energy:", x$Sims_with_higher_energy, "%\n")
  }

  if (!is.null(x$Amplitude_sum)) {
    cat("Amplitude Sum:", round(x$Amplitude_sum, 3), "\n")
    cat("Simulated Amplitude Sum: M =", round(x$Sim_Ampsum_M, 3), ", SD =", round(x$Sim_Ampsum_SD, 3), "\n")
    cat("Sims with \u2265 this Amplitude:", x$Sims_with_higher_Amplitude, "%\n")
  }

  invisible(x)
}


#' @export
#' @method print coe_rw
print.coe_rw <- function(x, ...) {
  cat("*** Change of Evidence (Random Walk) Results ***\n\n")
  if (!is.null(x$Num_Sims)) cat("Number of Simulations:", x$Num_Sims, "\n")
  if (!is.null(x$harmonic_p)) cat("CoE RW p-value:", round(x$harmonic_p, 6), "\n")
  cat("-----------------------------------\n")

  if (!is.null(x$maxrw)) {
    cat("Max |RW| Test\n\n")
    print(x$maxrw, header = FALSE)
    cat("-----------------------------------\n")
  }

  if (!is.null(x$energyrw)) {
    cat("RW Energy Test\n\n")
    print(x$energyrw, header = FALSE)
    cat("-----------------------------------\n")
  }

  if (!is.null(x$ffttest)) {
    cat("FFT Test (density.rw)\n\n")
    print.coe(x$ffttest, header = FALSE)
    cat("-----------------------------------\n")
  }

  invisible(x)
}
