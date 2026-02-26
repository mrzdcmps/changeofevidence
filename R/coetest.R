#' Maximum BF Analysis
#'
#' This function evaluates the likelihood of the highest reached Bayes Factor.
#'
#' The highest reached (Maximum) BF of the experimental data is compared to those of all simulations.
#' An unusual peak can be assumed if less than 5\% of simulations show a higher BF at any time.
#'
#' @param data A a seqbf object or vector containing Bayes Factors.
#' @param sims.df A dataframe containing simulations, including columns "simid" and "bf".
#' @return A list containing the Maximum BF of the experimental data, its position in the temporal order of data points, and the proportion of simulations with higher BFs at any time.
#' @examples
#' \dontrun{
#' r.maxbf <- maxbf(seqbf, sims)
#' r.maxbf <- maxbf(tbl$bf, sims.df = newsims)
#' }
#' @export

# Maximum BF analysis
maxbf <- function(data, sims.df = NULL) {
  if (inherits(data, "seqbf")) data <- data$BF
  if (is.null(sims.df)) stop("Please provide a simulation dataframe via sims.df.")

  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)

  if (!is.numeric(data)) stop("Data must be a numeric vector!")
  if (length(data) != nrow(sims.df[sims.df$simid == simids[1], ])) stop("Data is not the same length as simulations!")

  # Replace early NAs with BF=1
  data <- .replace_early_na(data)
  if ("bf" %in% names(sims.df)) {
    sims.df$bf <- unlist(tapply(sims.df$bf, sims.df$simid, .replace_early_na))
  }
  
  sim.maxbf <- tapply(sims.df$bf, sims.df$simid, max, na.rm = TRUE)
  max_data_bf <- max(data)
  max_data_bf_n <- which.max(data)
  sims_with_higher_bf <- (sum(sim.maxbf >= max_data_bf, na.rm = TRUE) / u.nsims) * 100  # Percentage
  
  maxbf.out <- list(
    MaxBF = max_data_bf,
    MaxBF_N = max_data_bf_n,
    Sims_with_higher_BFs = sims_with_higher_bf,
    Num_Sims = u.nsims,
    data = data
  )
  
  class(maxbf.out) <- "coe"
  return(maxbf.out)
}


### Helper function: Replace early NAs with BF=1
.replace_early_na <- function(bf_vector) {
  if (!any(is.na(bf_vector))) {
    return(bf_vector)
  }
  # Find the last position with NA
  last_na_pos <- max(which(is.na(bf_vector)))
  # Replace all values up to and including last NA with 1
  bf_vector[1:last_na_pos] <- 1
  return(bf_vector)
}

### Helper function: Calculate Energy
.energycount <- function(data, nullenergy = nullenergy){
  pracma::trapz(as.numeric(1:length(data)), data)-nullenergy
}

#' BF Energy Analysis
#'
#' This function calculates the Energy of the evidence over time and evaluates its likelihood.
#'
#' The Energy contains all areas above BF=1 minus all areas below BF=1. A positive energy therefore indicates an overall orientation of the evidence towards H1 over data collection.
#' Energy is calculated using the trapezoidal integration "trapz" of the pracma-library.
#' Mean and SD scores of the simulations' energies are provided.
#'
#' @param data A seqbf object or vector containing Bayes Factors.
#' @param sims.df A dataframe containing simulations, including columns "simid" and "bf".
#' @return A list containing the Energy of the experimental data, the mean and SD of the energy values of all simulations, and the proportion of simulations with a higher energy than the experimental data.
#' @examples
#' \dontrun{
#' r.energybf <- energybf(seqbf, sims)
#' r.energybf <- energybf(tbl$bf, sims.df = newsims)
#' }
#' @export

# Energy of BF
energybf <- function(data, sims.df = NULL) {
  if (inherits(data, "seqbf")) data <- data$BF
  if (is.null(sims.df)) stop("Please provide a simulation dataframe via sims.df.")

  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)
  nullenergy <- length(data) - 1

  if (!is.numeric(data)) stop("Data must be a numeric vector!")
  if (length(data) != nrow(sims.df[sims.df$simid == simids[1], ])) stop("Data is not the same length as simulations!")

  # Replace early NAs with BF=1
  data <- .replace_early_na(data)
  if ("bf" %in% names(sims.df)) {
    sims.df$bf <- unlist(tapply(sims.df$bf, sims.df$simid, .replace_early_na))
  }
  
  #sim.energy <- pbapply::pbtapply(sims.df$bf, sims.df$simid, .energycount, nullenergy)
  sim.energy <- tapply(sims.df$bf, sims.df$simid, .energycount, nullenergy)
  real_energy <- pracma::trapz(as.numeric(1:length(data)), data) - nullenergy

  sims_with_higher_energy <- (sum(sim.energy >= real_energy, na.rm = TRUE) / u.nsims) * 100  # Percentage
  
  energybf.out <- list(
    Energy = real_energy,
    Simenergy_M = mean(sim.energy),
    Simenergy_SD = sd(sim.energy),
    Sims_with_higher_energy = sims_with_higher_energy,
    Num_Sims = u.nsims,
    data = data
  )
  
  class(energybf.out) <- "coe"
  return(energybf.out)
}


#' Create a Fast Fourier transformed vector
#'
#' This function converts input signals into its Fast Fourier transform (FFT).
#'
#' The FFT shows the spectral density of the input. It can be understood as a harmonic analysis and indicates the strength of osciallative components comprising the signal.
#' It can therefore be used to assess the meaning of oscillation as a characteristic. N/2 Frequencies are used as sampling rate.
#'
#' @param data A seqbf object or vector to be transformed.
#' @examples
#' \dontrun{
#' density.bf <- fftcreate(seqbf)
#' tblFFT <- data.frame(density.bf = fftcreate(tbl$bf), density.rw = fftcreate(tbl$rw))
#' }
#' @export

# Create FFT
fftcreate <- function(data){
  
  if(inherits(data,"seqbf") == TRUE) data <- data$BF
  
  L <- length(data) # Length
  L2 <- ceiling(L/2) # Consider only first half, round up
  T <- 1/L # Rate
  Y <- fft(data) # Fast Fourier Transformation
  abs(Y/L)[1:L2]
}

### Helper-Function: Create FFT without cutting frequencies
.fftcreatefull <- function(data){
  L <- length(data) # Length
  T <- 1/L # Rate
  Y <- fft(data) # Fast Fourier Transformation
  abs(Y/L)
}

### Helper-Function: Count Frequencies
.fftcount <- function(data, sims.df = NULL, sims.df.col = "density.bf"){
  
  CHz <- sims.df[sims.df$index==as.numeric(data[2]),]
  output <- data.frame(H = as.numeric(data[2]), LowerSims = (sum(CHz[[sims.df.col]] < as.numeric(data[1])))/length(unique(sims.df$simid)))
  output
  
}

#' Frequency Analysis Test
#'
#' This function analyses FFTs and compares the amplitudes of all frequencies to those of simulations.
#'
#' If you want to use the old "Top5-Frequency" method instead of amplitude sums, indicate by setting 'top5 = TRUE'.
#'
#' @param data A seqbf-object or a vector containing Fourier transformed (spectral density) data (use `fftcreate` function).
#' @param sims.df A dataframe containing simulations, including columns "index" and "simid".
#' @param sims.df.col The column of the simulation dataframe that contains the comparison data.
#' @param top5 Logical. If set to TRUE, function will additionally return the Top5-Frequency method. For each frequency the function counts how many simulations show a higher amplitude. If no more than 5 percent of simulations are above the experimental value, it is considered a "Top5-Frequency". The proportion of Top5-Frequencies indicates the pronouncedness of oscillatory elements in the data.
#' @return A list containing a dataframe of all frequencies and the proportion of simulations with a lower amplitude, and information on the sum of amplitudes.
#' @examples
#' \dontrun{
#' ffttest(bf, sims)
#' r.fft <- ffttest(fftcreate(seqbf))
#' r.fftrw <- ffttest(tblFFT$density.rw, sims.df = newsims, sims.df.col = "density.rw")
#' }
#' @export

ffttest <- function(data, sims.df = NULL, sims.df.col = "density.bf", top5 = FALSE) {
  if (inherits(data, "seqbf")) data <- changeofevidence::fftcreate(data$BF)
  if (is.null(sims.df)) stop("Please provide a simulation dataframe via sims.df.")

  if (var(data[1:3]) == 0) stop("It seems like you specified a vector containing BFs. Please use fftcreate(bf) to Fourier transform first.")
  if (!is.numeric(sims.df[[sims.df.col]])) stop("Wrong sims data. Does sims.df.col exist?")
  if (ceiling(max(sims.df$index / 2)) != length(data)) stop("Lengths of FFT data and sims do not match! Consider using simredo()")
  
  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)
  
  data.df <- data.frame(data = data, H = seq_along(data))
  ampsum <- sum(data.df$data)
  
  sims.df <- sims.df[sims.df$index <= length(data), ]
  simampsum <- tapply(sims.df[[sims.df.col]], sims.df$simid, sum)

  sims_with_higher_amp <- (sum(ampsum <= simampsum, na.rm = TRUE) / u.nsims) * 100  # Percentage
  
  fftbf.out <- list(
    Amplitude_sum = ampsum,
    Sims_with_higher_Amplitude = sims_with_higher_amp,
    Sim_Ampsum_M = mean(simampsum),
    Sim_Ampsum_SD = sd(simampsum),
    Num_Sims = u.nsims,
    data = data
  )
  
  if (top5) {
    list.HB <- pbapply::pbapply(data.df, 1, .fftcount, sims.df = sims.df, sims.df.col = sims.df.col)
    list.HB <- do.call(rbind, list.HB)
    fftbf.out$FFTComparison <- list.HB
    fftbf.out$Top5_Frequencies_above_95pct <- sum(list.HB$LowerSims > 0.95) / length(data) * 100
  }
  
  class(fftbf.out) <- "coe"
  return(fftbf.out)
}



# Create a copy of simulations and only redo FFT
#' Create a copy of simulation-df with different amount of trials.
#'
#' This function copies the BFs of an already computated simulation dataframe and re-calculates only the density via a FFT.
#'
#' Computating BFs for monte-carlo-simulations can take a very long time. You can reuse an already calculated simulation-df and cut off the amount of trials per simulation. This function re-calculates only the densities of BF and RW via the Fast Fourier transformation.
#'
#' @param df A dataframe containing monte-carlo simulations. Must include the columns 'bf', 'density.bf', 'index', and 'simid'.
#' @param n The amount of trials per simulation. Must be smaller than in the original simulation df.
#' @param rw boolean. Set to TRUE if you also want to recalculate the density of the Random Walk. Needs the colums 'rw' and 'density.rw'.
#' @examples
#' \dontrun{
#' sim.new <- simredo(sims, length(bf.new), rw=F)
#' }
#' @export
simredo <- function(df, n, rw = TRUE){
  if(n > max(df$index)) stop("Amount of trials (n) is larger than in the provided simulation dataframe!")

  # Preserve attributes before subsetting
  attrs_to_preserve <- attributes(df)
  attrs_to_preserve$names <- NULL  # Will be set by subset
  attrs_to_preserve$row.names <- NULL  # Will be set by subset

  df.new <- subset(df, index <= n)
  cat("Recalculating BF FFT:\n")
  df.new$density.bf <- unlist(pbapply::pbtapply(df.new$bf, df.new$simid, .fftcreatefull, simplify = T))
  if(rw == T) {
    cat("Recalculating RW FFT:\n")
    df.new$density.rw <- unlist(pbapply::pbtapply(df.new$rw, df.new$simid, .fftcreatefull, simplify = T))
  }

  # Restore attributes (updating N to reflect new length)
  if (!is.null(attrs_to_preserve$N)) {
    attrs_to_preserve$N <- n
  }
  for (attr_name in names(attrs_to_preserve)) {
    if (!attr_name %in% c("names", "row.names", "class")) {
      attr(df.new, attr_name) <- attrs_to_preserve[[attr_name]]
    }
  }

  # Restore coe_sim class if it was present
  if ("coe_sim" %in% attrs_to_preserve$class) {
    class(df.new) <- c("coe_sim", "data.frame")
  }

  df.new
}


#' Change of Evidence Analysis
#'
#' This function performs a comprehensive analysis of Bayesian evidence over time by 
#' combining Maximum BF, Energy BF, and FFT analyses.
#'
#' The function validates that simulation data matches the experimental data length,
#' automatically adjusting simulations if they are longer than needed, and performs
#' all three core analyses: maximum BF evaluation, energy calculation, and frequency
#' analysis via FFT.
#'
#' @param data A seqbf object or vector containing Bayes Factors.
#' @param sims.df A dataframe containing simulations, must include column "simid".
#' @return A list containing results from maxbf, energybf, and ffttest analyses.
#' @examples
#' \dontrun{
#' result <- coe(seqbf, sims)
#' result <- coe(bf_vector, sims.df = newsims)
#' }
#' @export

coe <- function(data, sims.df) {
  # Step 1: Handle seqbf objects
  if (inherits(data, "seqbf")) {
    data <- data$BF
  }

  # Validate inputs
  if (!is.numeric(data)) {
    stop("data must be a numeric vector or seqbf object")
  }

  if (!is.data.frame(sims.df)) {
    stop("sims.df must be a dataframe")
  }

  if (!"simid" %in% names(sims.df)) {
    stop("sims.df must contain a 'simid' column")
  }

  # Step 1a: Replace early NAs with BF=1 in empirical data
  data <- .replace_early_na(data)

  # Step 1b: Replace early NAs with BF=1 in simulation data
  if ("bf" %in% names(sims.df)) {
    sims.df$bf <- unlist(tapply(sims.df$bf, sims.df$simid, .replace_early_na))
  }
  
  # Step 2: Check if sims.df matches the data length
  data_length <- length(data)
  
  # Get the number of rows for a single simulation
  sim_ids <- unique(sims.df$simid)
  if (length(sim_ids) == 0) {
    stop("No simulations found in sims.df")
  }
  
  # Check the length of the first simulation
  first_sim_length <- sum(sims.df$simid == sim_ids[1])
  
  # Validate all simulations have the same length
  sim_lengths <- sapply(sim_ids, function(id) sum(sims.df$simid == id))
  if (!all(sim_lengths == first_sim_length)) {
    stop("Not all simulations have the same length")
  }
  
  sim_length <- first_sim_length
  
  # Step 2a: If simulation is shorter than data, stop
  if (sim_length < data_length) {
    stop(paste("Simulation length (", sim_length, ") is shorter than data length (", 
               data_length, "). Cannot proceed.", sep = ""))
  }
  
  # Step 2b: If simulation is longer than data, use simredo
  if (sim_length > data_length) {
    message(paste("Simulation length (", sim_length, ") is longer than data length (", 
                  data_length, "). Adjusting simulations using simredo.", sep = ""))
    sims.df <- simredo(sims.df, data_length, rw = FALSE)
  }
  
  # Step 3: Perform all three analyses
  tryCatch({
    # Maximum BF analysis
    maxbf_result <- maxbf(data, sims.df)
    
    # Energy BF analysis  
    energybf_result <- energybf(data, sims.df)
    
    # FFT analysis
    fftdata <- fftcreate(data)
    ffttest_result <- ffttest(fftdata, sims.df)

    # Step 4: Calculate harmonic mean of empirical p-values
    # Convert percentages to proportions (p-values)
    p_maxbf <- maxbf_result$Sims_with_higher_BFs / 100
    p_energy <- energybf_result$Sims_with_higher_energy / 100
    p_fft <- ffttest_result$Sims_with_higher_Amplitude / 100

    # Calculate harmonic mean: n / sum(1/p_i)
    # Handle edge cases where any p-value is 0 or NA
    p_values <- c(p_maxbf, p_energy, p_fft)
    if (any(is.na(p_values))) {
      warning("One or more p-values are NA. Check your simulation data.")
      harmonic_p <- NA
    } else if (any(p_values == 0)) {
      harmonic_p <- 0
    } else {
      harmonic_p <- 3 / sum(1 / p_values)
    }

    # Step 5: Store all results in a list
    results <- list(
      Data_Length = data_length,
      Num_Sims = length(sim_ids),
      maxbf = maxbf_result,
      energybf = energybf_result,
      ffttest = ffttest_result,
      harmonic_p = harmonic_p
    )
    
    # Add class for potential future methods
    class(results) <- c("coe", "list")
    
    return(results)
    
  }, error = function(e) {
    stop(paste("Error during analysis:", e$message))
  })
}


#' @export
#' @method print coe
print.coe <- function(x, ..., header = TRUE) {
  # Only top-level prints the header
  if (header) {
    cat("*** Change of Evidence Results ***\n\n")
    if (!is.null(x$Num_Sims)) {
      cat("Number of Simulations:", x$Num_Sims, "\n")
    }
    if (!is.null(x$harmonic_p)) {
      cat("CoE p-value:", round(x$harmonic_p, 6), "\n")
    }
    cat("-----------------------------------\n")
  }
  
  # Print MaxBF info
  if (!is.null(x$MaxBF)) {
    cat("Max BF:", round(x$MaxBF, 3), "at N =", x$MaxBF_N, "\n")
    cat("Sims with \u2265 this BF:", x$Sims_with_higher_BFs, "%\n")
  }
  
  # Print Energy info
  if (!is.null(x$Energy)) {
    cat("Energy:", round(x$Energy, 3), "\n")
    cat("Simulated Energy: M =", round(x$Simenergy_M, 3), ", SD =", round(x$Simenergy_SD, 3), "\n")
    cat("Sims with \u2265 this Energy:", x$Sims_with_higher_energy, "%\n")
  }
  
  # Print FFT info
  if (!is.null(x$Amplitude_sum)) {
    cat("Amplitude Sum:", round(x$Amplitude_sum, 3), "\n")
    cat("Simulated Amplitude Sum: M =", round(x$Sim_Ampsum_M, 3), ", SD =", round(x$Sim_Ampsum_SD, 3), "\n")
    cat("Sims with \u2265 this Amplitude:", x$Sims_with_higher_Amplitude, "%\n")
    if (!is.null(x$Top5_Frequencies_above_95pct)) {
      cat("Top 5 Frequencies > 95% of Simulations:", x$Top5_Frequencies_above_95pct, "%\n")
    }
  }
  
  # Composite object \u2014 call sub-prints without headers
  if (!is.null(x$maxbf)) {
    cat("MaxBF Test\n\n")
    print(x$maxbf, header = FALSE)
    cat("-----------------------------------\n")
  }
  
  if (!is.null(x$energybf)) {
    cat("BF Energy Test\n\n")
    print(x$energybf, header = FALSE)
    cat("-----------------------------------\n")
  }
  
  if (!is.null(x$ffttest)) {
    cat("FFT Test\n\n")
    print(x$ffttest, header = FALSE)
    cat("-----------------------------------\n")
  }
  
  invisible(x)
}



#' @export
#' @method plot coe
plot.coe <- function(x, sims.df = NULL, ...) {
  if (!inherits(x, "coe")) {
    stop("Object must be of class 'coe' to plot.")
  }

  # Helper function to convert percentage to p-value
  percentile_to_pvalue <- function(percentile) {
    percentile / 100
  }

  # Aggregate coe object — build plots
  if (!is.null(x$maxbf) || !is.null(x$energybf) || !is.null(x$ffttest)) {

    common_theme <- ggplot2::theme(
      text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 10, face = "bold")
    )

    plots <- list()

    # Panel 1: BF trajectory with MaxBF marked
    if (!is.null(x$maxbf)) {
      p_maxbf <- percentile_to_pvalue(x$maxbf$Sims_with_higher_BFs)
      # Only add simulations if sims argument is provided
      # Filter out non-positive BF values to avoid log scale issues
      if (!is.null(sims.df)) {
        sims_subset <- subset(sims.df, simid <= 100 & !is.na(bf) & bf > 0)
      } else {
        sims_subset <- NULL
      }
      p1 <- plotbf(x$maxbf$data, sims.df = sims_subset, show_annotations = FALSE) +
        ggplot2::geom_point(aes(x = x$maxbf$MaxBF_N, y = x$maxbf$MaxBF),
                           color = "red", size = 3) +
        ggplot2::ggtitle(sprintf("Sequential BF (Max = %.2f at n = %d)",
                                x$maxbf$MaxBF, x$maxbf$MaxBF_N)) +
        ggplot2::theme(legend.position = "none") +
        common_theme
      plots <- c(plots, list(p1))
    }


    # Panel 4: Ridgeline plot of all three distributions
    if (!is.null(x$harmonic_p) && !is.null(sims.df)) {
      # Calculate simulation max BFs
      sim_maxbfs <- tapply(sims.df$bf, sims.df$simid, function(x) {
        valid_bfs <- x[!is.na(x) & is.finite(x) & x > 0]
        if (length(valid_bfs) > 0) max(valid_bfs) else NA
      })
      sim_maxbfs <- sim_maxbfs[!is.na(sim_maxbfs) & is.finite(sim_maxbfs) & sim_maxbfs > 0]

      # Calculate simulation energies
      nullenergy <- length(x$energybf$data) - 1
      sims_clean <- sims.df
      sims_clean$bf <- unlist(tapply(sims.df$bf, sims.df$simid, .replace_early_na))
      sim_energies <- tapply(sims_clean$bf, sims_clean$simid, function(bf) {
        pracma::trapz(as.numeric(1:length(bf)), bf) - nullenergy
      })
      sim_energies <- sim_energies[is.finite(sim_energies)]

      # Calculate simulation FFT amplitudes
      sim_amplitudes <- tapply(sims.df$density.bf, sims.df$simid, sum, na.rm = TRUE)
      sim_amplitudes <- sim_amplitudes[is.finite(sim_amplitudes)]

      # Convert to percentile ranks (0-100 scale)
      maxbf_percentiles <- ecdf(sim_maxbfs)(sim_maxbfs) * 100
      energy_percentiles <- ecdf(sim_energies)(sim_energies) * 100
      fft_percentiles <- ecdf(sim_amplitudes)(sim_amplitudes) * 100

      # Calculate observed percentiles
      obs_maxbf_pct <- ecdf(sim_maxbfs)(x$maxbf$MaxBF) * 100
      obs_energy_pct <- ecdf(sim_energies)(x$energybf$Energy) * 100
      obs_fft_pct <- ecdf(sim_amplitudes)(x$ffttest$Amplitude_sum) * 100

      # Create data frame for ridgeline plot
      ridge_data <- data.frame(
        percentile = c(maxbf_percentiles, energy_percentiles, fft_percentiles),
        test = factor(
          rep(c("MaxBF", "Energy", "FFT"),
              c(length(maxbf_percentiles), length(energy_percentiles), length(fft_percentiles))),
          levels = c("FFT", "Energy", "MaxBF")  # Reverse order for bottom-to-top display
        )
      )

      # Observed values for annotations
      obs_data <- data.frame(
        test = factor(c("MaxBF", "Energy", "FFT"), levels = c("FFT", "Energy", "MaxBF")),
        percentile = c(obs_maxbf_pct, obs_energy_pct, obs_fft_pct),
        p_value = c(
          percentile_to_pvalue(x$maxbf$Sims_with_higher_BFs),
          percentile_to_pvalue(x$energybf$Sims_with_higher_energy),
          percentile_to_pvalue(x$ffttest$Sims_with_higher_Amplitude)
        )
      )

      # Adjust text position based on percentile to avoid clipping
      obs_data$hjust_val <- ifelse(obs_data$percentile > 80, 1.15, -0.15)

      p4 <- ggplot2::ggplot(ridge_data, aes(x = percentile, y = test)) +
        ggridges::geom_density_ridges(
          aes(fill = test),
          alpha = 0.7,
          scale = 0.9,
          rel_min_height = 0.01
        ) +
        ggplot2::geom_segment(
          data = obs_data,
          aes(x = percentile, xend = percentile, y = as.numeric(test), yend = as.numeric(test) + 0.85),
          color = "red",
          linewidth = 1.2,
          linetype = "solid"
        ) +
        ggplot2::geom_text(
          data = obs_data,
          aes(x = percentile, y = as.numeric(test) + 0.5,
              label = sprintf("p=%.3f", p_value),
              hjust = hjust_val),
          size = 3,
          color = "red",
          fontface = "bold"
        ) +
        ggplot2::scale_fill_manual(
          values = c("MaxBF" = "#619CFF", "Energy" = "#00BA38", "FFT" = "#F8766D")
        ) +
        ggplot2::labs(
          x = "Percentile Rank",
          y = "",
          title = sprintf("CoE Distributions (Harmonic p = %.3f)", x$harmonic_p)
        ) +
        ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        common_theme

      plots <- c(plots, list(p4))
    } else if (!is.null(x$harmonic_p)) {
      # Fallback: bar chart if no sims available
      p_data <- data.frame(
        Test = c("MaxBF", "Energy", "FFT", "Harmonic\nMean"),
        p_value = c(
          percentile_to_pvalue(x$maxbf$Sims_with_higher_BFs),
          percentile_to_pvalue(x$energybf$Sims_with_higher_energy),
          percentile_to_pvalue(x$ffttest$Sims_with_higher_Amplitude),
          x$harmonic_p
        ),
        Type = c("Individual", "Individual", "Individual", "Combined")
      )
      p_data$Test <- factor(p_data$Test, levels = p_data$Test)

      p4 <- ggplot2::ggplot(p_data, aes(x = Test, y = p_value, fill = Type)) +
        ggplot2::geom_col(color = "black", width = 0.7) +
        ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed",
                           color = "red", linewidth = 0.8) +
        ggplot2::geom_text(aes(label = sprintf("%.3f", p_value)),
                          vjust = -0.5, size = 3) +
        ggplot2::scale_fill_manual(values = c("Individual" = "steelblue",
                                               "Combined" = "darkorange")) +
        ggplot2::labs(x = "", y = "p-value",
                     title = "CoE p-values") +
        ggplot2::ylim(0, max(0.15, max(p_data$p_value) * 1.2)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        common_theme
      plots <- c(plots, list(p4))
    }

    # Arrange plots (typically 2: BF trajectory + ridgeline)
    if (length(plots) == 2) {
      combined <- plots[[1]] / plots[[2]]
    } else if (length(plots) == 1) {
      combined <- plots[[1]]
    } else {
      stop("No plots available to display.")
    }

    print(combined)
    return(invisible())
  }
}



