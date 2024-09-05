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
#' r.maxbf <- maxbf(seqbf, sims)
#' r.maxbf <- maxbf(tbl$bf, sims.df = newsims)
#' @export

# Maximum BF analysis
maxbf <- function(data, sims.df = sims) {
  
  # Handle seqbf input
  if (inherits(data, "seqbf")) data <- data$BF
  
  # Precompute unique simulation IDs and their count
  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)
  
  # Input validation
  if (!is.numeric(data)) stop("Data must be a numeric vector!")
  if (length(data) != nrow(sims.df[sims.df$simid == simids[1],])) stop("Data is not the same length as simulations!")
  
  # Compute maximum BF for both the real data and simulations
  sim.maxbf <- tapply(sims.df$bf, sims.df$simid, max, na.rm = TRUE)
  max_bf <- max(data)
  max_bf_N <- which.max(data)
  sims_with_higher_bf <- (sum(sim.maxbf >= max_bf) / u.nsims) * 100
  
  # Return result as a list and assign the custom class 'coe'
  maxbf.out <- list(
    MaxBF = max_bf,
    MaxBF_N = max_bf_N,
    Sims_with_higher_BF = sims_with_higher_bf,
    Total_Sims = u.nsims
  )
  
  class(maxbf.out) <- "coe"
  return(maxbf.out)
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
#' r.energybf <- energybf(seqbf, sims)
#' r.energybf <- energybf(tbl$bf, sims.df = newsims)
#' @export

# Energy of BF
energybf <- function(data, sims.df = sims) {
  
  # Handle seqbf input
  if (inherits(data, "seqbf")) data <- data$BF
  
  # Precompute unique simulation IDs and their count
  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)
  
  # Input validation
  if (!is.numeric(data)) stop("Data must be a numeric vector!")
  if (length(data) != nrow(sims.df[sims.df$simid == simids[1],])) stop("Data is not the same length as simulations!")
  
  # Energy calculation for both real data and simulations
  nullenergy <- length(data) - 1
  sim.energy <- pbapply::pbtapply(sims.df$bf, sims.df$simid, .energycount, nullenergy)
  real.energy <- pracma::trapz(as.numeric(1:length(data)), data) - nullenergy
  sims_with_higher_energy <- (sum(sim.energy >= real.energy) / u.nsims) * 100
  
  # Return result as a list and assign the custom class 'coe'
  energybf.out <- list(
    Energy = real.energy,
    Simenergy_M = mean(sim.energy),
    Simenergy_SD = sd(sim.energy),
    Sims_with_higher_Energy = sims_with_higher_energy,
    Total_Sims = u.nsims
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
#' density.bf <- fftcreate(seqbf)
#' tblFFT <- data.frame(density.bf = fftcreate(tbl$bf), density.rw = fftcreate(tbl$rw))
#' @export

# Create FFT
fftcreate <- function(data){
  
  if(inherits(data,"seqbf") == TRUE) data <- data$BF
  
  L <- length(data) # Länge
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
.fftcount <- function(data, sims.df = sims, sims.df.col = "density.bf"){
  
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
#' ffttest(bf, sims)
#' r.fft <- ffttest(fftcreate(seqbf))
#' r.fftrw <- ffttest(tblFFT$density.rw, sims.df = newsims, sims.df.col = "density.rw")
#' @export

ffttest <- function(data, sims.df = sims, sims.df.col = "density.bf", top5 = FALSE) {
  
  # Handle seqbf input
  if (inherits(data, "seqbf")) data <- changeofevidence::fftcreate(data$BF)
  
  # Validate inputs
  if (var(data[1:3]) == 0) stop("It seems like you specified a vector containing BFs. Please use fftcreate(bf) to Fourier transform first.")
  if (!is.numeric(sims.df[[sims.df.col]])) stop("Wrong sims data. Does sims.df.col exist?")
  if (ceiling(max(sims.df$index / 2)) != length(data)) stop("Lengths of FFT data and sims do not match! Consider using simredo()")
  
  # Precompute unique simulation IDs and their count
  simids <- unique(sims.df$simid)
  u.nsims <- length(simids)
  
  # Subset sims.df to match data length
  sims.df <- sims.df[sims.df$index <= length(data), ]
  
  # Compute the amplitude sum of real data
  ampsum <- sum(data)
  
  # Compute the amplitude sum of simulations
  simampsum <- tapply(sims.df[[sims.df.col]], sims.df$simid, sum)
  
  # Optionally, analyze top 5 frequencies
  if (top5) {
    data.df <- data.frame(data = data, H = seq_along(data))
    list.HB <- pbapply::pbapply(data.df, 1, .fftcount, sims.df = sims.df, sims.df.col = sims.df.col)
    list.HB <- dplyr::bind_rows(list.HB)
    top5_above_95 <- sum(list.HB$LowerSims > 0.95) / length(data)
  }
  
  # Prepare output
  fftbf.out <- list(
    Amplitude_sum = ampsum,
    Sims_with_higher_Amplitude = sum(ampsum <= simampsum) / u.nsims * 100,
    Sim_Ampsum_M = mean(simampsum),
    Sim_Ampsum_SD = sd(simampsum),
    Total_Sims = u.nsims
  )
  
  if (top5) {
    fftbf.out$Frequencies_above_95 <- top5_above_95
  }
  
  class(fftbf.out) <- "coe"
  return(fftbf.out)
}



#' Frequency Analysis Likelihood Test
#'
#' This function estimates the likelihood of a specific proportion of Top5-Frequencies.
#' To serve as reference, the Top5-Frequencies of 10% of simulations in comparison to all simulations are calculated.
#'
#' @param df A dataframe containing the results of a ffttest.
#' @param proportion The porportion of simulations in sims.df that should be tested against all simulations.
#' @param sims.df A dataframe containing simulations, including columns "index" and "simid".
#' @param sims.df.col The column of the simulation dataframe that contains the comparison data.
#' @return A vector containing the number of Top5-Frequencies of simulations.
#' @examples
#' fftlikelihood(r.fft)


# Count likelihood and distribution of Top5 occurrences
fftlikelihood <- function(df, proportion = 100, sims.df = sims, sims.df.col = "density.bf"){
  warning("This function is deprecated and will be removed in a future version. Please use `ffttest()` instead.")
  require(foreach)
  require(doParallel)
  
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  cat(">> FREQUENCY ANALYSIS LIKELIHOOD << \n")
  cat("This test runs in parallel. See log.txt for status updates!")
  
  likelihoodlist <- foreach(i=1:(length(unique(sims.df$simid))*(proportion/100)), .combine=c) %dopar% {
    sink("log.txt", append=TRUE)  
    cat(paste(Sys.time(),"Starting iteration",i,"of",length(unique(sims.df$simid))*(proportion/100),"\n"))
    sink()
    tmpdat.r <- subset(sims.df, simid == i)
    tmptest <- changeofevidence::ffttest(tmpdat.r[[sims.dfl.col]], sims.df)
    likeresult <- sum(tmptest$LowerSims > 0.95)
    likeresult
  }
  stopCluster(cl)
  
  cat("\nLikelihood for",sum(df$LowerSims > 0.95),"(=",(sum(df$LowerSims > 0.95)/nrow(df))*100,"%) or more Top5-Frequencies is estimated",(sum(likelihoodlist >= sum(df$LowerSims > 0.95))/length(likelihoodlist))*100,"% \n")
  
  likelihood.out <- list("Likelihood" = sum(likelihoodlist >= sum(df$LowerSims > 0.95))/length(likelihoodlist), "Top5 of Sims" = likelihoodlist)
  return(likelihood.out)
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
#' sim.new <- simredo(sims, length(bf.new), rw=F)
#' @export
simredo <- function(df, n, rw = TRUE){
  if(n > max(df$index)) stop("Amount of trials (n) is larger than in the provided simulation dataframe!")
  df.new <- subset(df, index <= n)
  cat("Recalculating BF FFT:\n")
  df.new$density.bf <- unlist(pbapply::pbtapply(df.new$bf, df.new$simid, .fftcreatefull, simplify = T))
  if(rw == T) {
    cat("Recalculating RW FFT:\n")
    df.new$density.rw <- unlist(pbapply::pbtapply(df.new$rw, df.new$simid, .fftcreatefull, simplify = T))
  }
  df.new
}


# # Print Method for COE tests
# print.coe <- function(x, ...) {
#   if (!inherits(x, "coe")) stop("Object is not of class 'coe'")
#   
#   # Print Maximum BF Output
#   if (!is.null(x$MaxBF)) {
#     cat("COE: Maximum BF test\n")
#     cat("Highest BF:", x$MaxBF, "(at N =", x$MaxBF_N, ")\n")
#     cat("Sims with this or higher BFs:", x$Sims_with_higher_BFs, "%\n")
#     cat("Number of simulations used:", x$n_sims, "\n")
#   }
#   
#   # Print Energy BF Output
#   if (!is.null(x$Energy)) {
#     cat("COE: BF Energy test\n")
#     cat("Energy of BF of data:", x$Energy, "\n")
#     cat("Sims Energy: M =", x$SimEnergy_M, ", SD =", x$SimEnergy_SD, "\n")
#     cat("Sims with this or higher Energy:", x$Sims_with_higher_energy, "%\n")
#     cat("Number of simulations used:", x$n_sims, "\n")
#   }
#   
#   # Print FFT Test Output
#   if (!is.null(x$Amplitude_sum)) {
#     cat("COE: FFT Amplitude Sum test\n")
#     cat("Sum of Amplitudes:", x$Amplitude_sum, "\n")
#     cat("Sims Amplitude Sums: M =", x$Sim_Ampsum_M, ", SD =", x$Sim_Ampsum_SD, "\n")
#     cat("Simulations with this or higher amplitude sum:", x$Sims_with_higher_Ampsum, "%\n")
#     if (!is.null(x$FFTComparison)) {
#       cat("Frequencies above 95% of Simulations:", x$Frequencies_above_95 * 100, "%\n")
#     }
#     cat("Number of simulations used:", x$n_sims, "\n")
#   }
#   
#   invisible(x)
# }

# Custom print method for 'coe' objects
#' @export
#' @method print coe
print.coe <- function(x, ...) {
  cat("\nChange of Evidence Results \n")
  for (name in names(x)) {
    cat(name, ":", x[[name]], "\n")
  }
  cat("\n")
}

