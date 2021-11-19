#' Maximum BF Analysis
#'
#' This function evaluates the likelihood of the highest reached Bayes Factor.
#'
#' The highest reached (Maximum) BF of the experimental data is compared to those of all simulations.
#' An unusual peak can be assumed if less than 5\% of simulations show a higher BF at any time.
#'
#' @param data A vector containing Bayes Factors.
#' @param sims.df A dataframe containing simulations, including columns "simid" and "bf".
#' @return A list containing the Maximum BF of the experimental data, its position in the temporal order of data points, and the proportion of simulations with higher BFs at any time.
#' @examples
#' r.maxbf <- maxbf(tbl$bf)
#' r.maxbf <- maxbf(tbl$bf, sims.df = newsims)
#' @export

# Maximum BF analysis
maxbf <- function(data, sims.df=sims){
  cat(">> MAXIMUM BF << \n")
  u.nsims <- length(unique(sims.df$simid))
  sim.maxbf <- tapply(sims.df$bf, sims.df$simid, max, na.rm=TRUE)
  cat("Highest BF:",max(data),"( at N =",which.max(data),") \n")
  cat("Sims with this or a higher BFs:",(sum(sim.maxbf >= max(data))/u.nsims)*100,"% \n")
  maxbf.out <- list("MaxBF" = max(data), "MaxBF (N)" = which.max(data), "Sims with this or higher BFs" = (sum(sim.maxbf >= max(data))/u.nsims))
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
#' @param data A vector containing Bayes Factors.
#' @param sims.df A dataframe containing simulations, including columns "simid" and "bf".
#' @return A list containing the Energy of the experimental data, the mean and SD of the energy values of all simulations, and the proportion of simulations with a higher energy than the experimental data.
#' @examples
#' r.energybf <- energybf(tbl$bf)
#' r.energybf <- energybf(tbl$bf, sims.df = newsims)
#' @export

# Energy of BF
energybf <- function(data, sims.df=sims){
  cat(">> BF ENERGY << \n")
  
  simids <- sims.df[seq(1, nrow(sims.df), max(sims.df$index)),]$simid
  u.nsims <- length(simids)
  sim.energy <- numeric(length = u.nsims)
  
  if(!is.numeric(data)) stop("Data must be a numeric vector!")
  if(length(data) != nrow(sims.df[sims.df$simid == simids[1],])) stop("Data is not the same length as simulations!")
  
  cat("Calculating Energy of sims... \n")
  
  nullenergy <- length(data)-1
  
  sim.energy <- pbapply::pbtapply(sims.df$bf, sims.df$simid, .energycount, nullenergy)
  
  real.energy <- pracma::trapz(as.numeric(1:length(data)), data)-nullenergy
  
  cat("\n Energy of BF of data: ",real.energy,"\n")
  cat("Sims Energy: M =",mean(sim.energy),", SD =",sd(sim.energy),"\n")
  cat("Sims with this or higher Energy:",(sum(sim.energy >= real.energy)/u.nsims)*100,"% \n")
  energybf.out <- list("Energy" = real.energy, "Simenergy (M)" = mean(sim.energy), "Simenergy (SD)" = sd(sim.energy), "Sims with this or higher energy" = (sum(sim.energy >= real.energy)/u.nsims))
  return(energybf.out)
}


#' Create a Fast Fourier transformed vector
#'
#' This function converts input signals into its Fast Fourier transform (FFT).
#'
#' The FFT shows the spectral density of the input. It can be understood as a harmonic analysis and indicates the strength of osciallative components comprising the signal.
#' It can therefore be used to assess the meaning of oscillation as a characteristic. N/2 Frequencies are used as sampling rate.
#'
#' @param data A vector to be transformed.
#' @examples
#' density.bf <- fftcreate(tbl$bf)
#' tblFFT <- data.frame(density.bf = fftcreate(tbl$bf), density.rw = fftcreate(tbl$rw))
#' @export

# Create FFT
fftcreate <- function(data){
  L <- length(data) # Länge
  L2 <- ceiling(L/2) # Consider only first half, round up
  T <- 1/L # Tastrate
  Y <- fft(data) # Fast Fourier Transformation
  P <- abs(Y/L)[1:L2]
}

### Helper-Function: Create FFT without cutting frequencies
.fftcreatefull <- function(data){
  L <- length(data) # Länge
  T <- 1/L # Tastrate
  Y <- fft(data) # Fast Fourier Transformation
  P <- abs(Y/L)
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
#' @param data A vector containing Fourier transformed (spectral density) data (use 'fftcreate' function).
#' @param sims.df A dataframe containing simulations, including columns "index" and "simid".
#' @param sims.df.col The column of the simulation dataframe that contains the comparison data.
#' @param top5 Logical. If set to TRUE, function will additionally return the Top5-Frequency method. For each frequency the function counts how many simulations show a higher amplitude. If no more than 5 percent of simulations are above the experimental value, it is considered a "Top5-Frequency". The proportion of Top5-Frequencies indicates the pronouncedness of oscillatory elements in the data.
#' @return A list containing a dataframe of all frequencies and the proportion of simulations with a lower amplitude, and information on the sum of amplitudes.
#' @examples
#' r.fftbf <- ffttest(tblFFT$density.bf)
#' r.fftrw <- ffttest(tblFFT$density.rw, sims.df = newsims, sims.df.col = "density.rw")
#' @export

ffttest <- function(data, sims.df = sims, sims.df.col = "density.bf", top5 = FALSE){
  if(var(data[1:3]) == 0) stop("It seems like you specified a vector containing BFs. Please use fftcreate(bf) to Fourier transform first.")
  if(!is.numeric(sims.df[[sims.df.col]])) stop("Wrong sims data. Does sims.df.col exist?")
  if(ceiling(max(sims.df$index/2)) != length(data)) stop("Lengths of FFT data and sims do not match! Consider using simredo()")
  cat(">> FREQUENCY ANALYSIS << \n")
  u.nsims <- length(unique(sims.df$simid))
  data.df <- data.frame(data = data, H = seq_along(data))
  ampsum <- sum(data.df$data)
  
  sims.df <- sims.df[sims.df$index <= length(data),]
  
  if(top5 == TRUE){
    list.HB <- pbapply::pbapply(data.df,1,.fftcount,sims.df = sims.df, sims.df.col = sims.df.col)
    list.HB <- dplyr::bind_rows(list.HB)
  }
  
  simampsum <- tapply(sims.df[[sims.df.col]], sims.df$simid, sum)
  
  cat("\nNumber of Frequencies: ",length(data),"\n")
  if(top5 == TRUE){
  cat("Number of Frequencies above 95% of Simulations:",sum(list.HB$LowerSims > 0.95),"(",(sum(list.HB$LowerSims > 0.95)/length(data))*100,"% )\n")
  cat("-------\n")
  }
  cat("Sum of Amplitudes:",ampsum,"\n")
  cat("Sims Amplitude Sums: M =",mean(simampsum),"SD =",sd(simampsum),"\n")
  cat("Simulations with this or higher amplitude sum:",(sum(ampsum <= simampsum)/u.nsims)*100,"% \n")
  if(top5 == TRUE){
  fftbf.out <- list("FFTComparison" = list.HB, "Frequencies above 95% of Sims" = sum(list.HB$LowerSims > 0.95)/length(data), "Amplitude sum" = ampsum, "Sims with higher Ampsum" = sum(ampsum < simampsum)/u.nsims, "Sim Ampsum (M)" = mean(simampsum), "Sim Ampsum (SD)" = sd(simampsum))
  } else {
  fftbf.out <- list("Amplitude sum" = ampsum, "Sims with this or higher Ampsum" = sum(ampsum <= simampsum)/u.nsims, "Sim Ampsum (M)" = mean(simampsum), "Sim Ampsum (SD)" = sd(simampsum))
  }  
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
