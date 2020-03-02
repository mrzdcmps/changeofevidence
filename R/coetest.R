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
  u.nsims <- max(sims.df$simid)
  sim.maxbf <- tapply(sims.df$bf, sims.df$simid, max, na.rm=TRUE)
  cat("Highest BF:",max(data),"( at N =",which.max(data),") \n")
  cat("Percentage of Sims with higher BFs:",(sum(sim.maxbf > max(data))/u.nsims)*100," \n")
  maxbf.out <- list("MaxBF" = max(data), "MaxBF (N)" = which.max(data), "Sims with higher BFs (%)" = (sum(sim.maxbf > max(data))/u.nsims)*100)
  return(maxbf.out)
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
  u.nsims <- max(sims.df$simid)
  sim.energy <- numeric(length = u.nsims)
  
  cat("Calculating Energy of sims... \n")
  pb = txtProgressBar(min = 0, max = u.nsims, initial = 0, style = 3)
  for (sid in 1:u.nsims){
    sim.energy.data <- subset(sims.df, sims.df$simid == sid)
    sim.energy[sid] <- pracma::trapz(as.numeric(rownames(sim.energy.data)), sim.energy.data$bf)-pracma::trapz(as.numeric(rownames(sim.energy.data)), rep(1, nrow(sim.energy.data)))
    setTxtProgressBar(pb,sid)
  }
  
  real.energy <- pracma::trapz(as.numeric(1:length(data)), data)-pracma::trapz(as.numeric(1:length(data)), rep(1, length(data)))
  
  cat("\n Energy of BF of data: ",real.energy,"\n")
  cat("Sims Energy: M =",mean(sim.energy),", SD =",sd(sim.energy),"\n")
  cat("Percentage of Sims with higher Energy:",(sum(sim.energy > real.energy)/u.nsims)*100," \n")
  energybf.out <- list("Energy" = real.energy, "Simenergy (M)" = mean(sim.energy), "Simenergy (SD)" = sd(sim.energy), "Sims with more energy (%)" = (sum(sim.energy > real.energy)/u.nsims)*100)
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
  T <- 1/L # Tastrate
  Y <- fft(data) # Fast Fourier Transformation
  P <- abs(Y/L)[1:(L/2)]
}


#' Frequency Analysis
#'
#' This function analyses FFTs and compares the amplitudes of all frequencies to those of simulations.
#'
#' For each Frequency this function counts how many simulations show a higher amplitude.
#' If no more than 5% of simulations are above the experimental value, it is added to the "Top5"-list.
#' The proportion of Top5-Frequencies indicates the pronouncedness of oscillatory elements in the data.
#'
#' @param data A vector containing Fourier transformed (spectral density) data.
#' @param sims.df A dataframe containing simulations, including columns "index" and "simid".
#' @param sims.df.col The column of the simulation dataframe that contains the comparison data.
#' @return A list containing the proportions of sims with a lower amplitude than each frequency.
#' @examples
#' r.fftbf <- ffttest(tblFFT$density.bf)
#' r.fftrw <- ffttest(tblFFT$density.rw, sims.df = newsims, sims.df.col = "density.rw")
#' @export

# Check for high Amplitudes
ffttest <- function(data, sims.df = sims, sims.df.col = "density.bf"){
  cat(">> FREQUENCY ANALYSIS << \n")
  list.HB <- list()
  pb = txtProgressBar(min = 0, max = length(data), initial = 0, style = 3)
  for (H in 1:length(data)){
    
    CHz <- sims.df[sims.df$index==H,]
    #if (((sum(CHz[[sims.df.col]] > data[H]))/max(sims.df$simid))<0.05) {
      #if (H < 50) {
      #  cat(H, "Hz: ",((sum(CHz$density > data[H]))/u.nsims)*100,"% \n")
      #}
      list.HB[[length(list.HB)+1]] = data.frame(H = H, LowerSims = (sum(CHz[[sims.df.col]] < data[H]))/max(sims.df$simid))
    #}
    setTxtProgressBar(pb,H)
  }
  cat("\nNumber of Frequencies: ",length(data),"\n")
  cat("Number of Frequencies above 95% of Simulations:",sum(list.HB$LowerSims > 0.95),"\n")
  cat("Percentage of Frequencies above 95% of Simulations:",(sum(list.HB$LowerSims > 0.95)/length(data))*100,"% \n")
  #fftbf.out <- list("Top5 Frequencies (#)" = length(list.HB), "Top5 Frequencies (%)"=(length(list.HB)/length(data))*100, "Top5 Frequencies" = list.HB)
  fftbf.out <- list.HB
  return(fftbf.out)
}