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


### Helper-Function: Count Frequencies
.fftcount <- function(data, sims.df = sims, sims.df.col = "density.bf"){
  
  CHz <- sims.df[sims.df$index==as.numeric(data[2]),]
  output <- data.frame(H = as.numeric(data[2]), LowerSims = (sum(CHz[[sims.df.col]] < as.numeric(data[1])))/max(sims.df$simid))
  output
  
}

#' Frequency Analysis Test
#'
#' This function analyses FFTs and compares the amplitudes of all frequencies to those of simulations.
#'
#' For each Frequency this function counts how many simulations show a higher amplitude.
#' If no more than 5% of simulations are above the experimental value, it is considered a "Top5-Frequency".
#' The proportion of Top5-Frequencies indicates the pronouncedness of oscillatory elements in the data.
#'
#' @param data A vector containing Fourier transformed (spectral density) data.
#' @param sims.df A dataframe containing simulations, including columns "index" and "simid".
#' @param sims.df.col The column of the simulation dataframe that contains the comparison data.
#' @return A dataframe containing all frequencies and the proportion of simulations with a lower amplitude.
#' @examples
#' r.fftbf <- ffttest(tblFFT$density.bf)
#' r.fftrw <- ffttest(tblFFT$density.rw, sims.df = newsims, sims.df.col = "density.rw")
#' @export

ffttest <- function(data, sims.df = sims, sims.df.col = "density.bf"){
  if(!is.numeric(sims.df[[sims.df.col]])) stop("Wrong sims data. Does sims.df.col exist?")
  cat(">> FREQUENCY ANALYSIS << \n")
  data.df <- data.frame(data = data, H = seq_along(data))
  
  list.HB <- pbapply::pbapply(data.df,1,.fftcount,sims.df = sims.df, sims.df.col = sims.df.col)
  list.HB <- dplyr::bind_rows(list.HB)
  
  cat("\nNumber of Frequencies: ",length(data),"\n")
  cat("Number of Frequencies above 95% of Simulations:",sum(list.HB$LowerSims > 0.95),"\n")
  cat("Percentage of Frequencies above 95% of Simulations:",(sum(list.HB$LowerSims > 0.95)/length(data))*100,"% \n")
  #fftbf.out <- list("Top5 Frequencies (#)" = length(list.HB), "Top5 Frequencies (%)"=(length(list.HB)/length(data))*100, "Top5 Frequencies" = list.HB)
  #fftbf.out <- list.HB
  return(list.HB)
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
#' @export


# Count likelihood and distribution of Top5 occurrences
fftlikelihood <- function(df, proportion = 100, sims.df = sims, sims.df.col = "density.bf"){
  require(foreach)
  require(doParallel)
  
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  cat(">> FREQUENCY ANALYSIS LIKELIHOOD << \n")
  cat("This test runs in parallel. See log.txt for status updates!")
  
  likelihoodlist <- foreach(i=1:(max(sims.df$simid)*(proportion/100)), .combine=c) %dopar% {
    sink("log.txt", append=TRUE)  
    cat(paste(Sys.time(),"Starting iteration",i,"of",max(sims.df$simid)*(proportion/100),"\n"))
    sink()
    tmpdat.r <- subset(sims.df, simid == i)
    tmptest <- changeofevidence::ffttest(tmpdat.r$density.bf, sims.df)
    likeresult <- sum(tmptest$LowerSims > 0.95)
    likeresult
  }
  stopCluster(cl)
  
  cat("\nLikelihood for",sum(df$LowerSims > 0.95),"(=",(sum(df$LowerSims > 0.95)/nrow(df))*100,"%) or more Top5-Frequencies is estimated",(sum(likelihoodlist >= sum(df$LowerSims > 0.95))/length(likelihoodlist))*100,"% \n")
  
  likelihood.out <- list("Likelihood" = sum(likelihoodlist >= sum(df$LowerSims > 0.95))/length(likelihoodlist), "Top5 of Sims" = likelihoodlist)
  return(likelihood.out)
}
