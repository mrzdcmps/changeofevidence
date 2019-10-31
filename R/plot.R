# Plot Functions

#' Plot Random Walk
#'
#' This function allows to plot random walks.
#'
#' The Random Walk can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
#' A random walk typically goes +1 for every correct choice and -1 for every incorrect one.
#' Depending on the amount of simulations, drawing might take a while. It might be wise to chose a smaller simulation set for this purpose.
#'
#' @param data A vector containing the random walk to be drawn.
#' @param sims.df A dataframe containing simulations, including column "simid" and "index". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column of the simulation dataframe to compare to.
#' @examples
#' p.rw <- plotrw(tbl$rw)
#' p.rw
#'
#' plotrw(tbl$rw, sims.df = sims, sims.df.col = "rw")
#'
#' sims1000 <- subset(sims, simid <= 1000)
#' plotrw(tbl, sims.df = sims1000)
#' @export

# Plot Random Walk
plotrw <- function(data, sims.df = NULL, sims.df.col = "rw", color = 1){
  library(ggplot2)
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)

  # Data for p-parabel
  p.s <- data.frame(n = 1:length(data))
  p.s$p.up <- 1.96*(sqrt(as.numeric(rownames(p.s))))
  p.s$p.dn <- -1.96*(sqrt(as.numeric(rownames(p.s))))

  xrow <- as.numeric(row.names(p.s))

  absolutemax <- max(c(max(data),abs(min(data))))

  if(!is.null(sims.df)) print("Depending on the amount of simulations to be drawn, this might take a while!")

  p <- ggplot2::ggplot()
  if (!is.null(sims.df)){
    p <- p + ggplot2::geom_line(data=sims.df, aes(x=index, y=sims.df[[sims.df.col]], group=simid), color=greycol)
  }
  p + ggplot2::geom_line(data=p.s, aes(x=xrow, y=p.up), color = "black", linetype="dotted", size=1)+
    ggplot2::geom_line(data=p.s, aes(x=xrow, y=p.dn), color = "black", linetype="dotted", size=1)+
    ggplot2::geom_line(data=as.data.frame(data), aes(x=xrow, y=data), color=cbPalette[color], size=1)+
    ggplot2::geom_hline(yintercept = 0, linetype="dashed", color="grey60", size=1)+
    ggplot2::labs(x="Trials", y = "Random Walk")+
    ggplot2::coord_cartesian(ylim = c(-absolutemax,absolutemax))+
    ggplot2::scale_x_continuous(expand = c(0,0))+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(legend.position = 'none')

  #ggsave(paste0("randomwalk.png"), width = 9, height = 7, dpi = 300, limitsize = TRUE)
}


#' Plot Sequential Bayesian Analysis
#'
#' This function allows to plot a sequential Bayesian analysis.
#'
#' The BF analysis can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
#' BFs above 1 indicate evidence towards H1, BFs below 1 indicate evidence towards H0.
#' Depending on the amount of simulations, drawing might take a while. It might be wise to chose a smaller simulation set for this purpose.
#'
#' @param df A vector containing sequential Bayes Factors.
#' @param sims.df A dataframe containing simulations, including column "simid". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column of the simulation dataframe to compare to.
#' @examples
#' p.bf <- plotbf(tbl$bf)
#' p.bf
#'
#' plotbf(tbl$bf, sims.df = sims)
#'
#' sims1000 <- subset(sims, simid <= 1000)
#' plotbf(tbl$bf, sims.df = sims1000)
#' @export

# Plot Sequential BF
plotbf <- function(data, sims.df = NULL, sims.df.col = "bf", color = 2){
  library(ggplot2)
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)

  if(!is.null(sims.df)) print("Depending on the amount of simulations to be drawn, this might take a while!")

  p <- ggplot2::ggplot()
  if (!is.null(sims.df)){
    p <- p + ggplot2::geom_line(data=sims.df, aes(x=index, y=sims.df[[sims.df.col]], group=simid), color=greycol)
  }
  p + ggplot2::geom_hline(yintercept = 1, color='grey60', linetype = 'solid')+
    ggplot2::geom_hline(yintercept = c(1000,300,100,30,10,3,1/3,1/10,1/30), color='grey60', linetype='dashed')+
    ggplot2::geom_line(data=as.data.frame(data), aes(x=as.numeric(1:length(data)), y=data), color=cbPalette[color], size=1)+
    ggplot2::labs(x="Trials", y = "Evidence (BF)")+
    ggplot2::scale_y_continuous(trans='log10', breaks = c(1000,300,100,30,10,3,1,1/3,1/10,1/30), labels = c("1000","300","100","30","10","3","1","1/3","1/10","1/30"))+
    ggplot2::scale_x_continuous(expand = c(0,0))+
    ggplot2::coord_cartesian(ylim = c(min(data),2*max(data)))+
    ggplot2::theme_classic(base_size = 14)+
    ggplot2::theme(legend.position = "none")

  #ggsave(paste0("bf.png"), width = 9, height = 7, dpi = 300, limitsize = TRUE)
}


#' Plot Fast Fourier Transform
#'
#' This function allows to plot a Fast Fourier Transform.
#'
#' The FFT can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
#' An experimental 95% Confidence Interval is calculated. Values above the black dotted line indicate Top5-amplitudes.
#' Use the full simulation set to correcly identify the 5% border line.
#'
#' @param data A vector containing Fourier transformed spectral densities.
#' @param sims.df A dataframe containing simulations, including column "simid". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column of the simulation dataframe to compare to.
#' @param n.hz The amount of frequencies to be displayed.
#' @examples
#' p.fftbf <- plotfft(tblFFT$density.bf, color = 2)
#' p.fftbf
#'
#' p.fftrw <- plotfft(tblFFT$density.rw, sims.df = sims, sims.df.col = "density.rw", color = 1)
#' p.fftrw
#' @export

# Plot FFT
# Data for 95-CI ribbon FFT
plotfft <- function(data, sims.df = NULL, sims.df.col = "density.bf", n.hz = 50, color = 3){
  library(ggplot2)
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)

  p <- ggplot2::ggplot()

  if(!is.null(sims.df)){

    simci.fft <- numeric(n.hz)
    for (sindex in 1:n.hz){
      simci.fft[sindex] <- sort(sims.df[sims.df$index == sindex,sims.df.col])[max(sims.df$simid)*0.95]
    }
    p <- p+
      ggplot2::geom_line(data=sims.df, aes(x=index, y=sims.df[[sims.df.col]], group=simid), color=greycol)+
      ggplot2::geom_line(data=as.data.frame(simci.fft), aes(x=as.numeric(1:n.hz), y=simci.fft), linetype="dotted", size=1)
  }
  # Plot FFT
  p+
    ggplot2::geom_line(data=as.data.frame(data), aes(x=as.numeric(1:length(data)), y=data), color=cbPalette[color], size=1)+
    ggplot2::labs(title=paste0("Fast Fourier Transform"), x="Frequency (No of Cycles)", y = "Amplitude")+
    ggplot2::coord_cartesian(xlim = c(1,n.hz), ylim = c(0,sort(data,partial=length(data)-1)[length(data)-1]))+
    ggplot2::scale_x_continuous(breaks = seq(0, n.hz, by = (n.hz/(n.hz/2))), expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0, 50, by = 1))+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(legend.position = "none")

  #ggsave(paste0("fft_rw.png"), width = 9, height = 7, dpi = 300, limitsize = TRUE)
}
