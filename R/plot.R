# Plot Functions

#' Plot Random Walk
#'
#' This function allows to plot random walks.
#'
#' The Random Walk can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
#' A random walk typically goes +1 for every correct choice and -1 for every incorrect one.
#' Depending on the amount of simulations, drawing might take a while. It might be wise to chose a smaller simulation set for this purpose.
#'
#' @param data A vector containing the random walk to be drawn or a list containing multiple vectors.
#' @param sims.df A dataframe containing simulations, including column "simid" and "index". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column of the simulation dataframe to compare to.
#' @param color A color in which the Random Walk will be drawn.
#' @param coordy A vector containing the minimum and maximum value of the y-coordinates to be drawn.
#' @param mu If Random Walk is of summed up bits, indicate the expected sum per step.
#' @examples
#' p.rw <- plotrw(tbl$rw)
#' p.rw
#'
#' plotrw(tbl$rw, sims.df = sims, sims.df.col = "rw", coordy = c(-50,50))
#' 
#' plotrw(list(exp = exp$rw, con = con$rw), mu = 50)
#'
#' sims1000 <- subset(sims, simid <= 1000)
#' plotrw(tbl, sims.df = sims1000)
#' @export

# Plot Random Walk
plotrw <- function(data, sims.df = NULL, sims.df.col = "rw", color = "black", coordy = c(-absolutemax,absolutemax), mu = NULL){
  library(ggplot2)
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)
  
  # Show legend and decide length for pparabel
  if(is.list(data)){
    show.legend <- "bottom"
    nmax <- max(lengths(data))
    absolutemax <- max(c(max(unlist(data)),abs(min(unlist(data)))))
    
  } else{
    show.legend <- "none"
    nmax <- length(data)
    absolutemax <- max(c(max(data),abs(min(data))))
  }
  
  # Add data point "0" at the beginning
  data <- c(0,data)
  
  # Data for p-parabel
  p.s <- data.frame(n = 0:nmax)
  if(is.null(mu)){
    p.s$p.up <- 1.96*(sqrt(as.numeric(rownames(p.s))))
    p.s$p.dn <- -1.96*(sqrt(as.numeric(rownames(p.s))))
  } else{
    p.s$p.up <- 1.96*(sqrt(as.numeric(rownames(p.s))*2*mu))/2
    p.s$p.dn <- -1.96*(sqrt(as.numeric(rownames(p.s))*2*mu))/2
  }
  xrow <- as.numeric(row.names(p.s))
  
  if(!is.null(sims.df)) print("Depending on the amount of simulations to be drawn, this might take a while!")
  
  p <- ggplot2::ggplot()+
    ggplot2::geom_line(data=p.s, aes(x=xrow, y=p.up), color = "grey60", linetype="dotted", size=1)+
    ggplot2::geom_line(data=p.s, aes(x=xrow, y=p.dn), color = "grey60", linetype="dotted", size=1)
  
  if (!is.null(sims.df)){
    p <- p + ggplot2::geom_line(data=sims.df, aes(x=index, y=.data[[sims.df.col]], group=simid), color=greycol)
  }
  
  if(is.list(data)){
    df <- NULL
    for(i in 1:length(data)){
      ydat <- data[[i]]
      if(!is.null(names(data))){
        ndf <- data.frame(element=as.factor(names(data)[i]),x=1:length(ydat),y=ydat)
      } else{
        ndf <- data.frame(element=paste0("data ",i),x=1:length(ydat),y=ydat)
      }
      df <- rbind(df,ndf)
    }
    p <- p + ggplot2::geom_line(data=df, aes(x=x, y=y, color=element), size=1)+
      ggplot2::scale_color_brewer("Data", type="qualitative", palette="Set1")
  } else {
    p <- p + ggplot2::geom_line(data=as.data.frame(data), aes(x=xrow, y=data), color=color, size=1)
  }
  
  p + ggplot2::geom_hline(yintercept = 0, linetype="dashed", color="grey60", size=1)+
    ggplot2::labs(x="Trials", y = "Random Walk")+
    ggplot2::scale_x_continuous(expand = c(0,0))+
    ggplot2::coord_cartesian(ylim = coordy)+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(legend.position = show.legend)
  
}


#' Plot Sequential Bayesian Analysis
#'
#' This function allows to plot a sequential Bayesian analysis.
#'
#' The BF analysis can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
#' BFs above 1 indicate evidence towards H1, BFs below 1 indicate evidence towards H0.
#' Depending on the amount of simulations, drawing might take a while. It might be wise to chose a smaller simulation set for this purpose.
#'
#' @param data A vector containing sequential Bayes Factors or a list containing multiple vectors.
#' @param sims.df A dataframe containing simulations, including column "simid". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column of the simulation dataframe to compare to.
#' @param color A color in which the Seq BF-function will be drawn.
#' @param coordy A vector containing the minimum and maximum value of the y-coordinates to be drawn.
#' @param label.x A character that is used as label for the x-axis ("N" for sum scores, "Trials" for binomial data).
#' @examples
#' p.bf <- plotbf(tbl$bf)
#' p.bf
#'
#' plotbf(list(test1=bf1,test2=bf2))
#' 
#' plotbf(tbl$bf, sims.df = sims)
#'
#' sims1000 <- subset(sims, simid <= 1000)
#' plotbf(tbl$bf, sims.df = sims1000)
#' @export

# Plot Sequential BF
plotbf <- function(data, sims.df = NULL, sims.df.col = "bf", color = "black", coordy = NULL, label.x = "N"){
  library(ggplot2)
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)
  
  # Set y coordinates
  if(is.null(coordy)){
    if(is.list(data)){
      coordy <- c(Reduce(min,data),2*Reduce(max,data))
    }
    else {
      coordy <- c(min(data),2*max(data))
    }
  }
  
  #Set minimum coordinates to 1/10 and 10
  if(coordy[1] > 1/10) coordy[1] <- 1/10
  if(coordy[2] < 10) coordy[2] <- 10
  
  # Show legend?
  if(is.list(data)) show.legend <- "bottom"
  else show.legend <- "none"
  
  # Specify Text annotations
  annotation <- c("Extreme ~H[1]","paste(\"Very Strong\",~H[1])","Strong ~H[1]","Moderate ~H[1]","Anecdotal ~H[1]","Anecdotal ~H[0]","Moderate ~H[0]","Strong ~H[0]","paste(\"Very Strong\",~H[0])","Extreme ~H[0]")
  annobreaks <- exp(c(mean(c(log(300),log(100))),mean(c(log(100),log(30))),mean(c(log(30),log(10))),mean(c(log(10),log(3))),mean(c(log(3),log(1))),mean(c(log(1),log(1/3))),mean(c(log(1/3),log(1/10))),mean(c(log(1/10),log(1/30))),mean(c(log(1/30),log(1/100))),mean(c(log(1/100),log(1/300)))))
  
  if(coordy[2] > 96 && coordy[2] < 125){
    annotation <- annotation[-1]
    annobreaks <- annobreaks[-1]
  }
  if(coordy[1] < 1/96 && coordy[1] > 1/125){
    annotation <- head(annotation,-1)
    annobreaks <- head(annobreaks,-1)
  }
  if(coordy[2] > 32 && coordy[2] < 41){
    annotation <- annotation[-c(1:2)]
    annobreaks <- annobreaks[-c(1:2)]
  }
  if(coordy[1] < 1/32 && coordy[1] > 1/41){
    annotation <- head(annotation,-2)
    annobreaks <- head(annobreaks,-2)
  }
  if(coordy[2] > 12 && coordy[2] < 15){
    annotation <- annotation[-c(1:3)]
    annobreaks <- annobreaks[-c(1:3)]
  }
  if(coordy[1] < 1/12 && coordy[1] > 1/15){
    annotation <- head(annotation,-3)
    annobreaks <- head(annobreaks,-3)
  }
  
  #Scale y-Axis
  if(coordy[2] <= 1000 && coordy[1] >= 1/1000){
    breaks = c(1000,300,100,30,10,3,1,1/3,1/10,1/30,1/100,1/300,1/1000)
    labels = c("1,000","300","100","30","10","3","1","1/3","1/10","1/30","1/100","1/300","1/1,000")
  } else if(coordy[2] <= 1000000 && coordy[1] >= 1/1000000){
    breaks = c(1000000,100000,10000,1000,100,10,1,1/10,1/100,1/1000,1/10000,1/100000,1/1000000)
    labels = c("1,000,000","100,000","10,000","1,000","100","10","1","1/10","1/100","1/1,000","1/10,000","1/100,000","1/1,000,000")
    #annotation = annotation[c(1,3,5,6,8,10)]
    #annobreaks = exp(c(mean(c(log(300),log(100))),mean(c(log(100),log(30))),mean(c(log(30),log(10))),mean(c(log(10),log(3))),mean(c(log(3),log(1))),mean(c(log(1),log(1/3))),mean(c(log(1/3),log(1/10))),mean(c(log(1/10),log(1/30))),mean(c(log(1/30),log(1/100))),mean(c(log(1/100),log(1/300)))))
  } else {
    breaks = c(1e+12,1e+10,1e+8,1e+6,1e+4,1e+2,1,1/1e+2,1/1e+4,1/1e+6,1/1e+8,1/1e+10,1/1e+12)
    labels = c("1e+12","1e+10","1e+8","1e+6","1e+4","1e+2","1","1/1e+2","1/1e+4","1/1e+6","1/1e+8","1/1e+10","1/1e+12")
    #annotation = annotation[c(1,10)]
    #annobreaks = exp(c(mean(c(log(1000),log(100))),mean(c(log(1/100),log(1/1000)))))
  }
  
  if(!is.null(sims.df)) print("Depending on the amount of simulations to be drawn, this might take a while!")
  
  p <- ggplot2::ggplot()
  if(!is.null(sims.df)) p <- p + ggplot2::geom_line(data=sims.df, aes(x=index, y=.data[[sims.df.col]], group=simid), color=greycol)
  p <- p + ggplot2::geom_hline(yintercept = 1, color='grey60', linetype = 'solid')+
    ggplot2::geom_hline(yintercept = breaks, color='grey60', linetype='dotted')
  if(is.list(data)){
    df <- NULL
    for(i in 1:length(data)){
      ydat <- data[[i]]
      if(!is.null(names(data))){
        ndf <- data.frame(element=as.factor(names(data)[i]),x=1:length(ydat),y=ydat)
      } else{
        ndf <- data.frame(element=paste0("data ",i),x=1:length(ydat),y=ydat)
      }
      df <- rbind(df,ndf)
    }
    p <- p + ggplot2::geom_line(data=df, aes(x=x, y=y, color=element), size=1)+
      ggplot2::scale_color_brewer("Data", type="qualitative", palette="Set1")
      if(coordy[2] <= 1000 && coordy[1] >= 1/1000) p <- p + ggplot2::annotate("text", x=max(df$x)*1.2, y=annobreaks, label=annotation, hjust=1, parse = TRUE)
    
  } else {
    p <- p + ggplot2::geom_line(data=as.data.frame(data), aes(x=as.numeric(1:length(data)), y=data), color=color, size=1)
      if(coordy[2] <= 1000 && coordy[1] >= 1/1000) p <- p + ggplot2::annotate("text", x=length(data)*1.2, y=annobreaks, label=annotation, hjust=1, parse = TRUE)
  }
    p + ggplot2::labs(x=label.x, y = "Evidence (BF)")+
    ggplot2::scale_y_log10(breaks = breaks, labels = labels)+
    ggplot2::coord_cartesian(ylim = coordy)+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(legend.position = show.legend, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
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
#' @param color A color in which the FFT will be drawn.
#' @param coordy A vector containing the minimum and maximum value of the y-coordinates to be drawn.
#' @examples
#' p.fftbf <- plotfft(tblFFT$density.bf, color = "blue")
#' p.fftbf
#'
#' p.fftrw <- plotfft(tblFFT$density.rw, sims.df = sims, sims.df.col = "density.rw")
#' p.fftrw
#' @export

# Plot FFT
# Data for 95-CI ribbon FFT
plotfft <- function(data, sims.df = NULL, sims.df.col = "density.bf", n.hz = 50, color = "black", coordy = c(0,secondhighestval)){
  secondhighestval <- sort(data,partial=length(data)-1)[length(data)-1]
  library(ggplot2)
  #cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)
  
  if(length(data) < n.hz) n.hz <- length(data)
  
  p <- ggplot2::ggplot()
  
  if(!is.null(sims.df)){
    
    simci.fft <- numeric(n.hz)
    for (sindex in 1:n.hz){
      simci.fft[sindex] <- sort(sims.df[sims.df$index == sindex,sims.df.col])[max(sims.df$simid)*0.95]
    }
    
    secondhighestval <- max(secondhighestval,simci.fft[1])

    p <- p+
      ggplot2::geom_line(data=sims.df, aes(x=index, y=.data[[sims.df.col]], group=simid), color=greycol)+
      ggplot2::geom_line(data=as.data.frame(simci.fft), aes(x=as.numeric(1:n.hz), y=simci.fft), linetype="dotted", size=1)
  }
  # Plot FFT
  p+
    ggplot2::geom_line(data=as.data.frame(data), aes(x=as.numeric(1:length(data)), y=data), color=color, size=1)+
    ggplot2::labs(title=paste0("Fast Fourier Transform"), x="Frequency (No of Cycles)", y = "Amplitude")+
    ggplot2::coord_cartesian(xlim = c(1,n.hz), ylim = coordy)+
    ggplot2::scale_x_continuous(breaks = seq(0, n.hz, by = (n.hz/(n.hz/2))), expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0, 50, by = 1))+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(legend.position = "none")
  
}