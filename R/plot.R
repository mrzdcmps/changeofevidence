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
#' @param p Indicate the probability of success for a single step.
#' @param mu If bits are summed up for a single step, this is the mean of the sum.
#' @examples
#' p.rw <- plotrw(tbl$rw)
#' p.rw
#'
#' plotrw(tbl$rw, sims.df = sims, sims.df.col = "rw", coordy = c(-50,50))
#' 
#' plotrw(list(exp = exp$rw, con = con$rw), mu = 5)
#' 
#' df$rw <- cumsum(ifelse(df$correct == TRUE, 1, -1))
#' plotrw(df$rw, p = 0.2) # 1 in 5 probability of success
#'
#' sims1000 <- subset(sims, simid <= 1000)
#' plotrw(tbl, sims.df = sims1000)
#' @export

# Plot Random Walk
plotrw <- function(data, sims.df = NULL, sims.df.col = "rw", color = "black", coordy = c(-absolutemax,absolutemax), mu = NULL, p = 0.5){
  library(ggplot2)
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)
  
  # Adjust the random walk if p != 0.5
  adjust_rw <- function(rw, p) {
    if (p != 0.5) {
      # Correct adjustment: Expected value at each step is 2p-1
      steps <- rw * (2*p - 1)
      return(cumsum(steps))
    }
    return(rw)
  }
  
  # Show legend and decide length for p-parabel
  if (is.list(data)) {
    show.legend <- "bottom"
    nmax <- max(lengths(data))
    absolutemax <- max(c(max(unlist(data)), abs(min(unlist(data)))))
    for (i in 1:length(data)) {
      data[[i]] <- adjust_rw(c(0, data[[i]]), p)
    }
  } else {
    show.legend <- "none"
    nmax <- length(data)
    absolutemax <- max(c(max(data), abs(min(data))))
    data <- adjust_rw(c(0, data), p)
  }
  
  # Data for confidence bounds (p-parabel)
  z <- 1.96  # 95% confidence interval
  p.s <- data.frame(n = 0:nmax)
  if(is.null(mu)){
    # Corrected variance calculation for random walk
    # Var(Sn) = 4np(1-p) where Sn is the position after n steps
    p.s$p.up <- (2*p - 1) * p.s$n + z * sqrt(4 * p.s$n * p * (1-p))
    p.s$p.dn <- (2*p - 1) * p.s$n - z * sqrt(4 * p.s$n * p * (1-p))
  } else {
    # If mu is specified, use it for the drift
    p.s$p.up <- mu * p.s$n + z * sqrt(2 * p.s$n * mu)
    p.s$p.dn <- mu * p.s$n - z * sqrt(2 * p.s$n * mu)
  }
  
  xrow <- as.numeric(p.s$n)
  
  if (!is.null(sims.df)) print("Depending on the amount of simulations to be drawn, this might take a while!")
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = p.s, aes(x = xrow, y = p.up), color = "grey60", linetype = "dotted", linewidth = 1) +
    ggplot2::geom_line(data = p.s, aes(x = xrow, y = p.dn), color = "grey60", linetype = "dotted", linewidth = 1)
  
  if (!is.null(sims.df)) {
    p <- p + ggplot2::geom_line(data = sims.df, aes(x = index, y = .data[[sims.df.col]], group = simid), color = greycol)
  }
  
  if (is.list(data)) {
    df <- NULL
    for (i in 1:length(data)) {
      ydat <- data[[i]]
      if (!is.null(names(data))) {
        ndf <- data.frame(element = as.factor(names(data)[i]), x = 1:length(ydat), y = ydat)
      } else {
        ndf <- data.frame(element = paste0("data ", i), x = 1:length(ydat), y = ydat)
      }
      df <- rbind(df, ndf)
    }
    p <- p + ggplot2::geom_line(data = df, aes(x = x - 1, y = y, color = element), linewidth = 1) +
      ggplot2::scale_color_brewer("Data", type = "qualitative", palette = "Set1")
  } else {
    p <- p + ggplot2::geom_line(data = as.data.frame(data), aes(x = xrow, y = data), color = color, linewidth = 1)
  }
  
  p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1) +
    ggplot2::labs(x = "Trials", y = "Random Walk") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = coordy) +
    ggplot2::theme_bw(base_size = 14) +
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
#' @param data A seqbf object or a vector containing sequential Bayes Factors or a list containing multiple vectors.
#' @param sims.df A dataframe containing simulations, including column "simid". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column of the simulation dataframe to compare to.
#' @param color A color in which the Seq BF-function will be drawn.
#' @param coordy A vector containing the minimum and maximum value of the y-coordinates to be drawn.
#' @param label.x A character that overrides the label for the x-axis ("N" for sum scores, "Trials" for binomial data).
#' @examples
#' 
#' plot(seqbf)
#' 
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
  
  if(inherits(data,"seqbf") == TRUE){
    if(data$`test type` == "independent") {
      testtype <- "Independent Samples"
    } else if(data$`test type` == "paired") {
      testtype <- "Paired Samples"
    } else if(data$`test type` == "one-sample") {
      testtype <- "One-Sample"
    } else if(data$`test type` == "binomial") {
      testtype <- "Binomial"
      label.x <- "Trials"
    } else if(data$`test type` == "correlation") {
      testtype <- "Correlation"
    } else {
      testtype <- "Unknown"
    }
    tails <- ifelse(data$alternative == "two.sided", "two-tailed", "one-tailed")
    subtitle <- paste0("BF = ",round(tail(data$BF,n=1),3)," (N = ", sum(data$`sample size`),")")
    caption <- paste0(
      testtype, " test; ",
      tails, "; ",
      data$prior[[1]],"(",data$prior[[2]],", ",data$prior[[3]],")")
    data <- data$BF
    showinfo <- TRUE
  }
  
  library(ggplot2)
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)
  
  # Set y coordinates
  if(is.null(coordy)){
    if(is.list(data)){
      coordy <- c(min(unlist(data), na.rm=T),
                  2*max(unlist(data), na.rm=T))
    }
    else {
      coordy <- c(min(data, na.rm=T),2*max(data, na.rm=T))
    }
  }
  
  #Set minimum coordinates to 1/10 and 10
  if(coordy[1] > 1/10) coordy[1] <- 1/10
  if(coordy[2] < 10) coordy[2] <- 10
  
  # Show legend?
  if(is.list(data)) show.legend <- "bottom"
  else show.legend <- "none"
  
  # Specify Text annotations
  annotationlist <- .annotations(coordy)
  
  annotation <- annotationlist[[1]]
  annobreaks <- annotationlist[[2]]
  
  #Scale y-Axis
  breaks <- annotationlist[[3]]
  labels <- annotationlist[[4]]
  
  
  if(!is.null(sims.df)) print("Depending on the amount of simulations to be drawn, this might take a while!")
  
  # Draw plot
  p <- ggplot2::ggplot()
  
  # Add simulations
  if(!is.null(sims.df)) p <- p + ggplot2::geom_line(data=sims.df, aes(x=index, y=.data[[sims.df.col]], group=simid), color=greycol)
  
  # Add horizontal lines
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
      df <- df[!is.na(df$y),]
    }
    p <- p + ggplot2::geom_line(data=df, aes(x=x, y=y, color=element), linewidth=1)+
      ggplot2::scale_color_brewer("Data", type="qualitative", palette="Set1")
    if(coordy[2] <= 1000 && coordy[1] >= 1/1000) p <- p + ggplot2::annotate("text", x=max(df$x)*1.2, y=annobreaks, label=annotation, hjust=1, parse = TRUE)
    
  } else {
    df <- data.frame(y=data, x=as.numeric(1:length(data)))
    df <- df[!is.na(df$y),]
    
    p <- p + ggplot2::geom_line(data=df, aes(x=x, y=y), color=color, linewidth=1)
    if(coordy[2] <= 1000 && coordy[1] >= 1/1000) p <- p + ggplot2::annotate("text", x=length(data)*1.2, y=annobreaks, label=annotation, hjust=1, parse = TRUE)
  }
  
  if(exists("showinfo")) p <- p + ggplot2::labs(subtitle = subtitle, caption = caption)
  
  p + ggplot2::labs(x=label.x, y = "Evidence (BF)")+
    ggplot2::scale_y_log10(breaks = breaks, labels = labels)+
    ggplot2::coord_cartesian(ylim = coordy)+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(legend.position = show.legend, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
}

#' @export
plot.seqbf <- function(data, ...) {
  suppressWarnings(plotbf(data, ...))
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
  if(inherits(data,"seqbf") == TRUE) data <- changeofevidence::fftcreate(data$BF)
  
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
      ggplot2::geom_line(data=as.data.frame(simci.fft), aes(x=as.numeric(1:n.hz), y=simci.fft), linetype="dotted", linewidth=1)
  }
  # Plot FFT
  p+
    ggplot2::geom_line(data=as.data.frame(data), aes(x=as.numeric(1:length(data)), y=data), color=color, linewidth=1)+
    ggplot2::labs(title=paste0("Fast Fourier Transform"), x=sprintf("Frequency (\u2116 of Cycles)"), y = "Amplitude")+
    ggplot2::coord_cartesian(xlim = c(1,n.hz), ylim = coordy)+
    ggplot2::scale_x_continuous(breaks = seq(0, n.hz, by = (n.hz/(n.hz/2))), expand = c(0,0))+
    #scale_y_continuous(breaks = seq(0, 50, by = 1))+
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(legend.position = "none")
  
}

#' Plot BF Robustness Analyses
#'
#' This function allows to plot a BF robustness analysis (bfttestRobustness())
#'
#' Returns an overview showing the the BFs for different prior parameters (location and scale)
#'
#' @param data An object generated with the bfRobustness()-function.
#' @examples
#' plotrobust(bfRobustness(seqbf))
#' @export
plotrobust <- function(data){
  
  if(inherits(data,"bfRobustness") == FALSE) stop("Please provide a bfRobustness object")
  
  # one Location or more Locations?
  robustnessType <- ifelse(length(unique(data$BFMatrix$prior.loc))==1, "uninformed", "informed")
  
  #get min and max
  min_max_values <- aggregate(data$BFMatrix$bf, by = list(data$BFMatrix$prior.loc), FUN = function(x) c(min = min(x), max = max(x)))
  min_max_values <- data.frame(prior.loc = min_max_values$Group.1,
                               min = min_max_values$x[, "min"],
                               max = min_max_values$x[, "max"])
  
  # Best and Worst BF
  best <- data$BFMatrix[data$BFMatrix$bf==max(data$BFMatrix$bf),]
  worst <- data$BFMatrix[data$BFMatrix$bf==min(data$BFMatrix$bf),]
  user <- data$BFMatrix[round(data$BFMatrix$prior.loc,2)==data$prior[[2]] & round(data$BFMatrix$prior.r,2)==data$prior[[3]],]
  
  # subtitle <- paste("BF =",round(tail(data$BF,n=1),3),"// N =", sum(data$`sample size`))
  title <- "Bayes Factor Robustness Test"
  subtitle <- paste0(
    "Best BF=",round(best$bf,2)," at ",data$prior[[1]],"(",best$prior.loc,",",best$prior.r,")\n",
    "Worst BF=",round(worst$bf,2)," at ",data$prior[[1]],"(",worst$prior.loc,",",worst$prior.r,")\n",
    "User BF=",round(user$bf,2)," at ",data$prior[[1]],"(",user$prior.loc,",",user$prior.r,")"
  )
  
  if(robustnessType=="informed") caption <- paste0(data$prior[[1]]," (Widths: ",min(data$BFMatrix$prior.r)," - ",max(data$BFMatrix$prior.r),")")
  else caption <- paste0(data$prior[[1]]," (Location: ",unique(data$BFMatrix$prior.loc),")")
  
  # plot
  # Set y coordinates
  coordy <- c(min(min_max_values$min),2*max(min_max_values$max, na.rm=T))
  
  #Set minimum coordinates to 1/10 and 10
  if(coordy[1] > 1/10) coordy[1] <- 1/10
  if(coordy[2] < 10) coordy[2] <- 10
  
  # Specify Text annotations
  annotationlist <- .annotations(coordy)
  
  annotation <- annotationlist[[1]]
  annobreaks <- annotationlist[[2]]
  
  #Scale y-Axis
  breaks <- annotationlist[[3]]
  labels <- annotationlist[[4]]
  
  if(robustnessType=="informed"){
    # Draw plot
    p <- ggplot2::ggplot()
    
    # Add horizontal lines
    p <- p + ggplot2::geom_hline(yintercept = 1, color='grey60', linetype = 'solid')+
      ggplot2::geom_hline(yintercept = breaks, color='grey60', linetype='dotted')
    
    # Ribbon
    p <- p + ggplot2::geom_ribbon(data=min_max_values, ggplot2::aes(x=prior.loc, ymin=min, ymax=max), alpha=0.3)
    if(coordy[2] <= 1000 && coordy[1] >= 1/1000) p <- p + ggplot2::annotate("text", x=max(min_max_values$prior.loc)*1.2, y=annobreaks, label=annotation, hjust=1, parse = TRUE)
    
    # Best and Worst BF
    p <- p + 
      #ggplot2::geom_line(data = subset(data$BFMatrix, prior.r==best$prior.r), ggplot2::aes(x=prior.loc, y=bf), color="cornflowerblue", linewidth=1)+
      #ggplot2::geom_line(data = subset(data$BFMatrix, prior.r==worst$prior.r), ggplot2::aes(x=prior.loc, y=bf), color="coral2", linewidth=1)+
      ggplot2::geom_line(data = subset(data$BFMatrix, prior.r==user$prior.r), ggplot2::aes(x=prior.loc, y=bf), color="chartreuse4", linewidth=1)+
      ggplot2::geom_point(data = best, ggplot2::aes(x=prior.loc, y=bf), color="cornflowerblue")+
      ggplot2::geom_point(data = worst, ggplot2::aes(x=prior.loc, y=bf), color="coral2")+
      ggplot2::geom_point(data = user, ggplot2::aes(x=prior.loc, y=bf), color="chartreuse4", shape=18, size=3)+
      ggplot2::geom_label(data = best, ggplot2::aes(x=prior.loc, y=bf, label=round(bf,2)), color="cornflowerblue", vjust=-0.5)+
      ggplot2::geom_label(data = worst, ggplot2::aes(x=prior.loc, y=bf, label=round(bf,2)), color="coral2", vjust=1.5)+
      ggplot2::geom_label(data = user, ggplot2::aes(x=prior.loc, y=bf, label=round(bf,2)), color="chartreuse4", vjust=1.5)
    
    p <- p + ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
    
    p + ggplot2::labs(x="Prior Location", y = "Evidence (BF)")+
      ggplot2::scale_y_log10(breaks = breaks, labels = labels)+
      ggplot2::coord_cartesian(ylim = coordy)+
      ggplot2::theme_bw(base_size = 14)+
      ggplot2::theme(legend.position = "none", panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
    
  } else {
    # Draw plot
    p <- ggplot2::ggplot()
    
    # Add horizontal lines
    p <- p + ggplot2::geom_hline(yintercept = 1, color='grey60', linetype = 'solid')+
      ggplot2::geom_hline(yintercept = breaks, color='grey60', linetype='dotted')
    
    if(coordy[2] <= 1000 && coordy[1] >= 1/1000) p <- p + ggplot2::annotate("text", x=max(data$BFMatrix$prior.r)*1.2, y=annobreaks, label=annotation, hjust=1, parse = TRUE)
    
    # Best and Worst BF
    p <- p + 
      #ggplot2::geom_line(data = subset(data$BFMatrix, prior.r==best$prior.r), ggplot2::aes(x=prior.loc, y=bf), color="cornflowerblue", linewidth=1)+
      #ggplot2::geom_line(data = subset(data$BFMatrix, prior.r==worst$prior.r), ggplot2::aes(x=prior.loc, y=bf), color="coral2", linewidth=1)+
      ggplot2::geom_line(data = data$BFMatrix, ggplot2::aes(x=prior.r, y=bf), color="black", linewidth=1)+
      ggplot2::geom_point(data = best, ggplot2::aes(x=prior.r, y=bf), color="cornflowerblue")+
      ggplot2::geom_point(data = worst, ggplot2::aes(x=prior.r, y=bf), color="coral2")+
      ggplot2::geom_point(data = user, ggplot2::aes(x=prior.r, y=bf), color="chartreuse4", shape=18, size=3)+
      ggplot2::geom_label(data = best, ggplot2::aes(x=prior.r, y=bf, label=round(bf,2)), color="cornflowerblue", vjust=-0.5)+
      ggplot2::geom_label(data = worst, ggplot2::aes(x=prior.r, y=bf, label=round(bf,2)), color="coral2", vjust=1.5)+
      ggplot2::geom_label(data = user, ggplot2::aes(x=prior.r, y=bf, label=round(bf,2)), color="chartreuse4", vjust=1.5)
    
    p <- p + ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
    
    p + ggplot2::labs(x="Prior Width", y = "Evidence (BF)")+
      ggplot2::scale_y_log10(breaks = breaks, labels = labels)+
      ggplot2::coord_cartesian(ylim = coordy)+
      ggplot2::theme_bw(base_size = 14)+
      ggplot2::theme(legend.position = "none", panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
    
  }
}


#' Tile Plot BF Robustness Analyses
#'
#' This function allows to plot a BF robustness analysis (bfttestRobustness())
#'
#' Returns a heatmap showing the the BFs for different prior parameters (location and scale)
#'
#' @param data An object generated with the bfRobustness()-function.
#' @param limit A BF limit which is marked in the plot.
#' @examples
#' plotrobustTile(bfRobustness(seqbf))
#' @export

plotrobustTile <- function(data, limit=10){
  
  if(inherits(data,"bfRobustness") == FALSE) stop("Please provide a bfRobustness object")
  
  maxvalue <- max(c(10,max(data$BFMatrix$bf)))
  breaks <- c(0,1,3,6,10,maxvalue)
  min <- subset(data$BFMatrix, bf==min(bf))
  max <- subset(data$BFMatrix, bf==max(bf))
  
  data$BFMatrix$col <- ifelse(data$BFMatrix$bf >= limit, TRUE, FALSE)
  
  ggplot2::ggplot(data=data$BFMatrix, ggplot2::aes(x=prior.loc, y=prior.r, fill=bf))+
    ggplot2::geom_raster()+
    ggplot2::geom_tile(data=subset(data$BFMatrix, col==TRUE), color="black")+
    ggplot2::geom_text(data = min, ggplot2::aes(x=prior.loc, y=prior.r, label = round(bf,2)), color = "black", size = 3) +
    ggplot2::geom_text(data = max, ggplot2::aes(x=prior.loc, y=prior.r, label = round(bf,2)), color = "black", size = 3) +
    ggplot2::scale_fill_gradientn(colors = c("cornflowerblue",
                                             "white",
                                             "green3",
                                             "yellow",
                                             "coral2",
                                             "darkred"),
                                  breaks = breaks,
                                  labels = scales::label_number(accuracy = 1),
                                  limits = c(0,maxvalue),
                                  values = scales::rescale(breaks))+
    ggplot2::labs(x="Prior Location", y="Prior Width", fill="BF10")+
    ggplot2::theme_minimal()
  
}

#' @export
plot.bfRobustness <- function(data, ...) {
  suppressWarnings(plotrobust(data, ...))
}


# Helper functions
.annotations <- function(coordy){
  
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
  
  list(annotation,
       annobreaks,
       breaks,
       labels)
}

