# Suppress R CMD check notes for ggplot2 aesthetic variables
utils::globalVariables(c(
  ".data",
  "n", "p.dn", "p.up", "rw", "group",
  "index", "simid", "x", "y", "element",
  "bf", "prior.loc", "prior.r"
))

#' Plot Random Walk
#'
#' This function plots random walks with appropriate confidence intervals.
#'
#' The Random Walk can be plotted by itself or in comparison to simulated data sets.
#' The function supports both classic binary random walks (+1/-1 steps) and deviation-based walks (e.g., from summed bits).
#' Depending on the amount of simulations, drawing might take a while. It might be wise to choose a smaller simulation set for this purpose.
#'
#' @param ... One or more numeric vectors containing the random walk(s) to be drawn, or a single list of vectors.
#' @param labels Optional character vector of names for multiple random walks. If NULL, auto-generates names.
#' @param sims.df A dataframe containing simulations, including columns "simid" and "index". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column in the simulation dataframe to compare to.
#' @param color A color or vector of colors for the random walk line(s). A single color is used for one line;
#'   a vector of colors applies custom colors to multiple lines (overrides default palette).
#' @param coordy A vector containing the minimum and maximum value of the y-coordinates to be drawn.
#' @param p Probability of success for binary random walks (0 < p < 1). Default is 0.5.
#' @param n_bits Number of bits being summed per trial (used with deviation-based walks).
#' 
#' @examples
#' # Example 1: Classic binary random walk (+1/-1 steps)
#' responses <- sample(c(TRUE, FALSE), 100, replace = TRUE)
#' rw_binary <- cumsum(ifelse(responses, 1, -1))
#' plotrw(rw_binary)
#'
#' # Example 2: Binary walk with biased probability
#' responses_biased <- sample(c(TRUE, FALSE), 100, replace = TRUE, prob = c(0.6, 0.4))
#' rw_biased <- cumsum(ifelse(responses_biased, 1, -1))
#' plotrw(rw_biased, p = 0.6)
#'
#' # Example 3: Deviation-based walk from summed bits (10 bits per trial)
#' bit_sums <- rbinom(100, 10, 0.5)
#' rw_deviation <- cumsum(bit_sums - 5)
#' plotrw(rw_deviation, n_bits = 10)
#'
#' # Example 4: Multiple random walks comparison
#' condition_a <- rbinom(80, 10, 0.5)
#' condition_b <- rbinom(80, 10, 0.55)
#' rw_list <- list(control = cumsum(condition_a - 5), experimental = cumsum(condition_b - 5))
#' plotrw(rw_list, n_bits = 10)
#'
#' # Example 5: With custom y-axis limits
#' plotrw(rw_deviation, n_bits = 10, coordy = c(-20, 20))
#'
#' \donttest{
#' # Example 6: With simulation data for comparison (slow - 1000 simulations)
#' sims_data <- do.call(rbind, lapply(1:100, function(i) {
#'   sim_rw <- cumsum(rbinom(100, 10, 0.5) - 5)
#'   data.frame(simid = i, index = 1:100, rw = sim_rw)
#' }))
#' plotrw(rw_deviation, sims.df = sims_data, sims.df.col = "rw", n_bits = 10)
#' }
#'
#' # Example 7: Different bit sizes
#' bit_sums_20 <- rbinom(100, 20, 0.5)
#' rw_20bits <- cumsum(bit_sums_20 - 10)
#' plotrw(rw_20bits, n_bits = 20)
#'
#' # Example 8: Biased bits (e.g., p=0.2)
#' bit_sums_biased <- rbinom(100, 10, 0.2)
#' rw_biased_bits <- cumsum(bit_sums_biased - 2)
#' plotrw(rw_biased_bits, n_bits = 10, p = 0.2)
#' @export

plotrw <- function(..., labels = NULL, sims.df = NULL, sims.df.col = "rw",
                   color = "black", coordy = NULL, p = 0.5, n_bits = NULL) {
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 200, maxColorValue = 255)

  # Collect inputs from ...
  args <- list(...)
  if (length(args) == 0) stop("No data provided")

  if (length(args) == 1) {
    data <- args[[1]]
    # Single list of vectors passed directly - apply labels if given
    if (is.list(data) && !is.null(labels)) {
      names(data) <- labels[seq_along(data)]
    }
  } else {
    # Multiple vectors passed via ... - combine into a named list
    data <- args
    if (!is.null(labels)) {
      names(data) <- labels[seq_along(data)]
    } else {
      names(data) <- paste0("Data ", seq_along(data))
    }
  }

  # Adjust the random walk if p != 0.5 (only for binary walks)
  adjust_rw <- function(rw, p) {
    if (p != 0.5 && is.null(n_bits)) {
      n <- length(rw) - 1
      return(rw + (2 * p - 1) * (0:n))
    }
    return(rw)
  }

  # Normalise to list internally; track whether multiple series are shown
  if (is.list(data)) {
    show.legend <- "bottom"
    nmax <- max(lengths(data))
    absolutemax <- max(c(max(unlist(data)), abs(min(unlist(data)))))
    for (i in seq_along(data)) {
      data[[i]] <- c(0, data[[i]])
      data[[i]] <- adjust_rw(data[[i]], p)
    }
  } else {
    show.legend <- "none"
    nmax <- length(data)
    absolutemax <- max(c(max(data), abs(min(data))))
    data <- c(0, data)
    data <- adjust_rw(data, p)
  }

  # Set coordy if not provided
  if (is.null(coordy)) {
    margin_factor <- 1.2
    coordy <- c(-absolutemax * margin_factor, absolutemax * margin_factor)
  }

  # Confidence bounds (p-parabola)
  z <- 1.96
  p.s <- data.frame(n = 0:nmax)
  if (is.null(n_bits)) {
    p.s$p.up <- (2 * p - 1) * p.s$n + z * sqrt(4 * p.s$n * p * (1 - p))
    p.s$p.dn <- (2 * p - 1) * p.s$n - z * sqrt(4 * p.s$n * p * (1 - p))
  } else {
    variance_per_step <- n_bits * p * (1 - p)
    p.s$p.up <-  z * sqrt(variance_per_step * p.s$n)
    p.s$p.dn <- -z * sqrt(variance_per_step * p.s$n)
  }

  # Base plot
  plot <- ggplot() +
    geom_ribbon(data = p.s, aes(x = n, ymin = p.dn, ymax = p.up), fill = greycol, alpha = 0.3) +
    geom_line(data = p.s, aes(x = n, y = p.up), color = "grey50", linetype = "dashed") +
    geom_line(data = p.s, aes(x = n, y = p.dn), color = "grey50", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    coord_cartesian(ylim = coordy) +
    labs(x = "Trial", y = "Cumulative Deviation", color = "") +
    theme_minimal(base_size = 14) +
    theme(legend.position = show.legend)

  # Add simulation data if provided
  if (!is.null(sims.df)) {
    plot <- plot +
      geom_line(data = sims.df,
                aes(x = index, y = .data[[sims.df.col]], group = simid),
                color = "grey70", alpha = 0.2, linewidth = 0.3)
  }

  # Add actual data
  if (is.list(data)) {
    plot_data <- data.frame()
    for (i in seq_along(data)) {
      temp_df <- data.frame(
        n     = 0:(length(data[[i]]) - 1),
        rw    = data[[i]],
        group = names(data)[i]
      )
      plot_data <- rbind(plot_data, temp_df)
    }
    # Preserve factor order for legend
    plot_data$group <- factor(plot_data$group, levels = names(data))
    plot <- plot +
      geom_line(data = plot_data, aes(x = n, y = rw, color = group), linewidth = 1)
    if (length(color) > 1) {
      plot <- plot + scale_color_manual(values = color)
    } else {
      plot <- plot + scale_color_brewer(palette = "Set1")
    }
  } else {
    plot_data <- data.frame(n = 0:(length(data) - 1), rw = data)
    plot <- plot +
      geom_line(data = plot_data, aes(x = n, y = rw), color = color[1], linewidth = 1)
  }

  return(plot)
}




#' Plot Sequential Bayesian Analysis
#'
#' This function allows to plot a sequential Bayesian analysis.
#'
#' The BF analysis can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
#' BFs above 1 indicate evidence towards H1, BFs below 1 indicate evidence towards H0.
#' Depending on the amount of simulations, drawing might take a while. It might be wise to chose a smaller simulation set for this purpose.
#'
#' @param ... One or more seqbf objects, vectors containing sequential Bayes Factors, or a single list containing multiple vectors.
#'   When multiple seqbf objects or vectors are provided, they will be plotted together for comparison.
#' @param labels Optional character vector of names for multiple datasets. If NULL, auto-generates names.
#' @param sims.df A dataframe containing simulations, including column "simid". Set to NULL if you don't want to display simulations.
#' @param sims.df.col The name of the column of the simulation dataframe to compare to. Default is "bf".
#' @param color A color or vector of colors for the BF line(s). A single color is used for one line;
#'   a vector of colors applies custom colors to multiple lines (overrides the default palette).
#' @param coordy A vector containing the minimum and maximum value of the y-coordinates to be drawn. If NULL, automatically determined.
#' @param label.x A character that overrides the label for the x-axis. Default is "N" (automatically set to "Trials" for binomial data).
#' @param show_annotations Logical. If TRUE (default), displays evidence strength annotations (e.g., "Moderate H1", "Strong H0"). Set to FALSE to hide annotations and use full plot width.
#' 
#' @return A ggplot2 object.
#' 
#' @examples
#' \dontrun{
#' # Single seqbf object
#' plot(seqbf)
#' plotbf(seqbf)
#'
#' # Single BF vector
#' p.bf <- plotbf(tbl$bf)
#' p.bf
#'
#' # Multiple seqbf objects with custom labels
#' plotbf(bf1, bf2, labels = c("Experiment 1", "Experiment 2"))
#'
#' # Multiple seqbf objects with auto-generated labels
#' plotbf(bf1, bf2, bf3)
#'
#' # List of BF vectors (backward compatible)
#' plotbf(list(test1 = bf1$BF, test2 = bf2$BF))
#'
#' # With simulated data
#' plotbf(tbl$bf, sims.df = sims)
#'
#' # With subset of simulations
#' sims1000 <- subset(sims, simid <= 1000)
#' plotbf(tbl$bf, sims.df = sims1000)
#'
#' # Custom y-axis limits
#' plotbf(bf1, bf2, coordy = c(1/30, 30))
#' }
#' @export

# Plot Sequential BF
plotbf <- function(..., labels = NULL, sims.df = NULL, sims.df.col = "bf", color = "black", coordy = NULL, label.x = "N", show_annotations = TRUE){

  # Capture all arguments passed via ...
  args <- list(...)
  
  # Handle different input patterns
  if(length(args) == 0) {
    stop("No data provided")
  }
  
  # Single argument - could be seqbf, vector, or list
  if(length(args) == 1) {
    data <- args[[1]]
    # If it's a list but not a seqbf object, treat as list of BF vectors (backward compat)
    if(is.list(data) && !inherits(data, "seqbf")) {
      # This is the old list(bf1$BF, bf2$BF) syntax - keep as is
      data <- data
    }
    # else: single seqbf object or vector - will be handled below
  } else {
    # Multiple arguments - convert seqbf objects to BF vectors
    data <- list()
    for(i in seq_along(args)) {
      if(inherits(args[[i]], "seqbf")) {
        data[[i]] <- args[[i]]$BF
      } else {
        data[[i]] <- args[[i]]
      }
    }
    
    # Apply labels
    if(!is.null(labels)) {
      if(length(labels) != length(data)) {
        warning("Length of labels does not match number of data objects")
      }
      names(data) <- labels
    } else {
      # Auto-generate names
      names(data) <- paste0("Data ", 1:length(data))
    }
  }
  
  # Continue with existing logic for seqbf object
  if(inherits(data,"seqbf") == TRUE){
    
    # Determine test type and name
    if(data$`test type` == "independent") {
      if(!is.null(data$parametric) && !data$parametric) {
        testtype <- "Mann-Whitney U"
      } else {
        testtype <- "Independent Samples t-test"
      }
    } else if(data$`test type` == "paired") {
      if(!is.null(data$parametric) && !data$parametric) {
        testtype <- "Wilcoxon Signed-Rank (Paired)"
      } else {
        testtype <- "Paired Samples t-test"
      }
    } else if(data$`test type` == "one-sample") {
      if(!is.null(data$parametric) && !data$parametric) {
        testtype <- "Wilcoxon Signed-Rank (One-Sample)"
      } else {
        testtype <- "One-Sample t-test"
      }
    } else if(data$`test type` == "binomial") {
      testtype <- "Binomial"
      label.x <- "Trials"
    } else if(data$`test type` == "correlation") {
      testtype <- "Correlation"
    } else {
      testtype <- "Unknown"
    }
    
    # Determine tails
    tails <- ifelse(data$alternative == "two.sided", "two-tailed", "one-tailed")
    
    # Create subtitle with final BF and sample size
    subtitle <- paste0("BF = ", round(tail(data$BF, n = 1), 3), 
                       " (N = ", sum(data$`sample size`), ")")
    
    # Create caption with test info
    caption <- paste0(
      testtype, "; ",
      tails, "; ",
      data$prior$distribution, "(", 
      round(data$prior[[2]], 3), ", ", 
      round(data$prior[[3]], 3), ")")
    
    # Add delta to caption if available
    final_delta <- tail(na.omit(data$delta), n = 1)
    if(length(final_delta) > 0 && !is.na(final_delta)) {
      caption <- paste0(caption, "; \u03b4 = ", round(final_delta, 3))
    }
    
    data <- data$BF
    showinfo <- TRUE
  }
  
  greycol <- rgb(red = 190, green = 190, blue = 190, alpha = 150, maxColorValue = 255)

  # Set y coordinates
  if(is.null(coordy)){
    if(is.list(data)){
      coordy <- c(min(unlist(data), na.rm = T),
                  2 * max(unlist(data), na.rm = T))
    }
    else {
      coordy <- c(min(data, na.rm = T), 2 * max(data, na.rm = T))
    }
  }
  
  # Set minimum coordinates to 1/10 and 10
  if(coordy[1] > 1/10) coordy[1] <- 1/10
  if(coordy[2] < 10) coordy[2] <- 10
  
  # Show legend?
  if(is.list(data)) show.legend <- "bottom"
  else show.legend <- "none"
  
  # Specify Text annotations
  annotationlist <- .annotations(coordy)
  
  annotation <- annotationlist[[1]]
  annobreaks <- annotationlist[[2]]
  
  # Scale y-Axis
  breaks <- annotationlist[[3]]
  labels <- annotationlist[[4]]
  breaks_filtered <- annotationlist[[5]]
  
  # Print message if simulations are drawn
  if(!is.null(sims.df)) print("Depending on the amount of simulations to be drawn, this might take a while!")
  
  # Initialize plot
  p <- ggplot2::ggplot()
  
  # Add simulations
  if(!is.null(sims.df)) p <- p + ggplot2::geom_line(data = sims.df, aes(x = index, y = .data[[sims.df.col]], group = simid), color = greycol)
  
  # Add horizontal lines
  p <- p + ggplot2::geom_hline(yintercept = 1, color = 'grey60', linetype = 'solid') +
    ggplot2::geom_hline(yintercept = breaks_filtered, color = 'grey60', linetype = 'dotted')
  
  # Unlist data if multiple bfs are provided
  if(is.list(data)){
    df <- NULL
    for(i in 1:length(data)){
      ydat <- data[[i]]
      # Ensure ydat is a simple numeric vector
      if (is.data.frame(ydat)) {
        ydat <- as.vector(as.matrix(ydat))
      } else if (is.list(ydat) && !is.data.frame(ydat)) {
        ydat <- unlist(ydat)
      }
      ydat <- as.numeric(ydat)
      if(!is.null(names(data))){
        ndf <- data.frame(element = as.factor(rep(names(data)[i], length(ydat))), x = 1:length(ydat), y = ydat)
      } else{
        ndf <- data.frame(element = rep(paste0("data ", i), length(ydat)), x = 1:length(ydat), y = ydat)
      }
      df <- rbind(df, ndf)
      df <- df[!is.na(df$y), ]
    }
    p <- p + ggplot2::geom_line(data = df, aes(x = x, y = y, color = element), linewidth = 1)
    if (length(color) > 1) {
      p <- p + ggplot2::scale_color_manual("", values = color)
    } else {
      p <- p + ggplot2::scale_color_brewer("", type = "qualitative", palette = "Set1")
    }

    # Add annotations outside plot area
    if(show_annotations && coordy[2] <= 1000 && coordy[1] >= 1/1000) {
      p <- p + ggplot2::annotate("text",
                                 x = Inf,
                                 y = annobreaks,
                                 label = annotation,
                                 hjust = -0.05,
                                 size = 3.5,
                                 parse = TRUE,
                                 color = "grey40")
    }
    
    # Add single data  
  } else {
    df <- data.frame(y = data, x = as.numeric(1:length(data)))
    df <- df[!is.na(df$y), ]
    
    p <- p + ggplot2::geom_line(data = df, aes(x = x, y = y), color = color, linewidth = 1)

    # Add annotations outside plot area
    if(show_annotations && coordy[2] <= 1000 && coordy[1] >= 1/1000) {
      p <- p + ggplot2::annotate("text",
                                 x = Inf,
                                 y = annobreaks,
                                 label = annotation,
                                 hjust = -0.05,
                                 size = 3.5,
                                 parse = TRUE,
                                 color = "grey40")
    }
  }
  
  # Add subtitle and caption if seqbf object
  if(exists("showinfo")) p <- p + ggplot2::labs(subtitle = subtitle, caption = caption)
  
  # Finalize plot
  # Adjust margins based on whether annotations are shown
  right_margin <- if(show_annotations) 70 else 5

  p + ggplot2::labs(x = label.x, y = "Evidence (BF)") +
    ggplot2::scale_y_log10(breaks = breaks, labels = labels) +
    ggplot2::scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    ggplot2::coord_cartesian(ylim = coordy, clip = if(show_annotations) "off" else "on") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = show.legend,
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.margin = margin(5, right_margin, 5, 5, "pt"))
  
}

#' @export
plot.seqbf <- function(x, ...) {
  suppressWarnings(plotbf(x, ...))
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
#' \dontrun{
#' p.fftbf <- plotfft(tblFFT$density.bf, color = "blue")
#' p.fftbf
#'
#' p.fftrw <- plotfft(tblFFT$density.rw, sims.df = sims, sims.df.col = "density.rw")
#' p.fftrw
#' }
#' @export

# Plot FFT
# Data for 95-CI ribbon FFT
plotfft <- function(data, sims.df = NULL, sims.df.col = "density.bf", n.hz = 50, color = "black", coordy = c(0,secondhighestval)){
  if(inherits(data,"seqbf") == TRUE) data <- changeofevidence::fftcreate(data$BF)
  
  secondhighestval <- sort(data,partial=length(data)-1)[length(data)-1]
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
#' \dontrun{
#' plotrobust(bfRobustness(seqbf))
#' }
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
#' \dontrun{
#' plotrobustTile(bfRobustness(seqbf))
#' }
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
plot.bfRobustness <- function(x, ...) {
  suppressWarnings(plotrobust(x, ...))
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
  
  # Filter annotations to only include those within coordy range
  within_range <- annobreaks >= coordy[1] & annobreaks <= coordy[2]
  annotation <- annotation[within_range]
  annobreaks <- annobreaks[within_range]
  
  #Scale y-Axis
  if(coordy[2] <= 1000 && coordy[1] >= 1/1000){
    breaks = c(1000,300,100,30,10,3,1,1/3,1/10,1/30,1/100,1/300,1/1000)
    labels = c("1,000","300","100","30","10","3","1","1/3","1/10","1/30","1/100","1/300","1,1,000")
  } else if(coordy[2] <= 1000000 && coordy[1] >= 1/1000000){
    breaks = c(1000000,100000,10000,1000,100,10,1,1/10,1/100,1/1000,1/10000,1/100000,1/1000000)
    labels = c("1,000,000","100,000","10,000","1,000","100","10","1","1/10","1/100","1/1,000","1/10,000","1/100,000","1/1,000,000")
  } else {
    breaks = c(1e+12,1e+10,1e+8,1e+6,1e+4,1e+2,1,1/1e+2,1/1e+4,1/1e+6,1/1e+8,1/1e+10,1/1e+12)
    labels = c("1e+12","1e+10","1e+8","1e+6","1e+4","1e+2","1","1/1e+2","1/1e+4","1/1e+6","1/1e+8","1/1e+10","1/1e+12")
  }
  
  # Filter threshold breaks for horizontal lines
  breaks_filtered <- breaks[breaks >= coordy[1] & breaks <= coordy[2]]
  
  list(annotation,
       annobreaks,
       breaks,
       labels,
       breaks_filtered)  # Add filtered breaks for horizontal lines
}

