% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotrw}
\alias{plotrw}
\title{Plot Random Walk}
\usage{
plotrw(
  data,
  sims.df = NULL,
  sims.df.col = "rw",
  color = "black",
  coordy = c(-absolutemax, absolutemax),
  mu = NULL
)
}
\arguments{
\item{data}{A vector containing the random walk to be drawn or a list containing multiple vectors.}

\item{sims.df}{A dataframe containing simulations, including column "simid" and "index". Set to NULL if you don't want to display simulations.}

\item{sims.df.col}{The name of the column of the simulation dataframe to compare to.}

\item{color}{A color in which the Random Walk will be drawn.}

\item{coordy}{A vector containing the minimum and maximum value of the y-coordinates to be drawn.}

\item{mu}{If Random Walk is of summed up bits, indicate the expected sum per step.}
}
\description{
This function allows to plot random walks.
}
\details{
The Random Walk can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
A random walk typically goes +1 for every correct choice and -1 for every incorrect one.
Depending on the amount of simulations, drawing might take a while. It might be wise to chose a smaller simulation set for this purpose.
}
\examples{
p.rw <- plotrw(tbl$rw)
p.rw

plotrw(tbl$rw, sims.df = sims, sims.df.col = "rw", coordy = c(-50,50))

plotrw(list(exp = exp$rw, con = con$rw), mu = 50)

sims1000 <- subset(sims, simid <= 1000)
plotrw(tbl, sims.df = sims1000)
}
