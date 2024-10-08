% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotbf}
\alias{plotbf}
\title{Plot Sequential Bayesian Analysis}
\usage{
plotbf(
  data,
  sims.df = NULL,
  sims.df.col = "bf",
  color = "black",
  coordy = NULL,
  label.x = "N"
)
}
\arguments{
\item{data}{A seqbf object or a vector containing sequential Bayes Factors or a list containing multiple vectors.}

\item{sims.df}{A dataframe containing simulations, including column "simid". Set to NULL if you don't want to display simulations.}

\item{sims.df.col}{The name of the column of the simulation dataframe to compare to.}

\item{color}{A color in which the Seq BF-function will be drawn.}

\item{coordy}{A vector containing the minimum and maximum value of the y-coordinates to be drawn.}

\item{label.x}{A character that overrides the label for the x-axis ("N" for sum scores, "Trials" for binomial data).}
}
\description{
This function allows to plot a sequential Bayesian analysis.
}
\details{
The BF analysis can be plotted by itself or in comparison to simulated data sets. GGplot2 is used to draw an image.
BFs above 1 indicate evidence towards H1, BFs below 1 indicate evidence towards H0.
Depending on the amount of simulations, drawing might take a while. It might be wise to chose a smaller simulation set for this purpose.
}
\examples{

plot(seqbf)

p.bf <- plotbf(tbl$bf)
p.bf

plotbf(list(test1=bf1,test2=bf2))

plotbf(tbl$bf, sims.df = sims)

sims1000 <- subset(sims, simid <= 1000)
plotbf(tbl$bf, sims.df = sims1000)
}
