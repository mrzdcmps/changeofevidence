% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coetest.R
\name{energybf}
\alias{energybf}
\title{BF Energy Analysis}
\usage{
energybf(data, sims.df = sims)
}
\arguments{
\item{data}{A seqbf object or vector containing Bayes Factors.}

\item{sims.df}{A dataframe containing simulations, including columns "simid" and "bf".}
}
\value{
A list containing the Energy of the experimental data, the mean and SD of the energy values of all simulations, and the proportion of simulations with a higher energy than the experimental data.
}
\description{
This function calculates the Energy of the evidence over time and evaluates its likelihood.
}
\details{
The Energy contains all areas above BF=1 minus all areas below BF=1. A positive energy therefore indicates an overall orientation of the evidence towards H1 over data collection.
Energy is calculated using the trapezoidal integration "trapz" of the pracma-library.
Mean and SD scores of the simulations' energies are provided.
}
\examples{
r.energybf <- energybf(seqbf, sims)
r.energybf <- energybf(tbl$bf, sims.df = newsims)
}
