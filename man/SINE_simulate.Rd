% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SINEModel_functions.R
\name{SINE_simulate}
\alias{SINE_simulate}
\title{Simulate a SINE process}
\usage{
SINE_simulate(theta = 3.14, sigma2 = 1, x0 = 0, times = seq(0, 1,
  length.out = 10))
}
\arguments{
\item{theta}{Theta parameter of the model, default to pi}

\item{sigma2}{observation process variance}

\item{x0}{starting value}

\item{times}{vector of simulation times}
}
\value{
a vector same length as times
}
\description{
function that simulates a SINE model trajectory using exact algorithm
}
