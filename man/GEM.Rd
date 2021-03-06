% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GEM}
\alias{GEM}
\title{EM algortihm}
\usage{
GEM(observations, observationTimes, thetaStart, sigma2Start,
  nIterations = 20L, nModels = 5L, randomWalkParam = 2,
  particleSize = 100L, densitySampleSize = 30L,
  logDensitySampleSize = 30L, backwardSampleSize = 2L,
  backwardSamplingMaxTry = 100000000L,
  skeletonSimulationMaxTry = 10000000L, estimateTheta = TRUE,
  estimateSigma2 = TRUE)
}
\description{
EM algortihm
}
