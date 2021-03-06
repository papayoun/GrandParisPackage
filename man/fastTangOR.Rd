% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fastTangOR}
\alias{fastTangOR}
\title{Main program for tangent online estimation}
\usage{
fastTangOR(observations, observationTimes, thetaModel, sigma2,
  updateOrders, gradientSteps, randomWalkParam = 2,
  particleSize = 100L, densitySampleSize = 30L,
  logDensitySampleSize = 30L, backwardSampleSize = 2L,
  backwardSamplingMaxTry = 100000000L,
  skeletonSimulationMaxTry = 10000000L, estimateTheta = TRUE,
  estimateSigma2 = TRUE, all = FALSE)
}
\description{
Main program for tangent online estimation
}
