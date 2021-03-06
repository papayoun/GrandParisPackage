# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title Proposition model
#' @name ProposalSINEModel
#' @description POD with the unidimensional SDE sinus model dX_t = sin(X_t - theta)dt + dW_t
#' and the observation model Y_t = X_t + N(0, sigma2), the proposal model is a random walk 
#' @export ProposalSINEModel
NULL

#' @title Class of partially diffusion process ruled by a SINE model
#' @name SINE_POD
#' @description POD with the unidimensional SDE sinus model dX_t = sin(X_t - theta)dt + dW_t
#' and the observation model Y_t = X_t + N(0, sigma2) 
#' @export SINE_POD
NULL

#' @title Main program for tangent online estimation
#' @name fastTangOR
#' @export
fastTangOR <- function(observations, observationTimes, thetaModel, sigma2, updateOrders, gradientSteps, randomWalkParam = 2, particleSize = 100L, densitySampleSize = 30L, logDensitySampleSize = 30L, backwardSampleSize = 2L, backwardSamplingMaxTry = 100000000L, skeletonSimulationMaxTry = 10000000L, estimateTheta = TRUE, estimateSigma2 = TRUE, all = FALSE, IS = FALSE) {
    .Call('_GrandParisPackage_fastTangOR', PACKAGE = 'GrandParisPackage', observations, observationTimes, thetaModel, sigma2, updateOrders, gradientSteps, randomWalkParam, particleSize, densitySampleSize, logDensitySampleSize, backwardSampleSize, backwardSamplingMaxTry, skeletonSimulationMaxTry, estimateTheta, estimateSigma2, all, IS)
}

#' @title EM algortihm
#' @name GEM
#' @export
GEM <- function(observations, observationTimes, thetaStart, sigma2Start, nIterations = 20L, nModels = 5L, randomWalkParam = 2, particleSize = 100L, densitySampleSize = 30L, logDensitySampleSize = 30L, backwardSampleSize = 2L, backwardSamplingMaxTry = 100000000L, skeletonSimulationMaxTry = 10000000L, estimateTheta = TRUE, estimateSigma2 = TRUE) {
    .Call('_GrandParisPackage_GEM', PACKAGE = 'GrandParisPackage', observations, observationTimes, thetaStart, sigma2Start, nIterations, nModels, randomWalkParam, particleSize, densitySampleSize, logDensitySampleSize, backwardSampleSize, backwardSamplingMaxTry, skeletonSimulationMaxTry, estimateTheta, estimateSigma2)
}

#' @title EM algortihm
#' @name GEM_IS
#' @export
GEM_IS <- function(observations, observationTimes, thetaStart, sigma2Start, nIterations = 20L, nModels = 5L, randomWalkParam = 2, particleSize = 100L, densitySampleSize = 30L, logDensitySampleSize = 30L, backwardSampleSize = 2L, backwardSamplingMaxTry = 100000000L, skeletonSimulationMaxTry = 10000000L, estimateTheta = TRUE, estimateSigma2 = TRUE) {
    .Call('_GrandParisPackage_GEM_IS', PACKAGE = 'GrandParisPackage', observations, observationTimes, thetaStart, sigma2Start, nIterations, nModels, randomWalkParam, particleSize, densitySampleSize, logDensitySampleSize, backwardSampleSize, backwardSamplingMaxTry, skeletonSimulationMaxTry, estimateTheta, estimateSigma2)
}

#' @title Tracking algortihm
#' @name E_track
#' @export
E_track <- function(observations, observationTimes, ind_tracked, thetaStart, sigma2Start, nIterations = 20L, randomWalkParam = 2, particleSize = 100L, densitySampleSize = 30L, logDensitySampleSize = 30L, backwardSampleSize = 2L, backwardSamplingMaxTry = 100000000L, skeletonSimulationMaxTry = 10000000L, estimateTheta = TRUE, estimateSigma2 = TRUE) {
    .Call('_GrandParisPackage_E_track', PACKAGE = 'GrandParisPackage', observations, observationTimes, ind_tracked, thetaStart, sigma2Start, nIterations, randomWalkParam, particleSize, densitySampleSize, logDensitySampleSize, backwardSampleSize, backwardSamplingMaxTry, skeletonSimulationMaxTry, estimateTheta, estimateSigma2)
}

#' @title Tracking with IS
#' @name track_IS
#' @export
E_track_IS <- function(observations, observationTimes, ind_tracked, thetaStart, sigma2Start, nIterations = 20L, randomWalkParam = 2, particleSize = 100L, densitySampleSize = 30L, logDensitySampleSize = 30L, backwardSampleSize = 2L, backwardSamplingMaxTry = 100000000L, skeletonSimulationMaxTry = 10000000L, estimateTheta = TRUE, estimateSigma2 = TRUE) {
    .Call('_GrandParisPackage_E_track_IS', PACKAGE = 'GrandParisPackage', observations, observationTimes, ind_tracked, thetaStart, sigma2Start, nIterations, randomWalkParam, particleSize, densitySampleSize, logDensitySampleSize, backwardSampleSize, backwardSamplingMaxTry, skeletonSimulationMaxTry, estimateTheta, estimateSigma2)
}

