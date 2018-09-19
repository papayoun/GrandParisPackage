// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "DebugMethods.h"
#include "GenericFunctions.h"
#include "SINEModel.h"
#include "SINE_POD.h"
#include "ProposalSINEModel.h"
#include "SDEParticleSmoother.h"
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R

//' @title Main program for tangent online estimation
//' @name fastTangOR
//' @export
// [[Rcpp::export]]
Rcpp::List fastTangOR(const Rcpp::NumericVector& observations,
                const Rcpp::NumericVector& observationTimes,
                const double& thetaModel, const double& sigma2,
                const Rcpp::LogicalVector& updateOrders,
                Rcpp::NumericMatrix& gradientSteps,
                const double& randomWalkParam = 2,
                const unsigned int& particleSize = 100,
                const unsigned int& densitySampleSize = 30,
                const unsigned int& logDensitySampleSize = 30,
                const unsigned int& backwardSampleSize = 2,
                const unsigned int& backwardSamplingMaxTry = 100000000,
                const unsigned int& skeletonSimulationMaxTry = 10000000,
                const bool& estimateTheta = true, const bool& estimateSigma2 = true,
                const bool& all = false){
  SINE_POD trueModel(thetaModel, sigma2);
  ProposalSINEModel propModel(randomWalkParam, trueModel, estimateTheta, estimateSigma2);
  SDEParticleSmoother mySmoother(observations, observationTimes,  propModel,
                                 particleSize, densitySampleSize , logDensitySampleSize,
                                 backwardSampleSize, backwardSamplingMaxTry, skeletonSimulationMaxTry);
  Rcpp::NumericMatrix output  = mySmoother.tangentFilterEstimation(updateOrders, gradientSteps);
  if(all){
    return Rcpp::List::create(Rcpp::Named("Estimates") = output,
                              Rcpp::Named("Particles") = mySmoother.getParticles(),
                              Rcpp::Named("Weights") = mySmoother.getWeights());
  }
  else
    return Rcpp::List::create(Rcpp::Named("Estimates") = output);
}
