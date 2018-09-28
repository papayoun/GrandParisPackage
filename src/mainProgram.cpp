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

SINE_POD generateModel(double theta0, double sigma20, double sd, double thresh = 0.01){
  double newTheta = theta0 + Rcpp::rnorm(1, 0, sd)[0];
  double newSigma2 = sigma20 + Rcpp::rnorm(1, 0, sd)[0];
  if(newSigma2 < thresh)
    newSigma2 = thresh;
  SINE_POD newModel(newTheta, newSigma2);
  return newModel;
}

//' @title EM algortihm
//' @name GEM
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix GEM(const Rcpp::NumericVector& observations,
               const Rcpp::NumericVector& observationTimes,
               const double thetaStart, const double sigma2Start,
               const unsigned int nIterations = 20,
               const unsigned int nModels = 5,
               const double randomWalkParam = 2,
               const unsigned int particleSize = 100,
               const unsigned int densitySampleSize = 30,
               const unsigned int logDensitySampleSize = 30,
               const unsigned int backwardSampleSize = 2,
               const unsigned int backwardSamplingMaxTry = 100000000,
               const unsigned int skeletonSimulationMaxTry = 10000000,
               const bool estimateTheta = true, const bool estimateSigma2 = true){
  Rcpp::NumericMatrix output(nIterations + 1, 2);
  output.fill(sigma2Start);
  output(0, 0) = thetaStart;
  for(unsigned int iter = 1; iter < nIterations + 1; iter++){
    SINE_POD startModel(output(iter - 1, 0), output(iter - 1, 1));
    std::vector<SINE_POD> testedModels(nModels + 1);
    testedModels[0] = startModel;
    for(int i = 1; i < (nModels + 1); i++){
      SINE_POD tmp = generateModel(output(iter - 1, 0), output(iter - 1, 1), 1.0 / iter);
      testedModels[i] = tmp;
    }
    ProposalSINEModel propModel(randomWalkParam, startModel, estimateTheta, estimateSigma2);
    SDEParticleSmoother mySmoother(observations, observationTimes,  propModel,
                                   particleSize, densitySampleSize , logDensitySampleSize,
                                   backwardSampleSize, backwardSamplingMaxTry, skeletonSimulationMaxTry);
    Rcpp::NumericVector E_value = mySmoother.evalEStep(testedModels);
    unsigned int bestModelInd = Rcpp::which_max(E_value);
    output(iter, 0) = testedModels[bestModelInd].getParams()[0];
    output(iter, 1) = testedModels[bestModelInd].getParams()[1];
  }
  return output;
}

