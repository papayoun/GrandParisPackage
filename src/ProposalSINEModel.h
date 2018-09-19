#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef SINE_POD_H
#define SINE_POD_H

class ProposalSINEModel{
private:
  double randomWalkStandDev;
  SINE_POD trueModel;
  Rcpp::NumericVector toEstimate;// Vector of 0 and 1 having the same length as number of params
  double aroundObsVariance = 1;
  // In the SINE model case, length 2. Equals 1 if the parameter is to estimate, 0 if it is known
  double getKalmanVariance(double Delta) const{
    double randomWalkVariance = randomWalkStandDev * randomWalkStandDev * Delta;
    return randomWalkVariance * aroundObsVariance / (randomWalkVariance + aroundObsVariance);
  }
  Rcpp::NumericVector getKalmanMeans(const Rcpp::NumericVector& startingParticles, 
                                     const double& futureObservation, 
                                     const double& Delta, 
                                     const double& kalmanVariance,
                                     const bool& pureRandomWalk = false) const{
    Rcpp::NumericVector eulerSchemeMean = startingParticles;
    if(!pureRandomWalk){
      eulerSchemeMean += trueModel.getModel().drift(startingParticles) * Delta;
    }
    Rcpp::NumericVector mean1 = eulerSchemeMean / (randomWalkStandDev * randomWalkStandDev * Delta);
    double mean2 = futureObservation / aroundObsVariance;
    return kalmanVariance * (mean1 + mean2);
  };
public:
  ProposalSINEModel(const double& RW, const SINE_POD& TM, 
                    bool estimTheta = true, bool estimSigma2 = true)
    : randomWalkStandDev(RW) , trueModel(TM){
    Rcpp::NumericVector tmp(2); tmp.fill(0);
    if(estimTheta)
      tmp(0) = 1;
    if(estimSigma2)
      tmp(1) = 1;
    toEstimate = tmp;
  };
  Rcpp::NumericVector simulateInitialParticle(const int& particleSize, 
                                              const double& observation) const{
    return Rcpp::rnorm(particleSize, observation, randomWalkStandDev);};
  Rcpp::NumericVector evalInitialDensity(const Rcpp::NumericVector& particles, 
                                         const double& observation) const{
    return Rcpp::dnorm(particles, observation, randomWalkStandDev);
  };
  Rcpp::NumericVector simulateNextParticle(const Rcpp::NumericVector& oldParticles, 
                                     const double& newObservation, const double& Delta,
                                     const bool& pureRandomWalk = false) const{
    //using euler method
    
    unsigned int particleSize = oldParticles.size();
    Rcpp::NumericVector output(particleSize);
    double kalmanVariance = getKalmanVariance(Delta);
    Rcpp::NumericVector kalmanMeans = getKalmanMeans(oldParticles, newObservation, 
                                                     Delta, kalmanVariance, pureRandomWalk);
    for(unsigned int i = 0; i < particleSize; i++){
      output[i] = Rcpp::rnorm(1, kalmanMeans[i], sqrt(kalmanVariance))[0];
    }  
    return output;
  };
  Rcpp::NumericVector densityNextParticle(const Rcpp::NumericVector& oldParticles,
                                    const Rcpp::NumericVector& newParticles,
                                    const double& newObservation,
                                    const double& Delta) const{
    //using euler method
    Rcpp::NumericVector output(oldParticles.size());
    double kalmanVariance = getKalmanVariance(Delta);
    Rcpp::NumericVector kalmanMeans = getKalmanMeans(oldParticles, newObservation, 
                                                     Delta, kalmanVariance);
    for(unsigned int i = 0; i < oldParticles.size(); i++){
      output[i] = GenericFunctions::dnorm(newParticles[i], 
                                          kalmanMeans[i], sqrt(kalmanVariance));
    }  
    return output;
  }
  double evalTransitionDensity(const double& oldParticle, const double& newParticle,
                               const double& t0, const double& tF, 
                               const unsigned int& sampleSize, const bool& GPE2 = false) const{
    return trueModel.getModel().unbiasedDensityEstimate(oldParticle, newParticle, t0, tF, sampleSize, GPE2);
  };
  Rcpp::NumericVector evalTransitionDensity(const Rcpp::NumericVector& oldParticles, const Rcpp::NumericVector& newParticles,
                                      const double& t0, const double& tF, const unsigned int& sampleSize, const bool& GPE2 = false) const{
    unsigned int particleSize = oldParticles.size();
    Rcpp::NumericVector output(particleSize);
    for(unsigned int i = 0; i < particleSize; i++){
      output[i] = evalTransitionDensity(oldParticles[i], newParticles[i], t0, tF, sampleSize, GPE2);
    }
    return output;
  };
  Rcpp::NumericVector evalGradLogTransitionDensity(const double& oldParticle, const double& newParticle,
                                             const double& t0, const double& tF, 
                                             const unsigned int& sampleSize, const unsigned int& maximalTries = 1000){
    // toEstimate is 0 or 1 valued vector. The gradient is w.r.t to unknown parameterss
    return toEstimate * trueModel.gradLogTransitionDensity(oldParticle, newParticle, 
                                                           t0, tF, sampleSize,
                                                           maximalTries);
  };
  Rcpp::NumericVector evalObservationDensity(const Rcpp::NumericVector& newParticles, const double& observation) const{
    return trueModel.observationDensity(newParticles, observation);
  };
  Rcpp::NumericVector evalGradObservationDensity(const double& hiddenState, const double& observation) const{
    return toEstimate * trueModel.gradObservationDensity(hiddenState, observation);
  }
  Rcpp::NumericVector evalGradLogObservationDensity(const double& newParticle, const double& observation) const{
    return toEstimate * trueModel.gradLogObservationDensity(newParticle, observation);
  };
  Rcpp::NumericVector getFixedGPETerms(const Rcpp::NumericVector& oldParticles, const double& childParticle,
                                 const double& Delta) const{
    unsigned int particleSize = oldParticles.size();
    Rcpp::NumericVector output(particleSize);
    for(unsigned int i = 0; i < particleSize; i++){
      output[i] = trueModel.getModel().fixedGPETerm(oldParticles[i], childParticle, Delta, false);
    }
    return output;
  }
  Rcpp::NumericVector getParams() const{
    return trueModel.getParams();
  }
  void setParams(const double& newTheta, const double& newSigma){
    trueModel.setParams(newTheta, newSigma);
  }
  void setParams(const Rcpp::NumericVector& newParams){
    trueModel.setParams(newParams);
  }
};
#endif
