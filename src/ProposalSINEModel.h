#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef SINE_POD_H
#define SINE_POD_H

//' @title Proposition model
//' @name ProposalSINEModel
//' @description POD with the unidimensional SDE sinus model dX_t = sin(X_t - theta)dt + dW_t
//' and the observation model Y_t = X_t + N(0, sigma2), the proposal model is a random walk 
//' @export ProposalSINEModel
class ProposalSINEModel{
private:
  double randomWalkStandDev;
  SINE_POD trueModel;
  Rcpp::NumericVector toEstimate;// Vector of 0 and 1 having the same length as number of params
  double aroundObsVariance ;
  void updateAroundObsVariance(){
    aroundObsVariance = trueModel.getSigma2();
  }
  // In the SINE model case, length 2. Equals 1 if the parameter is to estimate, 0 if it is known
  double getKalmanVariance(double Delta) const{
    double randomWalkVariance = randomWalkStandDev * randomWalkStandDev * Delta;
    return randomWalkVariance * aroundObsVariance / (randomWalkVariance + aroundObsVariance);
  }
  Rcpp::NumericVector getKalmanMeans(const Rcpp::NumericVector& startingParticles, 
                                     double futureObservation, 
                                     double Delta, 
                                     double kalmanVariance,
                                     bool pureRandomWalk) const{
    Rcpp::NumericVector eulerSchemeMean = Rcpp::clone(startingParticles);
    if(!pureRandomWalk){
      Rcpp::NumericVector addTerm = trueModel.getModel().drift(startingParticles) * Delta;
      eulerSchemeMean += addTerm;
    }
    Rcpp::NumericVector mean1 = eulerSchemeMean / (randomWalkStandDev * randomWalkStandDev * Delta);
    double mean2 = futureObservation / aroundObsVariance;
    Rcpp::NumericVector output = kalmanVariance * (mean1 + mean2);
    return output;
  };
public:
  ProposalSINEModel(double RW, const SINE_POD& TM, 
                    bool estimTheta = true, bool estimSigma2 = true)
    : randomWalkStandDev(RW) , trueModel(TM){
    Rcpp::NumericVector tmp(2); tmp.fill(0);
    if(estimTheta)
      tmp(0) = 1;
    if(estimSigma2)
      tmp(1) = 1;
    toEstimate = tmp;
    updateAroundObsVariance();
  };
  ProposalSINEModel(double RW, double theta, double sigma2, 
                    bool estimTheta = true, bool estimSigma2 = true)
    : randomWalkStandDev(RW) , trueModel(theta, sigma2){
    Rcpp::NumericVector tmp(2); tmp.fill(0);
    if(estimTheta)
      tmp(0) = 1;
    if(estimSigma2)
      tmp(1) = 1;
    toEstimate = tmp;
    updateAroundObsVariance();
  };
  Rcpp::NumericVector simulateInitialParticle(unsigned int particleSize, 
                                              double observation) const{
    return Rcpp::rnorm(particleSize, observation, randomWalkStandDev);};
  Rcpp::NumericVector evalInitialDensity(const Rcpp::NumericVector& particles, 
                                         double observation) const{
    return Rcpp::dnorm(particles, observation, randomWalkStandDev);
  };
  Rcpp::NumericVector simulateNextParticle(const Rcpp::NumericVector& oldParticles, 
                                     double newObservation, double Delta,
                                     bool pureRandomWalk) const{
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
                                    double newObservation,
                                    double Delta,
                                    bool pureRandomWalk) const{
    //using euler method
    Rcpp::NumericVector output(oldParticles.size());
    double kalmanVariance = getKalmanVariance(Delta);
    Rcpp::NumericVector kalmanMeans = getKalmanMeans(oldParticles, newObservation, 
                                                     Delta, kalmanVariance, pureRandomWalk);
    for(unsigned int i = 0; i < oldParticles.size(); i++){
      output[i] = GenericFunctions::dnorm(newParticles[i], 
                                          kalmanMeans[i], sqrt(kalmanVariance));
    }  
    return output;
  }
  double evalTransitionDensityUnit(double oldParticle, double newParticle,
                               double t0, double tF, 
                               unsigned int sampleSize, bool GPE2 = false) const{
    return trueModel.getModel().unbiasedDensityEstimate(oldParticle, newParticle, t0, tF, sampleSize, GPE2);
  };
  Rcpp::NumericVector evalTransitionDensity(const Rcpp::NumericVector& oldParticles, const Rcpp::NumericVector& newParticles,
                                      double t0, double tF, unsigned int sampleSize, bool GPE2 = false){
    unsigned int particleSize = oldParticles.size();
    Rcpp::NumericVector output(particleSize);
    for(unsigned int i = 0; i < particleSize; i++){
      output[i] = evalTransitionDensityUnit(oldParticles[i], newParticles[i], t0, tF, sampleSize, GPE2);
    }
    return output;
  };
  Rcpp::NumericVector evalGradLogTransitionDensity(double oldParticle, double newParticle,
                                             double t0, double tF, 
                                             unsigned int sampleSize, unsigned int maximalTries = 1000){
    // toEstimate is 0 or 1 valued vector. The gradient is w.r.t to unknown parameterss
    return toEstimate * trueModel.gradLogTransitionDensity(oldParticle, newParticle, 
                                                           t0, tF, sampleSize,
                                                           maximalTries);
  };
  Rcpp::NumericVector evalObservationDensity(const Rcpp::NumericVector& newParticles, double observation) const{
    return trueModel.observationDensity(newParticles, observation);
  };
  Rcpp::NumericVector evalGradObservationDensity(double hiddenState, double observation) const{
    return toEstimate * trueModel.gradObservationDensity(hiddenState, observation);
  }
  Rcpp::NumericVector evalGradLogObservationDensity(double newParticle, double observation) const{
    return toEstimate * trueModel.gradLogObservationDensity(newParticle, observation);
  };
  Rcpp::NumericVector evalGradObsDensityGeometricMean(const Rcpp::NumericVector& particles,
                                                      double observation){
    return toEstimate * trueModel.gradObsDensityGeometricMean(particles, observation);
  }
  Rcpp::NumericVector getFixedGPETerms(const Rcpp::NumericVector& oldParticles, double childParticle,
                                 double Delta) const{
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
  void setParams(double newTheta, double newSigma){
    trueModel.setParams(newTheta, newSigma);
    updateAroundObsVariance();
  }
  void setParams(const Rcpp::NumericVector& newParams){
    trueModel.setParams(newParams);
    updateAroundObsVariance();
  }
};
// Exposes the class to Rcpp
RCPP_MODULE(ProposalModel_Module) {
  using namespace Rcpp;
  class_<ProposalSINEModel>("ProposalSINEModel")
    .constructor<double, double, double, bool, bool>("constructor") // This exposes the default constructor
    .method( "transDens" , (Rcpp::NumericVector (ProposalSINEModel::*)(const Rcpp::NumericVector&, const Rcpp::NumericVector&, double,double,unsigned int,bool) )( 
  &ProposalSINEModel::evalTransitionDensity)  )
    .method( "obsDens", &ProposalSINEModel::evalObservationDensity)
    .method("propDens", &ProposalSINEModel::densityNextParticle)
  ;
}
#endif