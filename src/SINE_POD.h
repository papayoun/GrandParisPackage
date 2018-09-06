#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef SINEMODEL_H
#define SINEMODEL_H

//' @title Class of partially diffusion process ruled by a SINE model
//' @name SINE_POD
//' @description POD with the unidimensional SDE sinus model dX_t = sin(X_t - theta)dt + dW_t
//' and the observation model Y_t = X_t + N(0, sigma2) 
//' @export SINE_POD
class SINE_POD{
private:
  SINEModel model;
  double observationVariance;
  double observationDensityUnit(const double& hiddenState, const double& observation) const{
    return GenericFunctions::dnorm(hiddenState, observation, pow(observationVariance, 0.5));};
public:
  SINE_POD(const SINEModel& SM, const double& variance)
    : model(SM), observationVariance(variance){}
  SINE_POD(const double& theta, const double& variance)
    : model(theta), observationVariance(variance){}
  Rcpp::List simulateHSandObs(const double& x0, const Rcpp::NumericVector& observationTimes){
    double standardDev = pow(observationVariance, 0.5);
    Rcpp::NumericVector hiddenStates = model.simulateTrajectory(x0, observationTimes);
    Rcpp::NumericVector observations = hiddenStates + Rcpp::rnorm(observationTimes.size(), 0, standardDev);
    return  Rcpp::List::create(Rcpp::Named("observations") = observations,
                               Rcpp::Named("hiddenStates") = hiddenStates,
                               Rcpp::Named("simulationTimes") = observationTimes);
    };
  SINEModel getModel() const{return model;};
  double getSigma2() const{return observationVariance;};
  void setSigma2(const double newSigma2){observationVariance = newSigma2;};
  
  // Density, log density, gradiend functions
  double unbiasedDensity(const double& x0, const double& xF, 
                         const double& t0, const double& tF, 
                         const unsigned int& sampleSize){
    return model.unbiasedDensityEstimate(x0, xF, t0, tF, sampleSize, false);
  }
  double unbiasedLogDensity(const double& x0, const double& xF, 
                             const double& t0, const double& tF,
                             const unsigned int& sampleSize){
    return model.unbiasedLogDensityEstimate(x0, xF, t0, tF, sampleSize, 1000);
  }
  double unbiasedGradLogDensity(const double& x0, const double& xF, 
                                const double& t0, const double& tF,
                                const unsigned int& sampleSize){
    return model.unbiasedGradLogDensityEstimate(x0, xF, t0, tF, sampleSize, 1000);
  }
  double observationDensity(const double& hiddenState, const double& observation) const{
    return observationDensityUnit(hiddenState, observation);
  };
  Rcpp::NumericVector observationDensity(const Rcpp::NumericVector& hiddenStates, const double& observation) const{
    unsigned int particleSize = hiddenStates.size();
    Rcpp::NumericVector output(particleSize);
    for(unsigned int i = 0; i < particleSize; i++){
      output[i] = observationDensityUnit(hiddenStates[i], observation);
    }
    return output;
  }
  Rcpp::NumericVector gradObservationDensity(const double& hiddenState, const double& observation) const{
    Rcpp::NumericVector output(2); output.fill(0);// The first term corresponds to theta and is null
    // -dnorm()/sigma
    double gaussTerm = observationDensityUnit(hiddenState, observation); 
    double quadrTerm = (hiddenState - observation) * (hiddenState - observation);
    output[1] =  gaussTerm * (quadrTerm / observationVariance - 1) / (2 * observationVariance);
    return output;
  }
  Rcpp::NumericVector gradLogObservationDensity(const double& hiddenState, const double& observation) const{
    Rcpp::NumericVector output(2); output.fill(0);// The first term corresponds to theta and is null
    double quadrTerm = (hiddenState - observation) * (hiddenState - observation);
    output[1] =  (quadrTerm / observationVariance - 1) / (2 * observationVariance);
    return output;
  }
  Rcpp::NumericVector gradLogTransitionDensity(const double& oldParticle, const double& newParticle,
                                         const double& t0, const double& tF, 
                                         const unsigned int& sampleSize, const unsigned int& maximalTries = 1000){
    Rcpp::NumericVector output(2); output.fill(0);// The second term corresponds to Sigma and is null
    output[0] = model.unbiasedGradLogDensityEstimate(oldParticle, newParticle, t0, tF, sampleSize, maximalTries);
    return output;
  }
  Rcpp::NumericVector getParams() const{
    Rcpp::NumericVector output(2);
    output[0] = model.getTheta();
    output[1] = getSigma2();
    return output;
  }
  void setParams(const double& newTheta, const double& newSigma2){
    model.setTheta(newTheta);
    setSigma2(newSigma2);
  }
  void setParams(const Rcpp::NumericVector& newParams){
    model.setTheta(newParams[0]);
    setSigma2(newParams[1]);
  }
};
// Exposes the class to Rcpp
RCPP_MODULE(SINEModel_Module) {
  using namespace Rcpp;
  class_<SINE_POD>("SINE_POD")
    .constructor<double, double>("constructor") // This exposes the default constructor
    .method("rSINE_POD", &SINE_POD::simulateHSandObs, 
  "Needs 2 args, a double (starting point), and a vector of increasing simul times")
    .method("density", &SINE_POD::unbiasedDensity)
    .method("logDensity", &SINE_POD::unbiasedLogDensity)
    .method("gradLogDensity", &SINE_POD::unbiasedGradLogDensity)
  ;
}

#endif