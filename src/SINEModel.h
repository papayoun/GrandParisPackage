#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef GenericFunctions_H
#define GenericFunctions_H

class SINEModel{  
private:
  ////////////// ATTRIBUTES /////////////////////////////////
  // The drift parameter
  double theta;
  // Useful bounds
  double psiLowerBound = -0.5;
  double psiUpperBound = 0.625;
  double phiUpperBound = psiUpperBound - psiLowerBound;
  double potentialUpperBound = 1;
  // An attribute skeleton for the simulation
  Rcpp::NumericVector skeleton;
  
  ////////////// METHODS /////////////////////////////////
  // Model methods
  // Potenial of a unique double
  double potentialUnit(const double& x) const {return -cos(x - theta);}
  // Drift of a unique double
  double driftUnit(const double& x) const{return sin(x -theta);};
  // Laplacian 
  double laplaciantUnit(const double& x) const{return cos(x -theta);};
  // Psi, i.e. 0.5 * (norm(Drift) + Trace(Hessian(Potential)))
  double psiUnit(const double& x) const{
    return 0.5  * (driftUnit(x) * driftUnit(x) + laplaciantUnit(x));
  };
  // Phi, i.e. Psi - LowerBound(Psi)
  double phiUnit(const double& x) const{return psiUnit(x) - psiLowerBound;};
  // Derivative of Phi w.r.t. theta
  // It is a double here, generally speaking should be a vector
  double gradThetaPhi(const double& x) const{
    return - 0.5 * sin(x - theta) * (-1 + 2 * cos(x - theta));
  };
  // Derivative of Potential w.r.t. theta
  // It is a double here, generally speaking should be a vector
  double gradThetaPotential(const double& x) const{return - sin(x - theta);};
  // Derivative of LowerBound(Psi) w.r.t. theta
  // It is a double here, generally speaking should be a vector
  double gradPsiLowerBound(const double& th) const{return 0.0;};
  // Non random term of the estimator I_{k+1, theta} as written in the paper
  // It is a double here, generally speaking should be a vector
  double fixedTermLogDensityEstimate(const double& x0, const double& y, 
                                     const double& Delta) const{
    double gradDeltaPot = gradThetaPotential(y) - gradThetaPotential(x0);
    return (gradDeltaPot - gradPsiLowerBound(theta) * Delta); 
  };
  // Simulation methods
  void initiateSkeleton(Rcpp::NumericVector skeletonTimes){
    Rcpp::NumericVector tmp(skeletonTimes.size()); 
    skeleton = tmp;
  }
  Rcpp::NumericVector simulSkelTimes(const unsigned int& skeletonSize, 
                                       const double& startingTime,
                                       const double& endingTime) const{
    Rcpp::NumericVector output = Rcpp::runif(skeletonSize, startingTime,
                                             endingTime);
    std::sort(output.begin(), output.end());
    return output;
  };
  // Do the skeleton is accepted
  // all skeleton or one without begin, min, end?
  bool skeletonAcceptance(const double& x0, const double& xF, const double& t0,
                          const double& tF, 
                          const Rcpp::NumericVector& skeletonTimes, 
                          bool saveSkeleton = false){
    unsigned int skeletonSize = skeletonTimes.size();
    double startingTime = t0;
    double startingPos = x0;
    bool acceptedSkeleton = true;
    unsigned int skeletonIndex = 0;// where am I in the skeleton?
    if(saveSkeleton)
      initiateSkeleton(skeletonTimes);
    while(acceptedSkeleton & (skeletonIndex < skeletonSize)){
      double skelValue = GenericFunctions::brownianBridgeSimulation(skeletonTimes(skeletonIndex),
                                                                    startingPos, xF,
                                                                    startingTime, tF);
      double Upsilon = Rcpp::runif(1, 0, phiUpperBound)[0];
      acceptedSkeleton = (phiUnit(skelValue) < Upsilon);
      // Update startingPos and startingTime
      startingPos = skelValue;
      startingTime = skeletonTimes(skeletonIndex);
      if(saveSkeleton)
        skeleton[skeletonIndex] = skelValue;
      skeletonIndex += 1;
    }
    return acceptedSkeleton;
  }
  // Exact conditional simulation at a random time
  double randomTimeConditionalExactSim(const double& x0, const double& xF, 
                                       const double& t0, const double& tF,
                                       const unsigned int maximalTries = 1000000){
    int tryNumber = 0;// reinitialize number of trials
    bool acceptedSkeleton = false;
    double simulatedXu;
    while((tryNumber < maximalTries) & (!acceptedSkeleton)){
      tryNumber += 1;
      double Delta = tF - t0;
      double skeletonSize = Rcpp::rpois(1, Delta * phiUpperBound)[0];
      if(skeletonSize == 0){
        acceptedSkeleton = true;
        double u = Rcpp::runif(1,t0, tF)[0];// random simulation time
        simulatedXu = GenericFunctions::brownianBridgeSimulation(u, x0, xF, t0, tF);
      } 
      else{
        Rcpp::NumericVector skeletonTimes = simulSkelTimes(skeletonSize, t0, tF);
        acceptedSkeleton = skeletonAcceptance(x0, xF, t0, tF, skeletonTimes, true);
        // The previous function returns a bool, but also updates the skeleton attribute with
        // tested skeleton values
        if(acceptedSkeleton){
          double u = Rcpp::runif(1, t0, tF)(0);// random simulation time
          // The skeleton do not include x0 and xF, a three possibilities treatment
          if(u < skeletonTimes[0]){// u is before the first skeleton time
            simulatedXu = GenericFunctions::brownianBridgeSimulation(u, x0, 
                                                                     skeleton[0], t0, 
                                                                     skeletonTimes[0]);
          } 
          else if(u > skeletonTimes[skeletonSize - 1]){// u is after the last skeleton time
            simulatedXu = GenericFunctions::brownianBridgeSimulation(u, skeleton[skeletonSize - 1],
                                                                     xF, skeletonTimes[skeletonSize - 1], tF);
          } 
          else{// u is within the skeleton times
            unsigned int segmentIndex = GenericFunctions::findInterval(u, skeletonTimes);
            simulatedXu = GenericFunctions::brownianBridgeSimulation(u, skeleton[segmentIndex - 1], 
                                                                     skeleton[segmentIndex],
                                                                     skeletonTimes[segmentIndex - 1],
                                                                     skeletonTimes[segmentIndex]);
          }
        }
      }// end if skeletonSize > 0
    }
    if(!acceptedSkeleton)
      std::cout << "Please note that a problem occured, no skeleton was accepted" << std::endl;
    return simulatedXu;
  }
  // Performing unconditional exact simulation, with a normal proposal
  double unconditionalExactSim(const double& x0, const double& t0, const double& tF,
                               const unsigned int maximalTries = 1000000){
    double Delta = tF - t0;
    double standardDev = sqrt(Delta);
    double candidate;
    bool acceptedCandidate = false;
    int tryNumber = 0;
    // The simulation is done by proposing candidates from a basic randomWalk
    // See thesis of Sermaidis, section 3.2.2 for details
    while((tryNumber < maximalTries) & (!acceptedCandidate)){
      // Simulating terminal point y according to h_t density
      tryNumber += 1;// Updating number of tries
      candidate = x0 + Rcpp::rnorm(1, 0, standardDev)[0];// Proposing a candidate
      double candidatePotential = potentialUnit(candidate);
      double acceptanceRatio = exp(candidatePotential - potentialUpperBound);
      // Rejection sampling procedure
      double u = Rcpp::runif(1, 0, 1)[0];
      if(u < acceptanceRatio){
        //simulation of the minimum
        bool acceptedSkeleton =false;
        // Homogeneous Poisson process generation
        int skeletonSize = Rcpp::rpois(1, Delta * phiUpperBound)[0];
        if(skeletonSize == 0){
          // Easiest case
          acceptedSkeleton = true;
        }
        else{
          Rcpp::NumericVector skeletonTimes = simulSkelTimes(skeletonSize, t0, tF);
          acceptedSkeleton = skeletonAcceptance(x0, candidate, t0, tF, 
                                                skeletonTimes, false);
        }
        acceptedCandidate = acceptedSkeleton;
      }
    }
    // The program stops if no candidate was accepted
    // assert(acceptedCandidate && "No candidate was accepted, increase maximalTries or decrease (tF- t0)");
    if(!acceptedCandidate)
      std::cout << "BAD SIMULATION";
    return candidate;
  };
  double getGPE2GammaTerm(const double& x0, const double& y, const double& Delta) const{
    double trigoTerm = 0.0;
    if(x0 == y){
      trigoTerm = cos(x0) - cos(2 * x0);// Use of the derivative as limit when y = x
    }
    else{
      trigoTerm = (sin(y) * (1 - cos(y)) - sin(x0) * (1 - cos(x0))) / (y - x0);
    }
    return Delta * (psiUpperBound - 0.25 * (1 + trigoTerm));
  }
public:
  // Construction method
  SINEModel(const double& hiddenStateParam)
    : theta(hiddenStateParam){};
  // getter
  double getTheta() const{return theta;};
  // setter
  void setTheta(const double newTheta){ theta = newTheta;};
  // Compute the potential for a vector
  Rcpp::NumericVector potential(const Rcpp::NumericVector& x) const{
    unsigned int n = x.size();
    Rcpp::NumericVector output(n);
    for(unsigned int i = 0; i < x.size(); i ++){
      output(i) = potentialUnit(x(i));
    }
    return output;
  };
  // double driftFun(const double& x) const{return sin(x - theta);};
  Rcpp::NumericVector drift(const Rcpp::NumericVector& x) const{
    unsigned int n = x.size();
    Rcpp::NumericVector output(n);
    for(unsigned int i = 0; i < x.size(); i ++){
      output(i) = driftUnit(x(i));
    }
    return output;
  };
  Rcpp::NumericVector psi(const Rcpp::NumericVector& x) const{
    unsigned int n = x.size();
    Rcpp::NumericVector output(n);
    for(unsigned int i = 0; i < x.size(); i ++){
      output(i) = psiUnit(x(i));
    }
    return output;
  };
  double fixedGPETerm(const double& x0, const double& y, const double& delta, bool logValue = false) const{
    double logGaussianTerm = - 0.5 * (log(2 * M_PI * delta)  + (y - x0) * (y - x0) / delta);
    double logModelTerm = potentialUnit(y) - potentialUnit(x0) - delta * psiLowerBound;
    if(logValue)
      return logGaussianTerm + logModelTerm;
    else
      return exp(logGaussianTerm + logModelTerm);
  };
  Rcpp::NumericVector simulateTrajectory(const double& x0,
                                         const Rcpp::NumericVector& simulationTimes,
                                         const unsigned int& maximalTries = 100000){
    unsigned int simulationSize = simulationTimes.size();
    Rcpp::NumericVector output(simulationSize); output.fill(x0);
    for(unsigned int i = 1; i < simulationSize; i++){
      output(i) = unconditionalExactSim(output(i - 1),
             simulationTimes(i - 1),
             simulationTimes(i), maximalTries);
    }
    return output;
  };// end of simulateTrajectory method
  //The rho_{\Delta_k} of Gloaguen et al. 2018
  
  double unbiasedDensityEstimate(const double& x0, const double& xF, 
                                 const double& t0, const double& tF, 
                                 const unsigned int& sampleSize, const bool& GPE2 = false) const{
    double Delta = tF-t0;
    double randomGPETerm = 0.0;
    for(unsigned int i = 0; i < sampleSize; i++){
      unsigned int skeletonSize;
      double logRandomProduct = 0;
      if(!GPE2){// then it is GPE 1
        skeletonSize = Rcpp::rpois(1, Delta * phiUpperBound)[0];//GPE1 uses a Poisson
        logRandomProduct -= skeletonSize * log(phiUpperBound);// Substracting the log normalization term
      }
      else{
        logRandomProduct += Delta * psiLowerBound; // This term is in the fixed GPE term, and must be
        // removed in case of the use of GPE2
        double betaTerm = 10; //arbitrary, has low influence
        double gammaTerm = getGPE2GammaTerm(x0, xF, Delta);//Gamma is the expectation
        skeletonSize = Rcpp::rnbinom_mu(1, betaTerm, gammaTerm)[0];
        logRandomProduct -= (psiUpperBound * Delta - skeletonSize * log(Delta)
                               + dnbinom_mu(skeletonSize, betaTerm, gammaTerm, true));
        for(unsigned int j = 0; j < skeletonSize; j++){
          logRandomProduct -= log(j + 1);//-log(j+1) is for factorial skeletonSize
        }
      }
      if(skeletonSize > 0){
        Rcpp::NumericVector skeletonTimes = simulSkelTimes(skeletonSize, t0, tF);
        double startingTime = t0;
        double startingPos = x0;
        unsigned int skeletonIndex = 0;
        while(skeletonIndex < skeletonSize){
          double skeletonValue = GenericFunctions::brownianBridgeSimulation(skeletonTimes(skeletonIndex),
                                                                            startingPos, xF,
                                                                            startingTime, tF);
          logRandomProduct += log(psiUpperBound - psiUnit(skeletonValue));
          startingTime = skeletonTimes(skeletonIndex);
          startingPos = skeletonValue;
          skeletonIndex += 1;
        }
      }
      randomGPETerm += exp(logRandomProduct) / sampleSize;
    }
    return fixedGPETerm(x0, xF, Delta) * randomGPETerm;
  };// end of unbiasedDensityEstimate method
  // Log density unbiased estimate, as defined by Olsson et al 2011
  double unbiasedLogDensityEstimate(const double& x0, const double& xF, 
                                    const double& t0, const double& tF,
                                    const unsigned int& sampleSize, 
                                    const unsigned int& maximalTries = 1000){
    double Delta = tF-t0;
    // We use here the equation B3 of  Olsson et al 2011, appendix B
    double fixedPart = fixedGPETerm(x0, xF, Delta, true);
    double randomPart = 0.0;
    for(unsigned int i =0; i < sampleSize; i++){
      double simulatedXu = randomTimeConditionalExactSim(x0, xF, t0, tF, maximalTries);
      randomPart += - Delta * phiUnit(simulatedXu) / sampleSize;
    }// end of loop over sampleSize
    return fixedPart + randomPart;
  };// end of unbiasedLogDensity method;
  double unbiasedGradLogDensityEstimate(const double& x0, const double& xF, 
                                        const double& t0, const double& tF,
                                        const unsigned int& sampleSize, 
                                        const unsigned int& maximalTries = 1000){
    double Delta = tF-t0;
    // We here the equation B3 of  Olsson et al 2011, appendix B
    double fixedPart = fixedTermLogDensityEstimate(x0, xF, Delta);
    double randomPart = 0.0;
    for(unsigned int i =0; i < sampleSize; i++){
      double simulatedXu = randomTimeConditionalExactSim(x0, xF, t0, tF, maximalTries);
      randomPart -= gradThetaPhi(simulatedXu) / sampleSize;
    }// end of loop over sampleSize
    return fixedPart + Delta * randomPart;
  };// end of unbiasedGradLogDensity method;
};
#endif