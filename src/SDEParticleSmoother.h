#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef ProposalSINEModel_H
#define ProposalSINEModel_H

class SDEParticleSmoother{
private:
  //Attributes
  // constant attributes
  const unsigned int particleSize;
  const Rcpp::IntegerVector particleIndexes;
  const Rcpp::NumericVector observations;
  const Rcpp::NumericVector observationTimes;
  const unsigned int observationSize;
  // to be updated attributes
  Rcpp::NumericVector observationDensityValues;
  Rcpp::NumericMatrix particleSet;
  Rcpp::NumericMatrix particleFilteringWeights;
  Rcpp::NumericMatrix tau_E_X; // matrices of dim particleSize * observationSize
  std::vector<Rcpp::NumericMatrix> tauEStep;// numberModels matrices of dim particleSize * observationSize
  std::vector<Rcpp::NumericMatrix> tauTangentFilter;
  ProposalSINEModel propModel;
  Rcpp::NumericVector zeta1, zeta2;// Conform to paper notations
  double sum_IS_weights;
  double zeta3;
  Rcpp::NumericVector oldGradient;
  Rcpp::NumericVector gradientStep;
  // Simulation settings attributes
  // constant ones
  const unsigned int densitySampleSize;
  const unsigned int logDensitySampleSize;
  const unsigned int backwardSampleSize;
  const unsigned int backwardSamplingMaxTry;
  const unsigned int skeletonSimulationMaxTry;
  const unsigned int numberOfSampledJs;
  // that could change
  unsigned int backwardSamplingCounter;
  Rcpp::IntegerVector  backwardIndexCandidates;
  double densityUpperBound;
  // Methods
  // Generic one const
  Rcpp::NumericVector normWeights(Rcpp::NumericVector unNormedWeights) const{return unNormedWeights / sum(unNormedWeights);};
  // Initialization methodes
  void initializeTauEStep(const unsigned int numberModels){
    std::vector<Rcpp::NumericMatrix>  initial(numberModels);
    for(unsigned int i = 0; i < numberModels; i++){
      Rcpp::NumericMatrix zeroMatrix(particleSize, observationSize); zeroMatrix.fill(0);
      initial[i] = zeroMatrix;
    }
    tauEStep = initial;
  };
  void initializeTauTangentFilter(){
    for(unsigned int i = 0; i < observationSize; i++){
      Rcpp::NumericMatrix zeroMatrix(particleSize, 2); zeroMatrix.fill(0);
      tauTangentFilter[i] = zeroMatrix;
    }
  };
  void initializeTau_E_X(){
    Rcpp::NumericMatrix zeroMatrix(particleSize, observationSize); zeroMatrix.fill(0);
    tau_E_X = zeroMatrix;
  };
  void initializeBackwardSampling(unsigned int ancestorIndex){
    backwardSamplingCounter = 0;
    Rcpp::NumericVector weights = particleFilteringWeights(Rcpp::_, ancestorIndex);
    backwardIndexCandidates = GenericFunctions::sampleReplace(particleIndexes, numberOfSampledJs, weights);
  }
  void setInitalParticles(){
    Rcpp::NumericVector firstParticles = propModel.simulateInitialParticle(particleSize, observations[0]);
    particleSet(Rcpp::_, 0) = firstParticles;
    Rcpp::NumericVector unNormedWeights = propModel.evalInitialDensity(firstParticles, observations[0]);
    particleFilteringWeights(Rcpp::_, 0) = normWeights(unNormedWeights);
  };
  // Propagation method
  void propagateParticles(const unsigned int& ancestorIndex, bool pureRW = false){
    double t0 = observationTimes[ancestorIndex];
    double tF = observationTimes[ancestorIndex + 1];
    double Delta = tF - t0;
    Rcpp::NumericVector currentWeights = particleFilteringWeights(Rcpp::_, ancestorIndex);
    Rcpp::IntegerVector selectedIndexes = GenericFunctions::sampleReplace(particleIndexes, 
                                                                          particleSize,
                                                                          currentWeights);
    Rcpp::NumericVector selectedParticles(particleSize);
    for(unsigned int i =0; i < particleSize; i++){
      selectedParticles[i] = particleSet(selectedIndexes[i], ancestorIndex);
    }
    double futureObs = observations[ancestorIndex + 1];
    Rcpp::NumericVector newParts = propModel.simulateNextParticle(selectedParticles, futureObs, Delta, pureRW);
    particleSet(Rcpp::_, ancestorIndex + 1)   = newParts;
    Rcpp::NumericVector proposalDensityValues = propModel.densityNextParticle(selectedParticles, 
                                                                              newParts, 
                                                                              futureObs, Delta, pureRW);
    Rcpp::NumericVector transitionDensityValues = propModel.evalTransitionDensity(selectedParticles, 
                                                                                  newParts,
                                                                            t0, tF, densitySampleSize, false);// the true is to use GPE2
    observationDensityValues = propModel.evalObservationDensity(newParts, futureObs);
    Rcpp::NumericVector unNormedWeights = transitionDensityValues * observationDensityValues / proposalDensityValues;
    particleFilteringWeights(Rcpp::_, ancestorIndex + 1) = normWeights(unNormedWeights);
  };
  void setDensityUpperBound(const unsigned int& childIndex, const unsigned int& particleIndex){
    double currentChild = particleSet(particleIndex, childIndex);
    Rcpp::NumericVector ancestorParts = particleSet(Rcpp::_, childIndex - 1);
    double delta = observationTimes[childIndex] - observationTimes[childIndex - 1]; 
    Rcpp::NumericVector fixedGPETerms = propModel.getFixedGPETerms(ancestorParts,
                                                                   currentChild, 
                                                                   delta);
    densityUpperBound = max(fixedGPETerms);
  };// By quadratic method
  unsigned int getBackwardIndex(const unsigned int& childIndex, const unsigned int& particleIndex){
    bool acceptedBackwardIndex = false;
    unsigned int tryBSCounter = 0;
    unsigned int ancestorIndexCandidate = -1;
    Rcpp::NumericVector backwardWeights = particleFilteringWeights(Rcpp::_, childIndex - 1);
    while((!acceptedBackwardIndex) & (tryBSCounter < backwardSamplingMaxTry)){
      if(backwardSamplingCounter == numberOfSampledJs){// IF Js values have been tested, try again
        backwardIndexCandidates = GenericFunctions::sampleReplace(particleIndexes, numberOfSampledJs, backwardWeights);
        backwardSamplingCounter = 0;
      }
      ancestorIndexCandidate = backwardIndexCandidates[backwardSamplingCounter];
      double u = Rcpp::runif(1, 0, 1)[0];
      double ancestor = particleSet(ancestorIndexCandidate, childIndex - 1);
      double child = particleSet(particleIndex, childIndex);
      double sampledQEstimate = propModel.evalTransitionDensityUnit(ancestor,//x0
                                                                child ,//xF
                                                                observationTimes[childIndex - 1], //t0
                                                                observationTimes[childIndex], densitySampleSize);
      acceptedBackwardIndex = (u < sampledQEstimate / densityUpperBound);
      tryBSCounter += 1;
      backwardSamplingCounter +=1;
    }
    if(!acceptedBackwardIndex)
      std::cout << "At the observation index " << childIndex << ", no ancestor was accepted for particle" << particleIndex << std::endl;
    return ancestorIndexCandidate;
  };
  unsigned int getBackwardIndex_IS(const unsigned int& childIndex){
    Rcpp::NumericVector backwardWeights = particleFilteringWeights(Rcpp::_,
                                                                   childIndex - 1);
    unsigned int ancestor = GenericFunctions::sampleReplace(particleIndexes,
                                                            1,
                                                            backwardWeights)(0);
    return ancestor;
  };
  void updateTauTracking(const unsigned int& childIndex, 
                         const unsigned int& childParticleIndex,
                      const unsigned int& ancestorParticleIndex, 
                      const bool& tracked_X){
    double h_term = 0;
    if(tracked_X){
      h_term = particleSet(ancestorParticleIndex, childIndex - 1);
    }
    tau_E_X(childParticleIndex, childIndex) += (tau_E_X(ancestorParticleIndex, childIndex - 1) + 
        h_term) / backwardSampleSize;
  };
  void updateTauTracking_IS(const unsigned int& childIndex, 
                         const unsigned int& childParticleIndex,
                         const unsigned int& ancestorParticleIndex,
                         const double IS_weight,
                         const bool& tracked_X){
    double h_term = 0;
    if(tracked_X){
      h_term = particleSet(ancestorParticleIndex, childIndex - 1);
    }
    tau_E_X(childParticleIndex, childIndex) += IS_weight * (tau_E_X(ancestorParticleIndex, childIndex - 1) + h_term);
  };
  void updateTauEStep(const unsigned int& childIndex, const unsigned int& childParticleIndex,
                      const unsigned int& ancestorParticleIndex, const std::vector<SINE_POD>& testedModels){
    for(unsigned int m = 0; m < testedModels.size(); m++){
      SINE_POD model = testedModels[m];
      double sampledLogQ = model.getModel().unbiasedLogDensityEstimate(particleSet(ancestorParticleIndex, childIndex - 1),
                                          particleSet(childParticleIndex, childIndex),
                                          observationTimes[childIndex - 1], observationTimes[childIndex], logDensitySampleSize, skeletonSimulationMaxTry);
      double logObsDensityTerm =  log(model.observationDensity(particleSet(childParticleIndex, childIndex), observations[childIndex]));
      tauEStep[m](childParticleIndex, childIndex) += (tauEStep[m](ancestorParticleIndex, childIndex - 1) + 
                                                      sampledLogQ + logObsDensityTerm) / backwardSampleSize;
    }
  };
  void updateTauEStep_IS(const unsigned int& childIndex, 
                         const unsigned int& childParticleIndex,
                        const unsigned int& ancestorParticleIndex,
                        const double IS_weight,
                        const std::vector<SINE_POD>& testedModels){
    for(unsigned int m = 0; m < testedModels.size(); m++){
      SINE_POD model = testedModels[m];
      double sampledLogQ = model.getModel().unbiasedLogDensityEstimate(particleSet(ancestorParticleIndex, childIndex - 1),
                                          particleSet(childParticleIndex, childIndex),
                                          observationTimes(childIndex - 1),
                                          observationTimes(childIndex),
                                          logDensitySampleSize,
                                          skeletonSimulationMaxTry);
      double logObsDensityTerm =  log(model.observationDensity(particleSet(childParticleIndex, childIndex),
                                                               observations(childIndex)));
      tauEStep[m](childParticleIndex, childIndex) += IS_weight * (tauEStep[m](ancestorParticleIndex, childIndex - 1) +
                                                                  sampledLogQ + logObsDensityTerm);
    }
  };
  void updateTauTangentFilter(const unsigned int& childIndex, const unsigned int& childParticleIndex,
                              const unsigned int& ancestorParticleIndex, bool debug = false){
    Rcpp::NumericVector sampledGradLogQ = propModel.evalGradLogTransitionDensity(particleSet(ancestorParticleIndex, childIndex - 1), 
                                                                                 particleSet(childParticleIndex, childIndex),
                                                                                 observationTimes[childIndex - 1], 
                                                                                 observationTimes[childIndex],
                                                                                 logDensitySampleSize, skeletonSimulationMaxTry);
    Rcpp::NumericVector gradLogObsDensityTerm =  propModel.evalGradLogObservationDensity(particleSet(childParticleIndex, childIndex), observations[childIndex]);
    // gradLogObsDensityTerm.fill(0);
    Rcpp::NumericVector updateTerm = (tauTangentFilter[childIndex - 1](ancestorParticleIndex, Rcpp::_)
      + sampledGradLogQ
      + gradLogObsDensityTerm) / backwardSampleSize;
    tauTangentFilter[childIndex](childParticleIndex, Rcpp::_) = (tauTangentFilter[childIndex](childParticleIndex, Rcpp::_)
        + updateTerm);
  };
  Rcpp::NumericVector getTauBar(const unsigned int&  childIndex){
    return(colMeans(tauTangentFilter[childIndex]));
  }
  void updateZetas(const unsigned int& ancestorIndex){
    zeta1.fill(0);
    zeta2.fill(0);
    Rcpp::NumericVector tauBar = getTauBar(ancestorIndex + 1);
    // DebugMethods::debugprint(tauBar, "tauBar", false);
    double currentObs = observations[ancestorIndex + 1];
    for(unsigned int j = 0; j < particleSize; j++){
      double currentParticle = particleSet(j, ancestorIndex + 1);
      zeta1 += propModel.evalGradObservationDensity(currentParticle, currentObs) / particleSize;
      zeta2 += (tauTangentFilter[ancestorIndex + 1](j, Rcpp::_) - tauBar) * observationDensityValues[j] / particleSize; //* particleFilteringWeights(j, ancestorIndex + 1);
    }
    // DebugMethods::debugprint(zeta1, "zeta1", false);
    // DebugMethods::debugprint(zeta2, "zeta2", false);
    zeta3 = Rcpp::mean(observationDensityValues);
  }
  void updateGradientStep(const Rcpp::NumericVector& currentGradient){
    Rcpp::LogicalVector checkSameSign = GenericFunctions::sameSign(currentGradient, oldGradient);
    for(unsigned int i = 0; i < currentGradient.size(); i++){
      if(oldGradient[i] == 0)
        gradientStep[i] = 0;
      else{
        double scaleFactor = 1.05;
        if(!checkSameSign[i]){
          // std::cout << "Hey, I changed sign";
          scaleFactor = 1 / scaleFactor;
        }
        double gradChangeFactor = 1;
        // double gradChangeFactor = abs(currentGradient[i] / oldGradient[i]);
        gradientStep[i] = scaleFactor * gradientStep[i] * gradChangeFactor;
      }
    }
    // std::cout << "Step value" << gradientStep[0] << std::endl;
    oldGradient = currentGradient;
  }
  Rcpp::NumericVector getGradientUpdate(const unsigned int& ancestorIndex, const unsigned int& updateNumber){
    updateZetas(ancestorIndex);
    Rcpp::NumericVector zetaVector = (zeta1 + zeta2) / zeta3;
    // Rcpp::NumericVector gradientDirection = getDirection(zetaVector);
    if(updateNumber == 0){
      oldGradient = zetaVector;
    }
    // else{
    //   updateGradientStep(zetaVector);// update the gradientStep attribute
    // }
    return gradientStep * zetaVector;
  }
  Rcpp::NumericVector getNewParam(Rcpp::NumericVector& gradientUpdate,
                                  double sigmaThresh = 0.0001) const{
    bool conditionSigma = false;
    Rcpp::NumericVector output(2);
    while(!(conditionSigma)){
      output = propModel.getParams() + gradientUpdate;
      // conditionSigma = (output[1] > 0) & (output[1] < 10);
      conditionSigma = (output[1] >= sigmaThresh);
      if(!(conditionSigma)){
        output[1] = sigmaThresh;
        conditionSigma = true;
      }
      //   
      // if(!(conditionTheta))
      //   gradientUpdate[0] = gradientUpdate[0] / 2;
    }
    return output;
  };
public:
  SDEParticleSmoother(const Rcpp::NumericVector& obs, const Rcpp::NumericVector& times, const ProposalSINEModel& mod,
                      unsigned int N, unsigned int densMonteCarloSize, unsigned int lDensMonteCarloSize,
                      unsigned int NTilde, unsigned int  BSMaxTry, unsigned int skelSimuMaxTry,
                      unsigned int maxNumberOfModels = 10, unsigned int numberJs = 10000)
    : particleSize(N), observations(obs), observationTimes(times), observationSize(obs.size()), propModel(mod),
      densitySampleSize(densMonteCarloSize), logDensitySampleSize(lDensMonteCarloSize), backwardSampleSize(NTilde),
      backwardSamplingMaxTry(BSMaxTry), skeletonSimulationMaxTry(skelSimuMaxTry), numberOfSampledJs(numberJs),
      particleIndexes(Rcpp::seq_len(N)  - 1), particleSet(Rcpp::NumericMatrix(N, obs.size())),
      oldGradient(mod.getParams().size()), gradientStep(mod.getParams().size()),
      zeta1(2), zeta2(2), zeta3(2),// Vectors f gradients
      particleFilteringWeights(Rcpp::NumericMatrix(N, obs.size())),  tauEStep(maxNumberOfModels),
      tauTangentFilter(observationSize), observationDensityValues(particleSize){
    for(unsigned int i = 0; i < maxNumberOfModels; i++){
      Rcpp::NumericMatrix zeroMatrix(particleSize, observationSize); zeroMatrix.fill(0);
      tauEStep[i] = zeroMatrix;
    }
    for(unsigned int i = 0; i < observationSize; i++){
      Rcpp::NumericMatrix zeroMatrix2(particleSize, 2); zeroMatrix2.fill(0);
      tauTangentFilter[i] = zeroMatrix2;
    }
    oldGradient.fill(1);
    gradientStep.fill(1);
  };
  Rcpp::NumericMatrix getParticles() const{return particleSet;};
  Rcpp::NumericMatrix getWeights() const{return particleFilteringWeights;};
  double eval_E_X_tracked(const int tracked_t){
    double output = 0;
    initializeTau_E_X();// Initialize matrix of 0
    setInitalParticles();
    for(int k = 0; k < (observationSize - 1);k++){
      bool update_X = (k == tracked_t);
      propagateParticles(k);
      initializeBackwardSampling(k);// Samples of ancestor index is made here
      for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
        setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
        sum_IS_weights = 0;
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          int chosenAncestorIndex = getBackwardIndex(k + 1, i);
          updateTauTracking(k + 1, i, chosenAncestorIndex, update_X);
          //k + 1 is the time index from which the backward is done, i is the corresponding particle of this generation
        }
      }
    }
    Rcpp::NumericVector lastWeights = particleFilteringWeights(Rcpp::_, observationSize - 1);
    output = sum(lastWeights * tau_E_X(Rcpp::_, observationSize - 1));
    return output;
  }
  double eval_E_X_tracked_IS(const int tracked_t){
    double output = 0;
    initializeTau_E_X();// Initialize matrix of 0
    setInitalParticles();
    for(int k = 0; k < (observationSize - 1);k++){
      bool update_X = (k == tracked_t);
      propagateParticles(k);
      // initializeBackwardSampling(k);// Samples of ancestor index is made here
      Rcpp::NumericVector currentWeights = particleFilteringWeights(Rcpp::_, k);
      for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
        // setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
        sum_IS_weights = 0;
        double curParticle = particleSet(i, k + 1);
        // Choosing ancestoir
        Rcpp::IntegerVector ancestInd = GenericFunctions::sampleReplace(particleIndexes,
                                                                        backwardSampleSize,
                                                                        currentWeights);
        Rcpp::NumericVector ancestPart(backwardSampleSize);
        Rcpp::NumericVector IS_weights(backwardSampleSize);
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          ancestPart(l) = particleSet(ancestInd(l), k);
          IS_weights(l) = propModel.evalTransitionDensityUnit(ancestPart(l), 
                     curParticle,
                     observationTimes(k), 
                     observationTimes(k + 1),
                     densitySampleSize, 
                     false);
          sum_IS_weights += IS_weights(l);
        }
        IS_weights = IS_weights / sum_IS_weights;
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          updateTauTracking_IS(k + 1, i, ancestInd(l), IS_weights(l), update_X);
          //k + 1 is the time index from which the backward is done, 
          //i is the corresponding particle of this generation
        }
        // std::cout << "sum of IS_Weights" << sum_IS_weights << std::endl;
      }
    }
    Rcpp::NumericVector lastWeights = particleFilteringWeights(Rcpp::_, observationSize - 1);
    output = sum(lastWeights * tau_E_X(Rcpp::_, observationSize - 1));
    return output;
  };// end of evalEstep method;
  
  Rcpp::NumericVector evalEStep(const std::vector<SINE_POD>& testedModels){
    unsigned int newParamSize = testedModels.size();
    Rcpp::NumericVector output(newParamSize);
    initializeTauEStep(newParamSize);// Initialize matrix of 0
    setInitalParticles();
    for(int k = 0; k < (observationSize - 1);k++){
      propagateParticles(k);
      initializeBackwardSampling(k);// Samples of ancestor index is made here
      for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
        setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
        sum_IS_weights = 0;
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          int chosenAncestorIndex = getBackwardIndex(k + 1, i);
          updateTauEStep(k + 1, i, chosenAncestorIndex, testedModels);
          //k + 1 is the time index from which the backward is done, i is the corresponding particle of this generation
        }
      }
    }
    Rcpp::NumericVector lastWeights = particleFilteringWeights(Rcpp::_, observationSize - 1);
    for(int m = 0; m < newParamSize; m++){
      output[m] = sum(lastWeights * tauEStep[m](Rcpp::_, observationSize - 1));
    }
    return output;
  }
  Rcpp::NumericVector evalEStep_IS(const std::vector<SINE_POD>& testedModels){
    unsigned int newParamSize = testedModels.size();
    Rcpp::NumericVector output(newParamSize);
    initializeTauEStep(newParamSize);// Initialize matrix of 0
    setInitalParticles();
    for(int k = 0; k < (observationSize - 1);k++){
      propagateParticles(k);
      // initializeBackwardSampling(k);// Samples of ancestor index is made here
      Rcpp::NumericVector currentWeights = particleFilteringWeights(Rcpp::_, k);
      for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
        // setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
        sum_IS_weights = 0;
        double curParticle = particleSet(i, k + 1);
        // Choosing ancestoir
        Rcpp::IntegerVector ancestInd = GenericFunctions::sampleReplace(particleIndexes,
                                                                        backwardSampleSize,
                                                                        currentWeights);
        Rcpp::NumericVector ancestPart(backwardSampleSize);
        Rcpp::NumericVector IS_weights(backwardSampleSize);
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          ancestPart(l) = particleSet(ancestInd(l), k);
          IS_weights(l) = propModel.evalTransitionDensityUnit(ancestPart(l), 
                                                            curParticle,
                                                            observationTimes(k), 
                                                            observationTimes(k + 1),
                                                            densitySampleSize, 
                                                            false);
          sum_IS_weights += IS_weights(l);
        }
        IS_weights = IS_weights / sum_IS_weights;
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          updateTauEStep_IS(k + 1, i, ancestInd(l), IS_weights(l), testedModels);
          //k + 1 is the time index from which the backward is done, 
          //i is the corresponding particle of this generation
        }
        // std::cout << "sum of IS_Weights" << sum_IS_weights << std::endl;
      }
    }
    // for(unsigned int m = 0; m < newParamSize; m++){
    //   std::cout << "m = " << m << std::endl;
    //   Rcpp::NumericMatrix taus = tauEStep[m];
    //   DebugMethods::debugprint(taus, "taus_MC");
    // }
    Rcpp::NumericVector lastWeights = particleFilteringWeights(Rcpp::_, observationSize - 1);
    for(int m = 0; m < newParamSize; m++){
      output[m] = sum(lastWeights * tauEStep[m](Rcpp::_, observationSize - 1));
    }
    return output;
  };// end of evalEstep method;
  Rcpp::NumericMatrix tangentFilterEstimation(const Rcpp::LogicalVector& updateOrders, 
                                              Rcpp::NumericMatrix& gradientSteps){
    Rcpp::NumericMatrix output(observationSize, 2);
    output(0, Rcpp::_) = propModel.getParams();
    initializeTauTangentFilter();// Initialize matrix of 0
    setInitalParticles();
    unsigned int updateCounter = 0;
    for(int k = 0; k < (observationSize - 1); k++){
      // std::cout << "iteration " << k << std::endl;
      // DebugMethods::debugprint(propModel.getParams(), "Current Param", false);
      gradientStep = gradientSteps(k, Rcpp::_);
      propagateParticles(k);
      initializeBackwardSampling(k);// Samples of ancestor index is made here
      for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
        setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
        // std::cout << "particle " << i << std::endl; 
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          int chosenAncestorIndex = getBackwardIndex(k + 1, i);
          updateTauTangentFilter(k + 1, i, chosenAncestorIndex, updateOrders[k]);
          //k + 1 is the time index from which the backward is done, i is the corresponding particle of this generation
        }
      }
      if(updateOrders[k]){
        // std::cout << "Iteration " << k << std::endl;
        Rcpp::NumericVector newGradientUpdate = getGradientUpdate(k, updateCounter);
        // DebugMethods::debugprint(newGradientUpdate, "Gradient", false);
        Rcpp::NumericVector newParams = getNewParam(newGradientUpdate);
        propModel.setParams(newParams);
        updateCounter += 1;
      }
      output(k + 1, Rcpp::_) = propModel.getParams();
    }
    return output;
  };// end of evalEstep method;
};
#endif