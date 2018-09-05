// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "DebugMethods.h"
#include "GenericFunctions.h"
#include "SINEModel.h"
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
//' @export GaussModLin
class GaussModLin{
private:
  double beta;
public:
  GaussModLin(double beta_):
    beta(beta_) {}
  double getBeta(){return beta;}
};


// Expose the classes
RCPP_MODULE(MyModule) {
  using namespace Rcpp;
  
  class_<GaussModLin>("GaussModLin")
    .constructor<double>("constructor") // This exposes the default constructor
    .method("getBeta", &GaussModLin::getBeta) // This exposes the estim method
  ;
  
}


