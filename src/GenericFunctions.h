#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#include <assert.h>

class GenericFunctions{
public:
  static Rcpp::IntegerVector findInterval(Rcpp::NumericVector x, Rcpp::NumericVector breaks) {
    // equivalent of findInterval function
    Rcpp::IntegerVector out(x.size());
    Rcpp::NumericVector::iterator it, pos;
    Rcpp::IntegerVector::iterator out_it;
    for(it = x.begin(), out_it = out.begin(); it != x.end();++it, ++out_it) {
      pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
      *out_it = std::distance(breaks.begin(), pos);
    }
    return out;
  }
  static unsigned int findInterval(double x, Rcpp::NumericVector breaks) {
    // equivalent of findInterval function, for a single value
    Rcpp::NumericVector vec(1);
    vec(0) = x;
    return findInterval(vec,breaks)[0];
  }
  // dnbinom function, but for a single int
  static double dnbinom(int x, double size, double p, bool lg){
    Rcpp::IntegerVector xtmp(1); xtmp(0) = x;
    return Rcpp::dnbinom(xtmp, size, p, lg)[0];
  }
  static double dnorm(double x, double mean,//Do not exist when "x" is a double
                      double sigma, bool log_val=false){// sigma is the standard deviation
    double log2pi = 1.8378770664093454835606594728112352797227949472755668;
    double logretval = -0.5*(log2pi + 2*log(sigma)
                               +(x-mean)*(x-mean)/(sigma*sigma));
    if(log_val){
      return logretval;
    } else {
      return exp(logretval);}
  }
  static Rcpp::LogicalVector sameSign(const Rcpp::NumericVector& x,
                                      const Rcpp::NumericVector& y){
    int n = x.size();
    Rcpp::LogicalVector output(n), signX(n), signY(n);
    for(int i = 0; i < n; i++){
      signX[i] = (x[i] >= 0);
      signY[i] = (y[i] >= 0);
      output[i] = !(signX[i] xor signY[i]);   
    }
    return output;
  }
  // Sampling with replacement
  static Rcpp::NumericVector sampleReplace(const Rcpp::NumericVector& x,
                                           const int& size,
                                           const Rcpp::NumericVector& probs){
    int nx=x.size();int np=probs.size();
    if(nx!=np){
      Rcpp::stop("Error, x and probs must have same length");
    }
    Rcpp::NumericVector cumprob(nx + 1);
    cumprob.fill(0);
    for(int i = 1; i < (nx + 1); i++){
      cumprob(i) = cumprob(i - 1) + probs(i - 1);
    }
    Rcpp::NumericVector us = Rcpp::runif(size, 0, cumprob(nx));
    Rcpp::IntegerVector inds = findInterval(us, cumprob) - 1;
    return x[inds];
  }
  static Rcpp::IntegerVector sampleReplace(const Rcpp::IntegerVector& x,
                                     const int& size,
                                     const Rcpp::NumericVector& probs){
    int nx=x.size();int np=probs.size();
    if(nx!=np){
      Rcpp::stop("Error, x and probs must have same length");
    }
    Rcpp::NumericVector cumprob(nx+1);
    cumprob.fill(0);
    for(int i=1;i<(nx+1);i++){
      cumprob(i) = cumprob(i-1) + probs(i-1);
    }
    Rcpp::NumericVector us = Rcpp::runif(size, 0, cumprob(nx));
    Rcpp::IntegerVector inds = findInterval(us, cumprob) - 1;
    return x[inds];
  }
  static double brownianBridgeSimulation(const double& simulationTime, 
                                         const double& startingPos, 
                                         const double& endingPos, 
                                         const double& startingTime,
                                         const double& endingTime){
    double totalTime = endingTime - startingTime;
    double deltaPos = endingPos - startingPos;
    double deltaTimeStart = simulationTime - startingTime;
    double deltaTimeEnd = (endingTime - simulationTime);
    double mean = startingPos + deltaTimeStart / totalTime * deltaPos;
    double variance =  deltaTimeEnd * deltaTimeStart / totalTime;
    return Rcpp::rnorm(1, mean , sqrt(variance))[0];
  };
  static Rcpp::NumericVector rInvGaussian(unsigned int sampleSize, double mu, 
                                   double lambda){
    // Sermaidis notations
    Rcpp::NumericVector output(sampleSize);
    double R; double J; double u2;
    for(unsigned int i = 0; i < sampleSize; i++){
      R = Rcpp::rchisq(1, 1)[0];
      J = mu + (mu * R - sqrt(4 * mu * lambda * R + mu * mu * R * R)) / lambda * mu * 0.5;
      u2 = Rcpp::runif(1, 0, 1)[0];
      if(u2 < mu / (mu + J))
        output[i] =J;
      else
        output[i] = mu * mu / J;
    }
    return output;
  }
};









