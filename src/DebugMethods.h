// Class of printing methods for object use in the main package

#include "RcppArmadillo.h"

class DebugMethods{// Useful functions for debugging
  public:
    static void debugprint(Rcpp::IntegerVector x, std::string name = ""){
      int n = x.size();
      std::cout << "Vecteur"<< name<< ": ";
      for(int i = 0; i < n-1; i++){
        std::cout << x(i) << ", ";
      }
      std::cout << x(n - 1) << std::endl;
    }
    static void debugprint(Rcpp::NumericMatrix x, std::string name = ""){
      int nr = x.nrow();
      int nc = x.ncol();
      std::cout << "Matrix "<< name<< ": " <<std::endl;
      for(int i = 0; i < nr; i++){
        std::cout << " | ";
        for(int j = 0; j < nc - 1; j++)
          std::cout << x(i, j) << " | ";
        std::cout << x(i, nc - 1) << " |" << std::endl;
      }
    }
    static void debugprint(Rcpp::LogicalVector x, std::string name = ""){
      int n = x.size();
      std::cout << "Vecteur"<< name<< ": ";
      for(int i = 0; i < n-1; i++){
        std::cout << x(i) << ", ";
      }
      std::cout << x(n - 1) << std::endl;
    }
    static void debugprint(Rcpp::NumericVector x, std::string name = "", bool summary = true){
      int n = x.size();
      std::cout << "Vecteur"<< name<< ": ";
      if(summary){
        std::cout << "mean "<< mean(x) << std::endl;
        std::cout << "var "<< var(x) << std::endl;
      }
      else{
        std::cout << "( ";
        for(int j = 0; j < n - 1; j++)
          std::cout << x(j) << ", ";
        std::cout << x(n - 1) << " )" << std::endl;
      }
    }
};


