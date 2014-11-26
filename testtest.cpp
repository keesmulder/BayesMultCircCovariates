

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
int tes() {
   arma::vec out = arma::vec(2);
   out(0) = 3;
   out(1) = 5;
   std::cout << pow(out, 2);
   return 0;
}
