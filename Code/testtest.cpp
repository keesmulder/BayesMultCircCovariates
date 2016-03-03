

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;
// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
int tes() {
   vec out = vec(2);
   out(0) = 3;
   out(1) = 5;
   cout << pow(out, 2);
   return 0;
}
