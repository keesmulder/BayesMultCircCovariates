
 // [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <iostream>
#include <math.h>

using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar



// [[Rcpp::export]]
vec tes() {
  rowvec rv = rowvec(2);
  rv(0)  = 3;
  rv(1)  = 2;

  colvec cv = colvec(2);
  cv(0)  = 30;
  cv(1)  = 23;

  return rv + cv;
}


