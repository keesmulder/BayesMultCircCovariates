
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
mat tes() {
  mat A = randu<mat>(5,10);
  return(A);
}

// [[Rcpp::export]]
mat tes2() {
  mat A = zeros<mat>(20, 0);
  mat B = zeros<mat>(0, 20);

  cout << A << endl;
  cout << B << endl;

  return(A * B);
}

// [[Rcpp::export]]
mat tes3() {
  mat A = zeros<mat>(20, 0);
  mat B = zeros<mat>(0, 20);

  cout << A << endl;
  cout << B << endl;

  return((double) B * A);
}
//
// // [[Rcpp::export]]
// double tes4() {
//   mat A = zeros<mat>(20, 0);
//   mat B = zeros<mat>(0, 20);
//
//   cout << A << endl;
//   cout << B << endl;
//
//   return(as_scalar(B * A));
// }
//

