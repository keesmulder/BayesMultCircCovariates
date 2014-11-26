
 // [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//using namespace arma;

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
arma::vec tes() {
  arma::vec out = arma::vec(1);
//  out.zeros();
  return out;
}



// [[Rcpp::export]]
arma::mat mt(arma::mat X, arma::vec bt) {
  arma::mat out = X;
  out.each_row() *= bt;
  return out;
}

