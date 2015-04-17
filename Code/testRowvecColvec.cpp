

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <iostream>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
bool vectest() {
   vec a = randu<vec>(5);
   cout << a;
   cout << "a colvec:" << a.is_colvec() << endl;
   cout << "a rowvec:" << a.is_rowvec() << endl;
   cout << "a vec:" << a.is_vec() << endl;

   rowvec at = a.t();

   cout << "at colvec:" << at.is_colvec() << endl;
   cout << "at rowvec:" << at.is_rowvec() << endl;
   cout << "at vec:" << at.is_vec() << endl;



   return a.is_colvec();
}
