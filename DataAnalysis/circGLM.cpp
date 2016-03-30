/*
# ----------------------------------------------------------
# circGLMc.cpp
# Runs an MCMC sampler for circular data.
#
# Kees Tim Mulder
#
#

 TODO
 Check if quantiles are obtained properly.

#
# This work was supported by a Vidi grant awarded to I. Klugkist from the
# Dutch Organization for Scientific research (NWO 452-12-010).
# ----------------------------------------------------------
 */



// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <iostream>
#include <math.h>
#include <typeinfo>


using namespace Rcpp;
using namespace std;
using namespace arma;

const double pi = boost::math::constants::pi<double>();

// [[Rcpp::export]]
double cvariance(vec x) {
  return var(x);
}




// [[Rcpp::export]]
vec atanLF(vec x, double r) {
  // arctangent-family link function, spreading range r, where r=1
  // corresponds to the semicircle, and r=2 is most common.
  return r * atan(x);
}
// [[Rcpp::export]]
double atanLFdouble(double x, double r) {
  // arctangent-family link function, spreading range r, where r=1
  // corresponds to the semicircle, and r=2 is most common.
  return r * atan(x);
}

// [[Rcpp::export]]
vec invAtanLF(vec x, double r) {
  // arctangent-family inverse link function, spreading range r, where r=1
  // corresponds to the semicircle, and r=2 is most common.
  return tan(x/r);
}

// [[Rcpp::export]]
double invAtanLFdouble(double x, double r) {
  // arctangent-family inverse link function, spreading range r, where r=1
  // corresponds to the semicircle, and r=2 is most common.
  return tan(x/r);
}


// [[Rcpp::export]]
NumericVector rvmc(int n, double mu, double kp) {
  /* FUNCTION rvmc -------------------------------------------
   Generate random variates from the von Mises distribution.

   n:      The number of random variates required.
   mu:     The required mean direction, mu.
   kp:     The required concentration, kappa.

   Returns: A vector of length n containing VM random variates.
   ------------------------------------------------------------ */

  // If kappa is very small, return a circular uniform draw, as otherwise the
  // algorithm will fail.
  if (kp < .0000001) {
    return runif(n, 0, 2*pi);
  }

  NumericVector th(n);
  int sn;
  double a, b, r, u1, u2, u3, z, f, c;
  bool cont;

  a = 1 + sqrt(1 + 4.0 * pow(kp, 2));
  b = (a - (sqrt(2.0*a)))/(2.0*kp);
  r = (1 + pow(b,2))/(2.0*b);

  for (int i=0; i<n; i++) {

    cont = TRUE;

    do {
      u1 = runif(1, 0, 1)[0];
      u2 = runif(1, 0, 1)[0];
      u3 = runif(1, 0, 1)[0];

      // STEP 1
      z = cos(pi*u1);
      f = (1 + r*z)/(r + z);
      c = kp * (r - f);

      // STEP 2
      if (c*(2-c) - u2 > 0) cont=FALSE;

      // STEP 3
      if (log(c/u2) + 1 - c >= 0) cont=FALSE;
    } while (cont);

    // STEP 4
    if (u3 - 0.5 > 0) {
      sn = 1;
    } else {
      sn = -1;
    }

    th[i] = fmod(sn * acos(f) + mu, 2.0*pi);
  }

  return th;
}

// [[Rcpp::export]]
double Wapprox (double t) {
  // Approximation of Lamberts W.
  return exp(1)*t /
  (1 + 1 / (pow(2 * exp(1) * t + 2, -0.5) + 1 / (exp(1) - 1) - pow(2, -0.5)));
}

// [[Rcpp::export]]
vec sampleKappa(double etag, int eta) {
  // Function to sample values for kappa in a von Mises distribution. Etag should
  // be - R * cos(mu - theta_bar). eta is the posterior n, which is n + c where c is
  // the number of observations contained in the conjugate prior. For
  // uninformative, c = 0 and eta = n.

  // beta_0 in Forbes & Mardia (2014) is renamed g here to avoid confusion with
  // the intercept in the GLM model.
  double g, kl, ku, c1, c2, c3, c4, c5, c6, beta, eps,
  alph, x, u, v1, v2, v;

  long double k0, i0, r, kp_can;

  int nAttempts = 0;

  // Boolean indicating whether the current candidate is accepted.
  bool cont;

  g = etag / eta;

  // Setup: Compute approximately optimal parameters for the rejection
  // sampler. kappa_0 is called k0.
  kl   = 2 / (etag + sqrt(2 * eta + etag * etag));
  ku   = (2 + (1 / eta) ) / (etag + g + sqrt(2 * eta + 1 + etag * etag));
  c1   = 0.5 + 0.5 * (1 - 0.5 / eta) / eta;
  k0   = (1 - c1) * kl + c1 * ku;
  i0   = boost::math::cyl_bessel_i(0, k0);
  r    = boost::math::cyl_bessel_i(1, k0) / i0;
  c2   = 0.25 / eta - 2 / (3 * sqrt(eta));
  if (g < c2) {
    beta = g + 1;
  } else {
    beta = g + r + (1 - r) / (1 + 40 * eta * pow(g - c2, 2));
  }
  c3   = (log(i0) / k0 - beta + g) / (beta - g - r);
  c4   = Wapprox(c3 * exp(c3));
  eps  = c4 * k0 / (c3 - c4);
  alph = (beta - g - r) * (k0 + eps);
  c5   = log(i0);

  // Apply rejection sampler.
  cont = TRUE;
  do {

    // Draw values from the gamma distribution with the
    // tweaked parameters.
    x = rgamma(1, eta * alph + 1, 1.0 / (eta * beta))[0];

    // Save the number of candidates.
    nAttempts++;

    if (x > eps) {

      // Compute kp_can and v.
      kp_can = x - eps;
      c6 = 0.5 * log(2 * pi * kp_can) - kp_can;
      u  = runif(1, 0, 1)[0];
      v1 = log(u) / eta - (beta - g) * (kp_can - k0);
      v2 = alph * log((kp_can+eps) / (k0 + eps)) - c5;
      v  = v1 + v2;

      // Break the loop if these tests are passed.
      if (kp_can < 0.258 || v < c6) {
        if ( v < c6 - log(1 + 1 / 2 * kp_can) ||
             v < -log(boost::math::cyl_bessel_i(0, kp_can))) {
          cont = FALSE;
        }
      }

    }
  } while (cont);

  vec out = vec(2);
  out(0) = kp_can;
  out(1) = nAttempts;
  return out;
}



// [[Rcpp::export]]
double computeMeanDirection (vec th) {
  // Compute the mean direction for some dataset th.

  double C = as_scalar(arma::sum(cos(th)));
  double S = as_scalar(arma::sum(sin(th)));
  return atan2(S, C);
}


// [[Rcpp::export]]
double computeResultantLength (vec th) {
  // Compute the resultant length for some dataset th.

  double C = as_scalar(arma::sum(cos(th)));
  double S = as_scalar(arma::sum(sin(th)));
  return sqrt(pow(C, 2) + pow(S, 2));
}


// [[Rcpp::export]]
vec quantile(vec x, vec q) {
  // Compute a quantile of vector x for each proportion in vector q.

  int nx = x.size();
  int nq = q.size();

  std::sort(x.begin(), x.end());

  vec out = vec(nq);

  for (int i = 0; i < nq; i++)
  {
    out(i) = x[nx*(q(i) - 0.000000001)];
  }

  return out;
}

// [[Rcpp::export]]
vec circQuantile(arma::vec th, vec q) {
  // Compute a circular quantile.

  double rotation = computeMeanDirection(th) - pi;

  return quantile(th + rotation, q) - rotation;
}


// [[Rcpp::export]]
double estimateMode(vec x, double cip) {
  // Compute the mode using interval cip%.

  int n, cil, chiv;
  double len, M;

  n = x.size();

  std::sort(x.begin(), x.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the lower bound of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  len = x[cil]-x[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (len > (x[i+cil]-x[i])) {
      len = (x[i+cil]-x[i]);
      chiv = i;
    }
  }

  M = (x[chiv+cil]+x[chiv])/2;

  return M;
}



// [[Rcpp::export]]
vec computeHDI(vec x, double cip) {
  // Compute the cip% Highest Density Interval for some vector x.

  int n, cil, chiv;
  double len;

  n = x.size();

  std::sort(x.begin(), x.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  len = x[cil]-x[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (len > (x[i+cil]-x[i])) {
      len = (x[i+cil]-x[i]);
      chiv = i;
    }
  }

  vec M(2);
  M(0) = x[chiv];
  M(1) = x[chiv+cil];

  return M;
}


// [[Rcpp::export]]
double estimateDensity(vec x, double x_0, double cip) {
  // A method to estimate the density of a random variable x (here, often an
  // MCMC sample) at a value x_0. It uses a 'histogram' solution.

  int n, cil, loc, xidx_lb, xidx_ub;;

  n = x.size();


  // cout << endl << "n: " << n;


  sort(x.begin(), x.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Theoretical location of x_0.
  loc = sum(x < x_0);

  // cout << endl << "loc: " << loc;

  // Upper and lower bounds of x_0' part's theoretical histogram bin.
  if (loc - cil/2 < 0) {
    xidx_lb = 0;
    xidx_ub = cil;
  } else if (loc + cil/2 > n) {
    xidx_lb = n - cil - 1;
    xidx_ub = n - 1;
  } else {
    xidx_lb = loc - cil/2;
    xidx_ub = loc + cil/2;
  }

  // cout << endl << "LB, UB: " << xidx_lb << "," << xidx_ub << endl;

  double len = x(xidx_ub) - x(xidx_lb);

  double p = cip/len;

  return p;
}







// [[Rcpp::export]]
double rhsll(double b0, double kp, vec bt, vec dt,
             vec th, mat X, mat D, double r) {
  // The log-likelihood for use in the MCMC-algorithm. For speed-up the constant
  // addition can be skipped for MH-steps for which the left-hand side is indeed
  // a constant (that is, everything but kappa). Returns the right-hand side of
  // the log-likelihood function.

  //   cout << endl << "D dim: "  << D.n_rows  << ", " << D.n_cols  << endl;
  //   cout << endl << "dt dim: " << dt.n_rows << ", " << dt.n_cols << endl;
  //   cout << endl << "X dim: "  << X.n_rows  << ", " << X.n_cols  << endl;
  //   cout << endl << "bt dim: " << bt.n_rows << ", " << bt.n_cols << endl;
  //   cout << endl << "D*dt: " << endl << (D * dt) << endl;
  //   cout << endl << "X*bt: " << endl << (X * bt) << endl;
  //   cout << endl << "atanLF(X * bt, r): " << endl << atanLF(X * bt, r) << endl;
  //   cout << endl << "b0: " << endl << b0 << endl;
  //   cout << endl << "th: " << endl << th << endl;

  return kp * arma::sum(cos(th - b0 - (D * dt) - atanLF(X * bt, r)));
}




// [[Rcpp::export]]
double ll(double b0, double kp, vec bt, vec dt,
          vec th, mat X, mat D, double r) {
  // Compute the full log-likelihood.

  int n = th.size();

  // The left-hand side of the likelihood function.
  double lhs = - n * log(2 * pi) - n * log(boost::math::cyl_bessel_i(0, kp));

  // The right-hand side of the likelihood function.
  double rhs = rhsll(b0, kp, bt, dt, th, X, D, r);

  return lhs + rhs;
}


// [[Rcpp::export]]
vec logProbNormal (vec x, vec mu, vec sd)
{
  // Compute the log of the probability of the normal distribution for a vector
  // of x's, means and sd's.
  return -log(sd) - log (sqrt(2*pi)) -
    ( pow(x - mu, 2) / (2 * pow(sd, 2) ) );
}

template <typename T>
T normal_pdf(T x, T m, T s)
{
  static const T inv_sqrt_2pi = 0.3989422804014327;
  T a = (x - m) / s;

  return inv_sqrt_2pi / s * std::exp(-T(0.5) * a * a);
}


// bwb: Bandwith for the proposal for beta
// conj_prior: A vector of length 3, containing, in that order, prior mean
// direction, prior resultant length, and prior sample size. Used for beta_0 and
// kappa.
// bt_prior: A matrix of size K*2, containing prior means in the first column,
// and prior variances in the second column, for each respective predictor in
// the model.
// kappaModeEstBandwith: The bandwith to use in estimation of the mode for the
// residual concentration, kappa. Reasonable values are roughly between .005 and
// .2, although lower values may be reasonable if Q is large.
// CIsize: What % credible intervals should be returned?
// debug, loopDebug: debug is basic progress debug, loopDebug also prints where
// the sampler is in the first run of the main MCMC-loop.
//
// [[Rcpp::export]]
Rcpp::List circGLMC(vec th, mat X, mat D,
                    vec conj_prior, mat bt_prior,
                    vec starting_values, int burnin, int lag,
                    vec bwb, double kappaModeEstBandwith, double CIsize,
                    int Q, double r, bool returnPostSample,
                    int bt_prior_type, bool reparametrize,
                    bool debug, bool loopDebug,
                    bool groupMeanComparisons) {


  if (debug) std::cout << "--- Start --- " << std::endl << " - Initialize: a, ";

  //  Measure current state of clock for time taken.
  clock_t begin = clock();

  double C_psi, S_psi, R_psi, psi_bar, bt_lograt, dt_lograt, etag;
  double piOver2 = pi/2;

  // Sample size n and number of continuous predictors K and categorical
  // predictors J.
  int n = th.n_elem;

  // Initialize Beta_0, conjugate prior
  vec    b0_chain = zeros<vec>(Q);
  double b0_cur   = starting_values(0);
  double b0_prior = conj_prior(0);
  double R_prior  = conj_prior(1);
  int    n_prior  = conj_prior(2);
  double C_prior  = R_prior * cos(b0_prior);
  double S_prior  = R_prior * sin(b0_prior);
  int    n_post   = n + n_prior;

  // Data after subtraction of prediction.
  vec psi    = zeros<vec>(n);

  if (debug) std::cout << "b, ";

  // Initialize Kappa
  int nkpcan    = 0; // Number of candidates for kappa
  vec kp_chain  = zeros<vec>(Q);
  double kp_cur = starting_values(1);
  vec sk_res    = zeros<vec>(2); // Results of sampleKappa function

  if (debug) std::cout << "c, ";

  // Initialize Beta
  int K         = X.n_cols;
  // bool useBt    = K > 0;

  // if (useBt) {
  vec bt_cur    = zeros<vec>(K);
  vec zt_cur    = atanLF(bt_cur, 1/piOver2);
  vec nbtacc    = zeros<vec>(K); // Vector with # accepted beta's.
  mat bt_chain  = zeros<mat>(Q, K);
  mat zt_chain  = zeros<mat>(Q, K);

  // Sum of beta*x of each predictor, before link function.
  mat X_bybt = zeros<mat>(n, K);

  // Candidate for vector of predictors.
  vec bt_can = zeros<vec>(K);
  vec zt_can = zeros<vec>(K);
  // }

  double bt_cur_prior = 0;
  double bt_can_prior = 0;

  if (debug) std::cout << "d, ";

  // Initialize Delta
  int J         = D.n_cols;
  // bool useDt    = J > 0;

  vec dt_cur    = zeros<vec>(J);
  vec dt_can    = zeros<vec>(J);
  vec ndtacc    = zeros<vec>(J); // Vector with # accepted delta's.
  mat dt_chain  = zeros<mat>(Q, J);
  mat D_bydt    = zeros<mat>(n, J);

  if (debug) std::cout << "e, ";
  if (debug) std::cout << "f, ";

  // Initialize matrix with the log-likelihood of each data point, given the
  // current parameters.
  mat ll_each_th_curpars = zeros<mat>(Q, n);

  // Compute number of iterations, taking lag and burn-in into account.
  int Qbylag = Q * lag + burnin;
  int isav = 0; // Keeps track of where to save those values not thinned out.

  if (debug) std::cout << "g - End Initialize. -" << std::endl;
  if (debug) std::cout << "- J=" << J << ", K=" << K << "- Loop:" << std::endl;

  //  Measure current state of clock for time taken.
  clock_t init_done = clock();

  for (int i = 0; i < Qbylag; i++)
  {

    if (loopDebug & (i==0)) std::cout << "a, ";


    ////////////
    // BETA_0 //
    ////////////

    // Obtain an n*K matrix with beta_k * x_{i, k} in each cell.
    X_bybt  = X * bt_cur;
    // if (useBt) X_bybt.each_row() %= bt_cur.t();
    D_bydt  = D * dt_cur;
    // if (useDt) D_bydt.each_row() %= dt_cur.t();

    if (loopDebug & (i==0)) std::cout << "b, ";

    // Obtain psi and its current properties.
    psi     = th - arma::sum(D_bydt, 1) - atanLF(arma::sum(X_bybt, 1), r);
    C_psi   = as_scalar(sum(cos(psi))) + C_prior;
    S_psi   = as_scalar(sum(sin(psi))) + S_prior;
    R_psi   = sqrt(pow(C_psi, 2) + pow(S_psi, 2));
    psi_bar = atan2(S_psi, C_psi);

    if (loopDebug & (i==0)) std::cout << "c, ";

    // Obtain a new value for Beta_0
    b0_cur = rvmc(1, psi_bar, R_psi * kp_cur)[0];

    if (loopDebug & (i==0)) std::cout << "d, ";

    ///////////
    // KAPPA //
    ///////////

    // Sample a new value for kappa.
    etag    = - R_psi * cos(b0_cur - psi_bar);
    sk_res  = sampleKappa(etag, n_post);
    kp_cur  = sk_res(0);

    if (loopDebug & (i==0)) std::cout << "e, ";

    // After we are done with the burn-in, start gathering the amount of
    // nAttempts that we need every time.
    if (i >= burnin)
    {
      nkpcan += sk_res(1);
    }

    if (loopDebug & (i==0)) std::cout << "f, [[";

    ///////////
    // DELTA //
    ///////////

    for (int j = 0; j < J; j++) {

      if (loopDebug & (i==0)) std::cout << "dt_" << j+1 << "{1";

      // The proposal for each delta is a von Mises distribution with mean
      // dt_cur[j] with and as R_psi * kappa the residual kappa. This is fairly
      // arbitrary, but reasonable. kp_cur might be too wide (which will cause
      // low acceptance), 0 will probably work but the acceptance rate will be
      // very low. This is likely slightly too narrow, but good in most cases.
      dt_can(j) = rvmc(1, dt_cur[j], R_psi * kp_cur)[0];

      dt_lograt = rhsll(b0_cur, kp_cur, bt_cur, dt_can, th, X, D, r) -
        rhsll(b0_cur, kp_cur, bt_cur, dt_cur, th, X, D, r);

      if (loopDebug & (i==0)) std::cout << ", 3";

      // Accept the candidate according to the MH-ratio.
      if (dt_lograt > log(runif(1, 0, 1)[0]))
      {

        // Set delta if chosen, but ensure that the new delta is in [-pi, pi).
        dt_cur(j) = fmod(dt_can(j) + pi, 2.0*pi) - pi;

        if (i >= burnin) ndtacc(j)++;
      }
      else
      {
        dt_can(j) = dt_cur(j);
      }

      if (loopDebug & (i==0)) std::cout << ", 4}, ";

    }


    //////////
    // BETA //
    //////////

    // Draw and possibly accept candidates, in sequence, for each of the K
    // predictors.
    for(int k = 0; k < K; k++) {

      if (loopDebug & (i==0)) std::cout << "bt_" << k+1 << "{1";

      if (!reparametrize) {
        bt_can(k) += runif(1, -bwb(k), bwb(k))[0];
        if (loopDebug & (i==0)) std::cout << "2a";
      } else {
        zt_can(k) += runif(1, -bwb(k), bwb(k))[0];

        if (loopDebug & (i==0)) std::cout << "2b";

        // Because zeta = -1 and zeta = 1 are essentially equal when r=2, we
        // allow the algorithm to "jump" from one side to the other.
        if (r == 2.0) {
          if (zt_can(k) < -1) zt_can(k) += 2;
          if (zt_can(k) >= 1) zt_can(k) -= 2;
        }


        // Move the new zeta-candidate back to the beta parametrization, as the
        // likelihoods are assessed on that scale.
        bt_can(k) = invAtanLFdouble(zt_can(k), 1/piOver2);
      }

      if (loopDebug & (i==0)) std::cout << "3";

      // Set the priors for bt constant prior or normal prior.
      if (bt_prior_type == 1)
      {
        bt_can_prior = arma::sum(logProbNormal(bt_can,
                                               bt_prior.col(0),
                                               bt_prior.col(1)));
        bt_cur_prior = arma::sum(logProbNormal(bt_cur,
                                               bt_prior.col(0),
                                               bt_prior.col(1)));
      }

      if (loopDebug & (i==0)) std::cout << "4";

      bt_lograt = rhsll(b0_cur, kp_cur, bt_can, dt_cur, th, X, D, r) + bt_can_prior -
        rhsll(b0_cur, kp_cur, bt_cur, dt_cur, th, X, D, r) - bt_cur_prior;

      // Accept the candidate according to the MH-ratio.
      if (bt_lograt > log(runif(1, 0, 1)[0]))
      {
        bt_cur(k) = bt_can(k);

        if (i >= burnin) nbtacc(k)++;
      }
      else
      {
        bt_can(k) = bt_cur(k);
      }

      if (loopDebug & (i==0)) std::cout << "5}, ";
    }
    if (loopDebug & (i==0)) std::cout << "]], g, ";

    // Obtain reparametrized estimates.
    zt_cur = atanLF(bt_cur, 1/piOver2);

    if (loopDebug & (i==0)) std::cout << "h, ";

    // For non-thinned out iterations, save the current values.
    if ( (i % lag == 0) & (i >= burnin))
    {
      if (loopDebug & (i==0)) std::cout << "save, ";

      b0_chain(isav)     = b0_cur;
      kp_chain(isav)     = kp_cur;
      bt_chain.row(isav) = bt_cur.t();
      dt_chain.row(isav) = dt_cur.t();
      zt_chain.row(isav) = zt_cur.t();

      for (int lli = 0; lli < n; lli++) {

        ll_each_th_curpars(isav, lli) = ll(b0_cur, kp_cur, bt_cur, dt_cur,
                           (vec) th.subvec(lli, lli),
                           (mat) X.row(lli),
                           (mat) D.row(lli),
                           r);
      }

      isav++;
    }

    if (loopDebug & (i==0)) std::cout << "done. ";
  }

  if (debug) std::cout << "End Loop. -" << std::endl << "Gather results: a, ";

  //  Measure current state of clock for time taken.
  clock_t loop_done = clock();

  //////////////////////
  // GATHER ESTIMATES //
  //////////////////////

  // Bounds for the CCI's
  vec qbounds = vec(2);
  qbounds(0)  = (1-CIsize)/2.0;
  qbounds(1)  = 1-((1-CIsize)/2.0);

  ////
  // Estimates Beta_0
  ////
  double b0_meandir = computeMeanDirection(b0_chain);
  vec b0_CCI  = circQuantile(b0_chain, qbounds);

  if (debug) std::cout << "b, ";

  ////
  // Estimates Kappa
  ////
  double kp_mean    = arma::mean(kp_chain);
  double kp_mode    = estimateMode(kp_chain, kappaModeEstBandwith);
  vec kp_HDI  = computeHDI(kp_chain, CIsize);

  if (debug) std::cout << "c, ";


  ////
  // Estimates Delta
  ////
  rowvec dt_mean    = zeros<rowvec>(J);
  dt_mean           = arma::mean(dt_chain, 0);

  // Get a mean direction for delta.
  rowvec dt_meandir = zeros<rowvec>(J);
  mat dt_CCI  = mat(2, J);

  // Obtain the mean direction and
  for (int j = 0; j < J; j++)
  {
    dt_meandir(j) = computeMeanDirection(dt_chain.col(j));

    dt_CCI.col(j) = circQuantile(dt_chain.col(j), qbounds);
  }


  if (debug) std::cout << "d, ";




  ////
  // Estimates Beta/Zeta
  ////

  rowvec bt_mean    = zeros<rowvec>(K);
  bt_mean           = arma::mean(bt_chain, 0);

  rowvec zt_mean    = zeros<rowvec>(K);
  zt_mean           = arma::mean(zt_chain, 0);

  if (debug) std::cout << "e, ";

  // Get a mean direction representation of the predictors (considering the
  // periodicity of their space of origin), which may or may not perform better
  // than the simple linear mean.
  rowvec zt_meandir = zeros<rowvec>(K);
  for (int predi = 0; predi < K; predi++)
  {
    zt_meandir(predi) = computeMeanDirection(zt_chain.col(predi)*pi)/pi;
  }

  if (debug) std::cout << "f, ";

  // Matrix with CCI's for beta and zeta.
  mat zt_CCI  = mat(2, K);
  mat bt_CCI  = mat(2, K);

  if (debug) std::cout << "g, ";

  for (int predi = 0; predi < K; predi++)
  {
    bt_CCI.col(predi) = quantile(bt_chain.col(predi), qbounds);
    zt_CCI.col(predi) = quantile(zt_chain.col(predi), qbounds);
  }

  if (debug) std::cout << "h, ";


  // Proportion of values for kappa that was accepted.
  double propacckp    = (double) Q * (double) lag / (double) nkpcan;
  vec propaccbt = nbtacc / ((double) Q * (double) lag);
  vec propaccdt = ndtacc / ((double) Q * (double) lag);

  if (debug) std::cout << "i, ";

  // Obtain predicted values for all data points.
  vec th_hat = (D * dt_meandir.t()) + atanLF(X * bt_mean.t(), r) + b0_meandir;

  if (debug) std::cout << "j, ";





  //// INFORMATION CRITERIA

  if (debug) std::cout << "(IC 1, ";

  // Log-likelihood of the data given the posterior estimates.
  double ll_th_estpars = ll(b0_meandir, kp_mode, bt_mean.t(), dt_meandir.t(),
                            th, X, D, r);

  // Log-likelihood of the data for each current posterior value in the chain.
  vec ll_th_curpars = sum(ll_each_th_curpars, 1);

  if (debug) std::cout << "IC 2, ";

  double lppd = accu(log(sum(exp(ll_each_th_curpars), 0))) - n * log(Q);

  if (debug) std::cout << "IC 3, ";

  // Calculate the amount of parameters.
  int n_par = 2 + J + K;

  // Obtain the Bayes-AIC, where the Bayesian estimates of the parameters are
  // used instead (in ll_data).
  double AIC_Bayes = - 2 * ll_th_estpars + 2 * n_par;

  // Obtain the two versions of DIC as in Gelman's BDA, 3rd ed.
  double p_DIC     = 2 * (ll_th_estpars - mean(ll_th_curpars));
  double p_DIC_alt = 2 * var(ll_th_curpars);

  double DIC     = - 2 * ll_th_estpars + 2 * p_DIC;
  double DIC_alt = - 2 * ll_th_estpars + 2 * p_DIC_alt;

  if (debug) std::cout << "IC 4, ";

  // Obtain the two versions of WAIC as in Gelman's BDA, 3rd ed.
  rowvec WAIC_logofmean = log(mean(exp(ll_each_th_curpars), 0));
  rowvec WAIC_meanoflog = mean(ll_each_th_curpars, 0);
  double p_WAIC1     = 2 * sum(WAIC_logofmean - WAIC_meanoflog);

  double p_WAIC2     = accu(var(ll_each_th_curpars, 0, 0));

  double WAIC1      = - 2 * (lppd - p_WAIC1);
  double WAIC2      = - 2 * (lppd - p_WAIC2);

  if (debug) std::cout << "IC 5)";





  if (debug) std::cout << endl << endl << "- BFs: - " << endl;

  //// BAYES FACTORS

  // Interval chosen for posterior density estimation.
  double cip = 1000.0/Q;

  // Prevent crashes for small Q runs. The BFs will be unreliable in this case.
  if (Q < 1000) cip = .1;

  if (debug) std::cout << "(Delta: ";

  //// Automatic Bayes Factors for delta's.
  vec dtIneqBFs = zeros<vec>(J);
  vec prop_pos_dts = zeros<vec>(J);

  for (int j = 0; j < J; j++) {
    prop_pos_dts(j) = ((double) sum(dt_chain.col(j) > 0)) / Q;

    if (debug) std::cout << "(" << prop_pos_dts(j) << ", ";

    // This is the fit/complexity Bayes Factor for delta_k > 0 vs. delta_k < 0.
    dtIneqBFs(j) = prop_pos_dts(j)/(1 - prop_pos_dts(j));

    if (debug) std::cout <<  dtIneqBFs(j) << "), ";
  }




  if (debug) std::cout << ") " << endl << "(Beta: ";



  //// Automatic Bayes Factors for beta's
  vec btIneqBFs    = zeros<vec>(K);
  vec btSDDBFs     = zeros<vec>(K);
  vec btSDDprior   = zeros<vec>(K);
  vec prop_pos_bts = zeros<vec>(K);

  for (int k = 0; k < K; k++) {
    prop_pos_bts(k) = ((double) sum(bt_chain.col(k) > 0)) / Q;

    if (debug) std::cout << "(" << prop_pos_bts(k) << ", ";

    // This is the fit/complexity Bayes Factor for beta_k > 0 vs. beta_k < 0.
    btIneqBFs(k) = prop_pos_bts(k)/(1 - prop_pos_bts(k));

    if (debug) std::cout <<  btIneqBFs(k) << ", ";

    // Null vs. Alternative Savage-Dickey Density BFs.
    btSDDprior(k) = normal_pdf((double) 0.0,
               (double) bt_prior(k, 0),
               (double) bt_prior(k, 1));


    btSDDBFs(k) = estimateDensity(bt_chain.col(k), 0, cip) / btSDDprior(k);

    if (debug) std::cout <<  btSDDBFs(k) << "), ";
  }




  if (debug) std::cout << ") " << endl << "(Mu_comp: ";







  //// Automatic Bayes Factors for mean comparisons
  int nmu = J + 1;

  // Number of comparisons
  int n_comp = (nmu*nmu - nmu)/2;

  vec muIneqBFs    = zeros<vec>(n_comp);
  vec muSDDBFs = zeros<vec>(n_comp);
  vec prop_lgr = zeros<vec>(n_comp);
  vec comp_chain = zeros<vec>(Q);

  // Create chains containing the means of each group in each iteration.
  mat mu_chain = zeros<mat>(Q, nmu);
  mu_chain.col(0) = b0_chain;


  for (int mui = 1; mui < nmu; mui++) {
    mu_chain.col(mui) = b0_chain + dt_chain.col(mui-1);
  }

  if (groupMeanComparisons & (J > 0)) {


    int comp_i = 0;

    for (int mi_first = 0; mi_first < nmu - 1; mi_first++) {
      for (int mi_last = mi_first + 1; mi_last < nmu; mi_last++) {


        if (debug) {cout<<endl<<"Comp:"<<mi_first<<","<<mi_last;}



        // Inequality Bayes Factors for mu_a > mu_b vs. mu_a < mu_b.
        prop_lgr(comp_i) = ((double) sum(mu_chain.col(mi_first) >
                                           mu_chain.col(mi_last))) / Q;

        muIneqBFs(comp_i) = prop_lgr(comp_i)/(1 - prop_lgr(comp_i));

        // Null vs. Alternative Savage-Dickey Density BFs.

        // Note that the the prior is chosen by default to be uniform on the
        // circle. This is defensible and even desired from an empirical Bayes
        // standpoint, and very easy in the circular case. It is questionable
        // from a subjective Bayes standpoint: After observing on group to lie
        // near a specific point, we generally do have some knowledge where the
        // other groups will lie.
        comp_chain = mu_chain.col(mi_first) - mu_chain.col(mi_last);

        muSDDBFs(comp_i) = estimateDensity(comp_chain, 0, cip)*2*pi;

        comp_i++;
      }
    }
  }



  if (debug) std::cout << ")" << endl << "Output list. ";






  //// OUTPUT LIST

  // Create a list of outputs.
  Rcpp::List out;

  out["b0_meandir"] = b0_meandir;
  out["b0_CCI"]     = b0_CCI;

  out["kp_mean"]    = kp_mean;
  out["kp_mode"]    = kp_mode;
  out["kp_HDI"]     = kp_HDI;
  out["kp_propacc"] = Rcpp::wrap(propacckp);

  out["bt_mean"]    = bt_mean;
  out["bt_CCI"]     = bt_CCI;
  out["bt_propacc"] = Rcpp::wrap(propaccbt.t());

  out["dt_meandir"] = dt_meandir;
  out["dt_CCI"]     = dt_CCI;
  out["dt_propacc"] = Rcpp::wrap(propaccdt.t());

  out["zt_mean"]    = zt_mean;
  out["zt_mdir"]    = zt_meandir;
  out["zt_CCI"]     = zt_CCI;


  out["lppd"]       = lppd;
  out["n_par"]      = n_par;

  out["ll_th_estpars"] = ll_th_estpars;

  if (returnPostSample) {
    out["ll_each_th_curpars"] = ll_each_th_curpars;
    out["ll_th_curpars"]      = ll_th_curpars;
    out["th_hat"]             = th_hat;
  }

  out["AIC_Bayes"]        = AIC_Bayes;

  out["p_DIC"]      = p_DIC;
  out["p_DIC_alt"]  = p_DIC_alt;
  out["DIC"]        = DIC;
  out["DIC_alt"]    = DIC_alt;

  out["p_WAIC1"]    = p_WAIC1;
  out["p_WAIC2"]    = p_WAIC2;
  out["WAIC1"]      = WAIC1;
  out["WAIC2"]      = WAIC2;

  // out["prop_pos_dts"] = prop_pos_dts;
  out["DeltaIneqBayesFactors"] = dtIneqBFs;
  out["BetaIneqBayesFactors"]  = btIneqBFs;
  out["BetaSDDBayesFactors"]  = btSDDBFs;
  if (groupMeanComparisons) {
    out["MuIneqBayesFactors"] = muIneqBFs;
    out["MuSDDBayesFactors"] = muSDDBFs;
  }

  out["SavedIts"]   = Rcpp::wrap(Q);
  out["TotalIts"]   = Rcpp::wrap(Qbylag);


  if (debug) std::cout << ", h. ";

  if (returnPostSample)
  {
    out["b0_chain"] = b0_chain;
    out["kp_chain"] = kp_chain;
    out["bt_chain"] = bt_chain;
    out["dt_chain"] = dt_chain;
    out["zt_chain"] = zt_chain;

    out["ll_each_th_curpars"] = ll_each_th_curpars;

    if (groupMeanComparisons) {out["mu_chain"] = mu_chain;}
  }

  if (debug) std::cout << "End Gather -" << std::endl << std::endl;


  // Save the time taken.
  clock_t end = clock();

  vec time_taken = zeros<vec>(4);
  time_taken(0) = double(init_done - begin) / CLOCKS_PER_SEC;
  time_taken(1) = double(loop_done - init_done) / CLOCKS_PER_SEC;
  time_taken(2) = double(end - loop_done) / CLOCKS_PER_SEC;
  time_taken(3) = double(end - begin) / CLOCKS_PER_SEC;

  out["TimeTaken"] = time_taken;

  return out;
}
