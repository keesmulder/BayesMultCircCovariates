
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//using namespace arma;

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <iostream>
#include <math.h>

using namespace Rcpp;


// TODO:
// Check if quantiles are obtained properly.

arma::vec atanlf(arma::vec x, double c) {
  // arctangent-family link function, spreading range c, where c=1 corresponds to
  // the semicircle.
  return c * atan(x);
}


NumericVector rvmc(int n, double mu, double kp) {
  /* FUNCTION rvmc -------------------------------------------
  Generate random variates from the von Mises distribution.

  n:      The number of random variates required.
  mu:     The required mean direction, mu.
  kp:     The required concentration, kappa.

  Returns: A vector of length n containing VM random variates.
  ------------------------------------------------------------ */

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
      z = cos(4*atan(1)*u1);
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

    th[i] = fmod(sn * acos(f) + mu, 8.0*atan(1));
  }

  return th;
}


double Wapprox (double t) {
  // Approximation of Lamberts W.
  return exp(1)*t / (1 + 1 / (pow(2 * exp(1) * t + 2, -0.5) + 1 / (exp(1) - 1) - pow(2, -0.5)));
}


arma::vec sampleKappa(double etag, int eta) {
  // beta_0 in Forbes & Mardia (2014) is renamed g here to avoid confusion with
  // the intercept in the GLM model.
  double g, kl, ku, c1, c2, c3, c4, c5, c6, k0, i0, r, beta, eps,
         alph, x, kp_can, u, v1, v2, v;

  int nAttempts = 0;
  // Indicating whether the candidate is accepted.
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

  // Apply rejection sampler
  cont = TRUE;
  do {

    // Draw values from the gamma distribution with the
    // tweaked parameters
    x = rgamma(1, eta * alph + 1, 1.0 / (eta * beta))[0];

    // Save the number of candidates
    nAttempts++;

    if (x > eps) {

      // Compute kp_can and v
      kp_can = x - eps;
      c6 = 0.5 * log(8 * atan(1) * kp_can) - kp_can;
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

  arma::vec out;
  out(0) = kp_can;
  out(0) = nAttempts;
  return out;
}



double computeMeanDirection (arma::vec th) {
// Compute the mean direction for some dataset th.
  double C = arma::as_scalar(arma::sum(cos(th)));
  double S = arma::as_scalar(arma::sum(sin(th)));
  return atan2(S, C);
}


double computeResultantLength (arma::vec th) {
// Compute the resultant length for some dataset th.
  double C = arma::as_scalar(arma::sum(cos(th)));
  double S = arma::as_scalar(arma::sum(sin(th)));
  return sqrt(pow(C, 2) + pow(S, 2));
}


// [[Rcpp::export]]
arma::vec quantile(arma::vec x, arma::vec q) {
  int nx = x.size();
  int nq = q.size();

  std::sort(x.begin(), x.end());

  arma::vec out = arma::vec(nq);

  for (int i = 0; i < nq; i++)
  {
    out(i) = x[nx*(q(i) - 0.000000001)];
  }

  return out;
}

// [[Rcpp::export]]
arma::vec circularQuantile(arma::vec th, arma::vec q) {

  double rotation = computeMeanDirection(th) - 4 * atan(1);

  return quantile(th + rotation, q) - rotation;
}


// [[Rcpp::export]]
double estimateMode(arma::vec x, double cip) {

  int n, cil, chiv;
  double len, M;

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

  M = (x[chiv+cil]+x[chiv])/2;

  return M;
}



// [[Rcpp::export]]
arma::vec computeHDI(arma::vec x, double cip) {

  int n, cil, chiv;
  double len;

  n = x.size();

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

  arma::vec M(2);
  M(0) = x[chiv];
  M(1) = x[chiv+cil];

  return M;
}


// [[Rcpp::export]]
double rhsll(double b0, double kp, arma::vec bt,
             arma::vec th, arma::mat X, double c) {
  // The log-likelihood for use in the MCMC-algorithm. For speed-up the constant
  // addition can be skipped for MH-steps for which the left-hand side is indeed a
  // constant (ie., not kappa)
  // The right-hand side of the log-likelihood function.

  arma::mat X_bybt = X;
  X_bybt.each_row() %= bt.t();

  double rhs = arma::sum(cos(th - b0 - atanlf(arma::sum(X_bybt, 1), c)));

  return rhs * kp;
}


// [[Rcpp::export]]
double ll(double b0, double kp, arma::vec bt,
          arma::vec th, arma::mat X, double c) {
// Full log-likelihood

  int n = th.size();

  // The left-hand side of the likelihood function.
  double lhs = - n * log(boost::math::cyl_bessel_i(0, kp));

  // The right-hand side of the likelihood function.
  double rhs = rhsll(b0, kp, bt, th, X, c);

  return lhs + rhs;
}


// [[Rcpp::export]]
arma::vec logProbNormal (arma::vec x, arma::vec mu, arma::vec sd)
{
  return -log(sd) - log (sqrt(8.0 * atan(1))) -
          ( pow(x - mu, 2) / (2 * pow(sd, 2) ) );
}

// [[Rcpp::export]]
double computeLogBtPost (double b0, double kp,
                arma::vec bt, arma::vec th, arma::mat X, double c,
                arma::vec mu_prior, arma::vec sd_prior) {

  double loglik = ll(b0, kp, bt, th, X, c);

  double logprior = 0;

  logprior = arma::sum(logProbNormal(bt, mu_prior, sd_prior));

  return loglik + logprior;
}



//# bwb: Bandwith for the proposal for beta
//circGLM <- function(th, X, linkfun, invlinkfun,
//                    b0_start=0, kp_start=1, bt_start=rep(0, ncol(X)),
//                    bwb=rep(.5, ncol(X)), Q=1000) {

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
// [[Rcpp::export]]
Rcpp::List circGLM(arma::vec th, arma::mat X,
                   arma::vec conj_prior, arma::mat bt_prior,
                   arma::vec starting_values, int burnin, int lag,
                   arma::vec bwb, double kappaModeEstBandwith, double CIsize,
                   int Q, double c, bool returnPostSample, char bt_prior_type) {


  double C_psi, S_psi, R_psi, psi_bar, bt_lograt, etag;
  int nkpcan = 0; // Number of candidates for kappa

  // Sample size n and number of predictors K.
  int n = th.n_elem;
  int K = X.n_cols;

  arma::vec b0_chain, kp_chain = arma::vec(Q);
  arma::mat bt_chain = arma::mat(Q, K);

  double b0_cur    = starting_values(0);
  double kp_cur    = starting_values(1);
  arma::vec bt_cur = starting_values(arma::span(2, K + 2));

  double b0_prior = conj_prior(0);
  double R_prior  = conj_prior(1);
  int    n_prior  = conj_prior(2);
  double C_prior  = R_prior * cos(b0_prior);
  double S_prior  = R_prior * sin(b0_prior);
  int    n_post   = n + n_prior;

  b0_chain(0)     = b0_cur;
  kp_chain(0)     = kp_cur;
  bt_chain.col(0) = bt_cur;

  arma::vec sk_res(2); // To hold the results of the sampleKappa-function
  arma::vec X_bybt(n); // Sum of beta*x of each predictor, before link function.
  arma::vec psi(n);    // Data after subtraction of prediction.
  arma::vec bt_can(K); // Candidate for vector of predictors.

  // Compute number of iterations, taking lag and burn-in into account.
  int Qbylag = Q * lag + burnin;
  int isav = 0; // Keeps track of where to save those values not thinned out.

  for (int i = 1; i < Qbylag; i++)
  {
    // Obtain an n*K matrix with beta_k * x_{i, k} in each cell.
    X_bybt  = X;
    X_bybt.each_row() %= bt_cur;

    // Obtain psi and its current properties.
    psi     = th - atanlf(sum(X_bybt, 1), c);
    C_psi   = arma::as_scalar(arma::sum(cos(psi))) + C_prior;
    S_psi   = arma::as_scalar(arma::sum(sin(psi))) + S_prior;
    R_psi   = sqrt(pow(C_psi, 2) + pow(S_psi, 2));
    psi_bar = atan2(S_psi, C_psi);

    b0_cur = rvmc(1, psi_bar, R_psi * kp_cur)[0];

    bt_can = bt_cur;

    for(int k = 0; k < K; k++) {

      bt_can(k) += runif(-1, -bwb(k), bwb(k))[0];

      if (bt_prior_type == "constant") {
        bt_lograt = rhsll(b0_cur, kp_cur, bt_can, th, X, c) -
                    rhsll(b0_cur, kp_cur, bt_cur, th, X, c);
      }
      else if (bt_prior_type == "normal")
      {
        bt_lograt = rhsll(b0_cur, kp_cur, bt_can, th, X, c) +
                    logProbNormal(bt_can, bt_prior.col(0), bt_prior.col(1)) -
                    rhsll(b0_cur, kp_cur, bt_cur, th, X, c) -
                    logProbNormal(bt_cur, bt_prior.col(0), bt_prior.col(1));
      }

      if (bt_lograt > log(runif(1, 0, 1)[0]))
      {
        bt_cur(k) = bt_can(k);
      }
      else
      {
        bt_can(k) = bt_cur(k);
      }
    }

    etag = - R_psi * cos(b0_cur - psi_bar);

    sk_res  = sampleKappa(etag, n_post);
    kp_cur  = sk_res(0);

    // After we are done with the burn-in, start gathering the amount of nAttempts
    // that we need every time.
    if (i >= burnin)
    {
      nkpcan += sk_res(1);
    }

    // For non-thinned out iterations, save the current values.
    if ( (i % lag == 0) & (i >= burnin))
    {
      b0_chain(isav)     = b0_cur;
      kp_chain(isav)     = kp_cur;
      bt_chain.row(isav) = bt_cur;

      isav++;
    }
  }



  //  Estimated values
  double b0_meandir = computeMeanDirection(b0_chain);
  double kp_mean    = arma::mean(kp_chain);
  double kp_mode    = estimateMode(kp_chain, kappaModeEstBandwith);
  arma::vec bt_mean = arma::mean(bt_chain, 0);


  // Bounds for the CCI's, which
  arma::vec qbounds = arma::vec(2);
  qbounds(0) = (1-CIsize)/2.0;
  qbounds(1) = 1-((1-CIsize)/2.0);


  arma::vec b0_CCI  = circularQuantile(b0_chain, qbounds);
  arma::vec kp_HDI  = computeHDI(kp_chain, CIsize);

//  matrix with CCI's for beta.
  arma::mat bt_CCI  = arma::mat(K, 2);

  for (int predi = 0; predi < K; predi++)
  {
     bt_CCI.row(predi) = quantile(bt_chain.col(predi), qbounds);
  }

  // Proportion of values for kappa that was accepted.
  float propacc = Q * lag / nkpcan;

  Rcpp::List out;

  out["b0_meandir"] = b0_meandir;
  out["kp_mean"]    = kp_mean;
  out["kp_mode"]    = kp_mode;
  out["b0_CCI"]     = b0_CCI;
  out["kp_HDI"]     = kp_HDI;
  out["bt_CCI"]     = bt_CCI;
  out["PropAcc"]    = Rcpp::wrap(propacc);



  if (returnPostSample)
  {
////  Try _
    out["b0_chain"]   = b0_chain;
    out["kp_chain"]   = kp_chain;
    out["bt_chain"]   = bt_chain;
  }


  return out;
}
