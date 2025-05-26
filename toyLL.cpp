//  Likelihood Computation for toySWIFT
//  (Version 3.0, January 3, 2023)
//  (c) Ralf Engbert & Maximilian M. Rabe, Universität Potsdam
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <Rcpp.h>
#include <iostream>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
// [[Rcpp::NumericVector]]
double toyLL_cpp(NumericVector lfreq, IntegerVector fpos,NumericVector fdur,double nu,double r,double mt,double iota,double beta) {

  double  kappa = 0.0;
  double  eta = -2.0;
  double  gamma = 1.0;

  double  sigma = 1.0/(1.0+2.0*nu+nu*nu);     // normalization of processing span
  int     NW = lfreq.length();          // number of words
  int     Nfix = fpos.length();         // number of fixations
  NumericVector amax = 1.0-beta*lfreq;   // word frequency dependent maximum activation
  NumericVector a(NW);        // word activation
  NumericVector s(NW);        // word saliency
  NumericVector p(NW);        // selection probability
  NumericVector lambda(NW);   // processing rate
  double  shape = 9.0;        // shape parameter of the gamma distribution
  double  factorial = 8.0*7.0*6.0*5.0*4.0*3.0*2.0;
  double  rate = shape/mt;    // rate parameter of the gamma distribution
  double  rate2 = 0.0;        // modulated rate
  double  dgamma = 0.0;       // gamma density value
  double  leftact = 0.0;      // activation left of fixation
  double  loglik = 0.0;       // total likelihood
  double  LLtime = 0.0;       // temporal likelihood
  double  LLspat = 0.0;       // spatial likelihood
  double  tfix = 0.0;         // current fixation duration
  int     k = 1;              // currently fixated word

  // simulation loop
  tfix = fdur[0];    // fixation duration
  k = fpos[0];       // fixated word
  for ( int j=1; j<=(Nfix-1); j++ ) {

    // 1. Update processing rates
    lambda.fill(0.0); // reset all lambda values
    if ( k-1>=1 ) lambda[k-1 -1] = nu*sigma;
    lambda[k -1] = sigma;
    if ( k+1<=NW ) lambda[k+1 -1] = nu*sigma;
    if ( k+2<=NW ) lambda[k+2 -1] = nu*nu*sigma;

    // 2. Evolve activations
    a = a + r*lambda*tfix/1000.0;
    for ( int l=1; l<=NW; l++ )
      if ( a[l -1]>amax[l -1] )  a[l -1] = amax[l -1];

      // Compute word saliencies
      s = amax*sin(M_PI*a/amax) + pow(10.0,eta);

      // Compute probability for target selection
      p = pow(s,gamma)/sum(pow(s,gamma));

      // 3. Spatial loglik
      k = fpos[j+1 -1];   // fixation position observed in the data
      LLspat = LLspat + log(p[k -1]);

      // 4. Temporal loglik
      tfix = fdur[j+1 -1];   // fixation duration observed in the data
      leftact = 1.0;
      if ( k>1 )  for ( int l=1; l<=k-1; l++ )  leftact *= 1.0+kappa*(s[l -1] - pow(10.0,eta));
      rate2 = rate*(1.0+iota*(a[k -1]))/leftact;
      dgamma = shape*log(rate2) - log(factorial) + (shape-1.0)*log(tfix) - tfix*rate2;
      LLtime = LLtime + dgamma;
  }
  loglik = LLtime + LLspat;
  return loglik;
}


