#include <TMB.hpp>

// Local 3D State-space model
 //you have to add the multivariate normal negative likelihood and you are done with this template
//also try and figure out which value to use.

//using namespace Rcpp;
using namespace density;

template<class Type>
Type objective_function<Type>::operator() ()
{
  // data:
  DATA_MATRIX(y);
  
  
  // parameters:
  PARAMETER_VECTOR(logit_Rho); // autocorrelation, logit of Rho
  PARAMETER_VECTOR(log_sigma_proc); // log(process SD)
  PARAMETER_VECTOR(log_sigma_obs); // log(observation SD)
  PARAMETER_MATRIX(states); // unobserved state vector
  
  // procedures: (transformed parameters)
  Type sigmaProA = exp(log_sigma_proc(0));
  Type sigmaObsA = exp(log_sigma_obs(0));
  Type RhoA = 2.0/(1+exp(-logit_Rho(0)))-1;
  
  Type sigmaProB = exp(log_sigma_proc(1));
  Type sigmaObsB = exp(log_sigma_obs(1));
  Type RhoB = 2.0/(1+exp(-logit_Rho(1)))-1;
  
  Type sigmaProC = exp(log_sigma_proc(2));
  Type sigmaObsC = exp(log_sigma_obs(2));
  Type RhoC = 2.0/(1+exp(-logit_Rho(2)))-1;
  
  Type sigmaProD = exp(log_sigma_proc(3));
  Type sigmaObsD = exp(log_sigma_obs(3));
  Type RhoD = 2.0/(1+exp(-logit_Rho(3)))-1;
  
  // reports on transformed parameters:
  ADREPORT(sigmaProA);
  ADREPORT(sigmaProB);
  ADREPORT(sigmaProC);
  ADREPORT(sigmaProD);
  ADREPORT(sigmaObsA);
  ADREPORT(sigmaObsB);
  ADREPORT(sigmaObsC);
  ADREPORT(sigmaObsD);
  ADREPORT(RhoA);
  ADREPORT(RhoB);
  ADREPORT(RhoC);
  ADREPORT(RhoD);
  
  // create and initialize covariance matrix for proc and obs
  matrix<Type> proCov(4,4);
    proCov << sigmaProA*sigmaProA, 0, 0, 0,
            0, sigmaProB*sigmaProB, 0, 0,
            0, 0, sigmaProC*sigmaProC, 0,
            0, 0, 0, sigmaProD*sigmaProD;
  
  matrix<Type> obsCov(4,4);
    obsCov << sigmaObsA*sigmaObsA, 0, 0, 0,
              0, sigmaObsB*sigmaObsB, 0, 0,
              0, 0, sigmaObsC*sigmaObsC, 0,
              0, 0, 0, sigmaObsD*sigmaObsD;
  
  matrix<Type> InitCov(4, 4);
    InitCov << 10 * 10, 0, 0, 0,
               0, 10 * 10, 0, 0,
               0, 0, 10 * 10, 0,
               0, 0, 0, 10 * 10;
  
  //create and initialize MVN object for each process and the initial values.
  MVNORM_t<Type> proDensity(proCov);
  MVNORM_t<Type> obsDensity(obsCov);
  MVNORM_t<Type> initDensity(InitCov);
  
  Type nll = 0.0; // initialize negative log likelihood
  vector<Type> mu(4);
  mu(0) = 0.0;
  mu(1) = 0.0;
  mu(2) = 0.0;
  mu(3) = 0.0;
  
  nll += initDensity(mu); //Initial negative likelihood
  
  // process model and negative likelihood:
  for (int i = 1; i < states.rows(); ++i) {
    mu(0) = states(i, 0) - RhoA * states(i - 1, 0);
    mu(1) = states(i, 1) - RhoB * states(i - 1, 1);
    mu(2) = states(i, 2) - RhoC * states(i - 1, 2);
    mu(3) = states(i, 3) - RhoD * states(i - 1, 3);
    nll += proDensity(mu); //process ll
  }
  //observation model and negative likelihood:
  for(int i = 0; i < y.rows(); ++i){
    mu = y.row(i) - states.row(i+1);
    nll += obsDensity(mu); // observation ll
  }
  //REPORT(states);
  return nll;
 }
 
