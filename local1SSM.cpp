//1D local SSM
#include <TMB.hpp> //https://kaskr.github.io/adcomp/_book/Toolbox.html#time-series

template<class Type>
Type objective_function<Type>::operator() ()
  {
  // data:
  DATA_VECTOR(y);
  
  // parameters:
  PARAMETER(logit_Rho); // density dependence parameter
  PARAMETER(log_sigma_proc); // log(process SD)
  PARAMETER(log_sigma_obs); // log(observation SD)
  PARAMETER_VECTOR(states); // unobserved state vector
  
  // procedures: (transformed parameters)
  Type Rho = 2.0/(1+exp(-logit_Rho))-1;
  Type sigma_proc = exp(log_sigma_proc);
  Type sigma_obs = exp(log_sigma_obs);
  
  // reports on transformed parameters:
  ADREPORT(Rho)
  ADREPORT(sigma_proc)
  ADREPORT(sigma_obs)
  ADREPORT(states)
    
  int n = y.size(); // get time series length
  
  Type nll = 0.0; // initialize negative log likelihood
  
  for(int i = 1; i < n; i++){
    // process model:
    Type m = Rho * states[i - 1];
    nll -= dnorm(states[i], m, sigma_proc, true);
    
    //Observation model
    nll -= dnorm(y[i-1], states[i], sigma_obs, true);
  }
  REPORT(states);
  return nll;
}
