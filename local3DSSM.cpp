// Local 3D State-space model
#include <TMB.hpp> //you have to add the multivariate normal negative likelihood and you are done with this template
//also try and figure out which value to use.

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

	// reports on transformed parameters:
	ADREPORT(sigmaProA);
	ADREPORT(sigmaProB);
	ADREPORT(sigmaObsA);
	ADREPORT(sigmaObsB);
	ADREPORT(RhoA);
	ADREPORT(RhoB);
	
	// create and initialize covariance matrix for proc and obs
	matrix<Type> proCov(2,2);
		proCov << sigmaProA*sigmaProA, 0,
			0, sigmaProB*sigmaProB;

	matrix<Type> obsCov(2,2);
		obsCov << sigmaObsA*sigmaObsA, 0,
			0, sigmaObsB*sigmaObsB;

	matrix<Type> InitCov(2, 2);
		InitCov << 10 * 10, 0, 0, 10 * 10;

	//create and initialize MVN object for each process and the initial values.
	MVNORM_t<Type> proDensity(proCov);
	MVNORM_t<Type> obsDensity(obsCov);
	MVNORM_t<Type> initDensity(InitCov);

	Type nll = 0.0; // initialize negative log likelihood
	vector<Type> mu(2);
	mu(0) = 0.0;
	mu(1) = 0.0;

	nll += initDensity(mu); //Initial negative likelihood
	
	// process and observation model and negative likelihood:
	for (int i = 1; i < states.rows(); ++i) {
		mu(0) = states(i, 0) - RhoA * states(i - 1, 0);
		mu(1) = states(i, 1) - RhoB * states(i - 1, 1);
		nll += proDensity(mu); //process ll
	}
	for(int i = 0; i < y.rows(); ++i){
		//mu(0) = y(i,0)-states(i+1,0);
		//mu(1) = y(i,1)-states(i+1,1);
		mu = y.row(i) - states.row(i);
		nll += obsDensity(mu); // observation ll
	}
	REPORT(states);
	return nll;/**/
	//return 0.0;
}
