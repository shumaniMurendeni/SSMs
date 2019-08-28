#state space modelling
setwd("Documents/Can State Space Models Be Problematic?")
#Load TMB library and compile the model
####################################### 
library(TMB)
library(tidyquant)
compile("local1SSM.cpp",flags = "-Wno-unused-variable", package = "TMB")
dyn.load(dynlib("local1SSM"))

#functions to help with simulation
#######################################
#setup parameters and import data.

FANG %>% group_by(symbol)
Y <- as.matrix(FANG[FANG$symbol=="GOOG",6]); #acf(Y) The data's mean changes with time.
Y <-diff(Y); #acf(Y) The data seems to have constant mean after differencing.

logit <- function(x){return(-log(2/(x+1) - 1))}
pars <- function(log_sigma_proc = 0.1, log_sigma_obs = 0.2,rho =0.5 , Y = Y){
  params <- list(logit_Rho = log(rho), log_sigma_proc = log_sigma_proc, log_sigma_obs =log_sigma_obs, u = rep(mean(Y),length(Y)+1))
  return(params)
}
Rho <- function(x) return(2/(1+exp(-x))-1)

params <- pars(Y=Y)

obj <- MakeADFun(data = list(y = Y), parameters = params,map=list(u=factor(params$u)), 
                 random = "u", DLL = "local1SSM")

opt1 <- optim(obj$par,obj$fn,obj$gr)
opt1
