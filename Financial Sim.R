
library(TMB)
library(tidyquant)

Y <- cbind(FANG$adjusted[FANG$symbol=="FB"],
           FANG$adjusted[FANG$symbol=="NFLX"],
           FANG$adjusted[FANG$symbol=="AMZN"],
           FANG$adjusted[FANG$symbol=="GOOG"])
rownames(Y) <- FANG$date[FANG$symbol=="FB"]
colnames(Y) <- c("FB","NFLX","AMZN","GOOG")
Y <- Y[-nrow(Y),]/Y[-1,]

simData <- function(N = 1000,P = 4,
                    Rho = c(0.75,0.75,0.75,0.75),
                    sigmaPro = c(0.1,0.1,0.1,0.1),
                    sigmaObs = c(0.2,0.2,0.2,0.2),
                    init = c(0,0,0,0)){
  states <- matrix(data = init[1:P],ncol = P,nrow = N+1,byrow = TRUE)
  Rho <- Rho[1:P]
  sigmaPro <- sigmaPro[1:P]
  sigmaObs <- sigmaObs[1:P]
  errPro <- matrix(data = 0, nrow = N, ncol = P)
  for (i in 1:P){errPro[,i] <- rnorm(N,0,sigmaPro[i])}
  
  for (i in 1:N+1){states[i,] = Rho*states[i-1,]+errPro[i-1,]}
    
  errObs <- matrix(data = 0, nrow = N, ncol = P)
  for (i in 1:P){
    errObs[,i] <- rnorm(N,0,sigmaObs[i])
  }
  Y <- states[-1,]+errObs
  return(list(Y=Y, states = states, Rho = Rho, 
              sigmaPro = sigmaPro, sigmaObs = sigmaObs))
}
simData(N = 10, P = 2)
TMBPars <- function(No.Cols = 4,No.Rows = 100){
  return(list(logit_Rho = rep(0.95,No.Cols), log_sigma_proc = rep(1,No.Cols),
              log_sigma_obs = rep(1,No.Cols),states = matrix(0, No.Rows+1,No.Cols)))
}

local1SSM <- dyn.load(dynlib("local1SSM"))
local3SSM <- dyn.load(dynlib("local3DSSM"))
local4SSM <- dyn.load(dynlib("local4SSM"))
Optimisation <- function(N = 1000, P = 4, sdObs, sdPro, Rho, init, DLL = "local4SSM"){
  x <- simData(N = N, P = P, Rho = Rho, sigmaPro = sdPro, sigmaObs = sdObs, init = init); #print("Shumani Rocks")
  obj1 <- MakeADFun(list(y = x$Y), TMBPars(No.Cols = P,No.Rows = N), random="states",DLL=DLL)
  # Minimise the nll
  #print("Shumani Rocks")
  opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr)
  
  est <- tryCatch(summary(sdreport(obj1)), error=function(e) NA)
  return(tryCatch(est[rownames(est)%in%c("sigmaObsA", "sigmaObsB","sigmaObsC", "sigmaObsD",
                                    "sigmaProA", "sigmaProB","sigmaProC", "sigmaProD",
                                    "RhoA","RhoB","RhoC","RhoD"),1],
           error=function(e) NA))
}

estimates <- matrix(0,300,12)
colnames(estimates) <- c("sigmaObsA", "sigmaObsB","sigmaObsC", "sigmaObsD",
                         "sigmaProA", "sigmaProB","sigmaProC", "sigmaProD",
                         "RhoA","RhoB","RhoC","RhoD")
for (i in 1:300){
  estimates[i,] <- Optimisation(N = 500, P = 4, sdObs = c(0.35, 0.4, 0.3,0.44),
                                sdPro = c(0.175,0.2, 0.15, 0.22), Rho = c(0.87,0.99, 0.98,0.89),
                                init = c(1,1,1,1), DLL = "local4SSM")}
estimates
write.csv(estimates, file = "estimates.csv")
#tryCatch(est[rownames(est)%in%"states",1],error = function(e) NA)
#obj1$env$parList()$states
