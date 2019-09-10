
suppressMessages(library(TMB))
suppressMessages(library(tidyquant))
suppressMessages(library(doParallel))
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
#simData(N = 10, P = 2)
TMBPars <- function(No.Cols = 4,No.Rows = 100){
  return(list(logit_Rho = rep(0.95,No.Cols), log_sigma_proc = rep(1,No.Cols),
              log_sigma_obs = rep(1,No.Cols),states = matrix(0, No.Rows+1,No.Cols)))
}

local1SSM <- dyn.load(dynlib("local1SSM"))
local2SSM <- dyn.load(dynlib("local3DSSM"))
local4SSM <- dyn.load(dynlib("local4SSM"))
Optimisation <- function(N = 1000, P = 4, sdObs, sdPro, Rho, init, DLL = "local4SSM"){
  x <- simData(N = N, P = P, Rho = Rho, sigmaPro = sdPro, sigmaObs = sdObs, init = init); #print("Shumani Rocks")
  obj1 <- MakeADFun(list(y = x$Y), TMBPars(No.Cols = P,No.Rows = N), random="states",DLL=DLL)
  # Minimise the nll
  #print("Shumani Rocks")
  suppressMessages(opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr))
  
  est <- tryCatch(summary(sdreport(obj1)), error=function(e) NA)
  return(tryCatch(est[rownames(est)%in%c("sigmaObsA", "sigmaObsB","sigmaObsC", "sigmaObsD",
                                    "sigmaProA", "sigmaProB","sigmaProC", "sigmaProD",
                                    "RhoA","RhoB","RhoC","RhoD"),1],
           error=function(e) NA))
}

##create matrix to store summary results from each run
simResEst <- matrix(0,6,12); simResSd <- matrix(0,6,12) #matrix to store estimations and their SD
colnames(simResEst) <- c("sigmaObsA", "sigmaObsB","sigmaObsC", "sigmaObsD",
                      "sigmaProA", "sigmaProB","sigmaProC", "sigmaProD",
                      "RhoA","RhoB","RhoC","RhoD")
colnames(simResSd) <- c("sigmaObsA", "sigmaObsB","sigmaObsC", "sigmaObsD",
                      "sigmaProA", "sigmaProB","sigmaProC", "sigmaProD",
                      "RhoA","RhoB","RhoC","RhoD")

sdObs <- c(0.175,0.2, 0.15, 0.22)
Rho <- c(0.87,0.99, 0.98,0.89)
sdPro <- matrix(data = rbind(c(0.165,0.19, 0.14, 0.20),c(0.175,0.2, 0.15, 0.22),
                             c(0.185,0.21, 0.16, 0.23),5*c(0.175,0.2, 0.15, 0.22),
                             8*c(0.175,0.2, 0.15, 0.22),10*c(0.175,0.2, 0.15, 0.22))
                , 6, 4)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

results <- foreach (j = 1:6, .combine = "c") %dopar%{
  library(TMB)
  local4SSM <- dyn.load(dynlib("local4SSM"))
  estimates <- matrix(0, 100, 12)
  colnames(estimates) <- paste(c("sigmaObsA", "sigmaObsB","sigmaObsC", "sigmaObsD",
                           "sigmaProA", "sigmaProB","sigmaProC", "sigmaProD",
                           "RhoA","RhoB","RhoC","RhoD"),j)
  set.seed(19950311)
  for (i in 1:100){
    estimates[i,] <- Optimisation(N = 1000, P = 4, sdObs = sdObs,
                                  sdPro = sdPro[j,], Rho = Rho,
                                  init = c(1,1,1,1), DLL = "local4SSM")}
  #estimates
  #write.csv(estimates, file = "estimates.csv")
  #tryCatch(est[rownames(est)%in%"states",1],error = function(e) NA)
  #obj1$env$parList()$states
  
  #for (i in 1:ncol(estimates)){
  #  hist(estimates[,i],
  #       xlab = colnames(estimates)[i],
  #       main = paste0("histogram of ",colnames(estimates[i])))
  #}
  #simResEst[j,] <- apply(estimates,2,mean)
  #simResSd[j,] <- apply(estimates,2,sd)
  list(est = apply(estimates, 2, mean), se = apply(estimates,2,sd))
}
#write.csv(simResEst, "paramsEst with seed.csv")
#write.csv(simResSd, "paramsSd with seed.csv")
stopCluster(cl)
results