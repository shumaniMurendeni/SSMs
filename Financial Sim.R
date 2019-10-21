
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

#Load TMB Template
dyn.load(dynlib("local4SSM"))

#Optimisation Function
Optimisation <- function(N = 1000, P = 4, sdObs, sdPro, Rho, init,
			Y = 1, DLL = "local4SSM"){
  x <- simData(N = N, P = P, Rho = Rho, sigmaPro = sdPro, sigmaObs = sdObs, init = init); #print("Shumani Rocks")
  if (length(Y)<2){Y=x$Y}
  obj1 <- MakeADFun(list(y = Y), TMBPars(No.Cols = P,No.Rows = N), random="states",DLL=DLL)
  # Minimise the nll
  #print("Shumani Rocks")
  #suppressMessages(
  opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr)#)
  
  est <- tryCatch(summary(sdreport(obj1)), error=function(e) NA)
  return(tryCatch(est[rownames(est)%in%c("sigmaObsA", "sigmaObsB","sigmaObsC", "sigmaObsD",
                                    "sigmaProA", "sigmaProB","sigmaProC", "sigmaProD",
                                    "RhoA","RhoB","RhoC","RhoD"),1],
           error=function(e) NA))
}

#set seed
set.seed(20190616)

##get parameter estimates using real data
parsEst <- Optimisation(N=1007,4,c(1,1,1,1),c(1,1,1,1),c(0.5,0.5,0.5,0.5),
			c(1,1,1,1), Y = Y)

sdPro <- parsEst[5:8]
Rho <- parsEst[9:12]
sdObs <- parsEst[1:4]
parsEst

cl <- makeCluster(detectCores())
registerDoParallel(cl)

results <- foreach (j = 1:200, .combine = "rbind") %dopar%{
  library(TMB)
  dyn.load(dynlib("local4SSM"))
  Optimisation(N = 1007, P = 4, sdObs = sdObs,
               sdPro = sdPro, Rho = Rho,
               init = c(1,1,1,1), DLL = "local4SSM")
}
stopCluster(cl)

results
save.image("1.RData")
#par(mfrow = c(3,4))

for (i in 1:nrow(results)){
	hist(results[,i], freq = F, xlab = names(results)[i])}
