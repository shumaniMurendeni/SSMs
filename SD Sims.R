
#Load R packages
suppressMessages(library(TMB))
#suppressMessages(library(tidyquant))
suppressMessages(library(doParallel))
#suppressMessages(library(devtools))

#find_rtools()

## Load the TMB Template
#compile("local1SSM.cpp")
dyn.load(dynlib("local1SSM"))

## Function to simulate qualitative data
simData <- function(N = 1000,
                    Rho = 0.75,
                    sigmaPro = 0.1,
                    sigmaObs = 0.2,
                    init = 1){
  states <- matrix(data = init,ncol = 1,nrow = N+1,byrow = TRUE)
  Rho <- Rho
  sigmaPro <- sigmaPro
  sigmaObs <- sigmaObs
  errPro <- matrix(data = 0, nrow = N, ncol = 1)
  
  for (i in 1:N+1){states[i,] = Rho*states[i-1,]+rnorm(1,0,sigmaPro)}
  
  errObs <- rnorm(N,0,sigmaObs)
  Y <- states[-1,]+errObs
  return(list(Y=Y, states = states, Rho = Rho, 
              sigmaPro = sigmaPro, sigmaObs = sigmaObs))
}

## Create initial parameters
TMBPars <- function(No.Rows = 100, Y = 1){
  return(list(logit_Rho = 0.95, log_sigma_proc = 1,
              log_sigma_obs = 1,states = matrix(mean(Y), No.Rows+1,1)))
}

## Optim function
Optimisation <- function(N = 1000, sdObs, sdPro, Rho, init, DLL = "local1SSM"){
  x <- simData(N = N, Rho = Rho, sigmaPro = sdPro, sigmaObs = sdObs, init = init);
  TMBpars <- TMBPars(No.Rows = N, x$Y)
  xmap <- 1:length(TMBpars$states)
  xmap[1] <- NA
  
  obj <- MakeADFun(list(y = x$Y),TMBpars, map=list(states=factor(xmap)),
                    random="states",DLL=DLL)
  # Minimise the nll
  #print("Shumani Rocks")
  suppressMessages(opt1 <- nlminb(obj$par,obj$fn,obj$gr))
  
  est <- tryCatch(summary(sdreport(obj)), error=function(e) NA)
  params = tryCatch(est[rownames(est)%in%c("sigma_proc", "sigma_obs","Rho"),1],
                    error=function(e) NA)
  rmse <- sqrt(sum((obj$env$parList()$states-x$states)^2)/N)
  return(list(pars = params, err = rmse))
}

cl <- makeCluster(detectCores())
registerDoParallel(cl)

sdObs <- 0.2
sdPro <- c(0.02,0.04,0.1,0.2,0.4,1,2)
Rho <- 0.50
par(mfrow = c(4,4))
results <- foreach (j = sdPro) %dopar%{
  library(TMB)
  dyn.load(dynlib("local1SSM"))
  res <- matrix(0, ncol = 4, nrow = 200)
  colnames(res) <- c("Rho", "SD Proc","SD Obs", "RMSE")
  for (i in 1:200){
    re <- Optimisation(N = 1000, sdObs = sdObs, sdPro = j, Rho = Rho, init = 0, DLL = "local1SSM")
    res[i,1:3] <- re$pars
    res[i,4] <- re$err
  }
  
  hist(res[,2], xlab = "sd Proc", freq = F, main = "Process sd")
  abline(v = j, col = "blue")
  
  hist(res[,3], xlab = "sd Obs", freq = F, main = "Observation sd")
  abline(v = sdObs, col = "blue")
  
  hist(res[,1], xlab = expression(rho), freq = F, main = expression(rho))
  abline(v = Rho, col = "blue")
  
  hist(res[,4], xlab = "RMSE", freq = F, main = "RMSE")
  return(res)
}

stopCluster(cl)
save.image("sdSim.RData")
