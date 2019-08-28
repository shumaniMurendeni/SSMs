
library(TMB)
library(tidyquant)

Y <- cbind(FANG$adjusted[FANG$symbol=="FB"],
           FANG$adjusted[FANG$symbol=="NFLX"],
           FANG$adjusted[FANG$symbol=="AMZN"],
           FANG$adjusted[FANG$symbol=="GOOG"])
rownames(Y) <- FANG$date[FANG$symbol=="FB"]
colnames(Y) <- c("FB","NFLX","AMZN","GOOG")

simData <- function(N = 1000,
                    Rho = c(0.75,0.75,0.75,0.75),
                    sigmaPro = c(0.1,0.1,0.1,0.1),
                    sigmaObs = c(0.2,0.2,0.2,0.2),
                    init = c(0,0,0,0)){
  states <- matrix(data = init,ncol = 4,nrow = N+1,byrow = TRUE)
  for (i in 2:(N+1)){
    err <- c(rnorm(1,0,sigmaPro[1]),rnorm(1,0,sigmaPro[2]),
             rnorm(1,0,sigmaPro[3]),rnorm(1,0,sigmaPro[4]))
    states[i,] = Rho*states[i-1,]+err
  }
  err <- c(rnorm(N,0,sigmaObs[1]),rnorm(N,0,sigmaObs[2]),
           rnorm(N,0,sigmaObs[3]),rnorm(N,0,sigmaObs[4]))
  Y <- states[-1,]+err
  return(list(Y=Y, states = states, Rho = Rho, 
              sigmaPro = sigmaPro, sigmaObs = sigmaObs))
}

TMBPars <- function(No.Cols = 4,No.Rows = 100){
  return(list(logit_Rho = rep(0.95,No.Cols), log_sigma_proc = rep(1,No.Cols),
              log_sigma_obs = rep(1,No.Cols),states = matrix(0, No.Rows+1,No.Cols)))
}

local1SSM <- dyn.load(dynlib("local1SSM"))
local3SSM <- dyn.load(dynlib("local3DSSM"))
local4SSM <- dyn.load(dynlib("local4SSM"))

x <- simData(N = 50)
obj1 <- MakeADFun(list(y = Y), TMBPars(No.Cols = 4,No.Rows = nrow(Y)), random="states",DLL="local4SSM")
# Minimise the nll
opt1 <- optim(obj1$par,obj1$fn,obj1$gr)

tryCatch(summary(sdreport(obj1)), error=function(e) NA)

#print("Wavhudi Gloria Matshevha")