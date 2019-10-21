#Adapted from Auger....

# Simulation parameters
# equivalent of rho =1
sdObs <- 0.2
sdProSeq <- c(0.02,0.04,0.1,0.2,0.4,1,2)
sd0 <- 0.01
n <- 200

# TMB parameters and functions
library(TMB)
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

#############################
# Make a loop to look at how the distribution of the parameter estimates
nSim <- 200

# Matrix that will keep track of information
parEstTMB <- matrix(nrow=nSim*length(sdProSeq), ncol=3)
colnames(parEstTMB) <- c("sdObsEst", "sdProEst", "sdProSim")

rmseEstTMB <- matrix(nrow=nSim*length(sdProSeq), ncol=2)
colnames(rmseEstTMB) <- c("estPar", "simPar")

set.seed(19961012)

for(k in seq_along(sdProSeq)){
  for(j in 1:nSim){
    indexSim <- (k-1)*nSim+j
    parEstTMB[indexSim,3] <- sdProSeq[k]
    
    # Simulate the data
    dat <- simData(N = n, Rho = 0.25, sigmaPro = sdProSeq[k], sigmaObs = sdObs,init = 0)
    
    # fit using TMB
    # Put the data into the TMB format
    dataTMB <- list(y = dat$Y)
    
    # Fix x1 to 0
    xmap <- 1:length(dat$states)
    xmap[1] <- NA
    # Making the objective function for our TMB model
    # Note that the true unobserved locations are considered random variables
    obj1 <- MakeADFun(dataTMB, TMBPars(No.Rows = n, Y = dat$Y), 
                      map=list(logit_Rho=factor(NA),states=factor(xmap)), random="states",DLL="local1SSM")
    # Minimise the nll
    opt1 <- nlminb(obj1$par,obj1$fn)
    
    # Parameter estimates
    srep <- tryCatch(summary(sdreport(obj1)), error=function(e) NA)
    parEstTMB[indexSim,1:2] <- tryCatch(srep[rownames(srep)%in%c("sigma_obs", "sigma_proc"),1],
                                        error=function(e) rep(NA,2))
    
    # RMSE
    stsParEst <- obj1$env$parList()$states
    rmseEstTMB[indexSim,1] <- sqrt(sum((stsParEst - dat$states)[-1]^2)/n)
    # With sim values
    val <- c(log(sdObs),log(sdProSeq[k]))
    obj1$fn(val)
    stsSimVal <- obj1$env$parList(val)$states
    rmseEstTMB[indexSim,2] <- sqrt(sum((stsSimVal - dat$states)[-1]^2)/n)
  }
}

quartz(width=6, height=6)
layout(matrix(1:(3*length(sdProSeq)),nrow=length(sdProSeq),byrow=TRUE))
par(mar=c(1.2,1.4,0.3,0.5), mgp=c(0.8,0.3,0), tck=-0.03, las=1, oma=c(1.5,1.6,0,0.3))
sum(is.na(parEstTMB[,2]))
yleg <- 0.1
xleg <- -0.7
placl <- "topleft"
lt <- c(LETTERS, "AA", "BB")
for(i in 1:length(sdProSeq)){
  sdProIndex <- parEstTMB[,3] == sdProSeq[i]
  hh <- hist(na.omit(parEstTMB[sdProIndex,1]),
             breaks=seq(0,max(parEstTMB[,1],na.rm=TRUE), length.out=50), plot=FALSE)
  plot(hh, ylim =c(0,round(max(hh$counts)*1.24)),
       xlab="", main="", ylab="",
       border=FALSE, col="darkgrey")
  abline(v=sdObs)
  ys <- seq(0.982,0.04,-0.157)
  title(ylab="Frequency", line=0.4, cex.lab=1.2, outer=TRUE, 
        adj=ys[i])
  rat <- sdObs/sdProSeq[i]
  legend("top",
         legend = substitute(sigma[epsilon] == rat*~sigma[eta],list(rat=rat)), bty="n")
  box()
  legend(placl, lt[i*3-2], bty="n", x.intersp=xleg, y.intersp=yleg)
  
  hist(parEstTMB[sdProIndex,2],
       breaks=seq(0,max(parEstTMB[sdProIndex,2],na.rm=TRUE),length.out=50), xlab="", main="",
       ylab="",
       border=FALSE, col="darkgrey")
  abline(v=unique(parEstTMB[sdProIndex,3]), col="red")
  box()
  legend(placl, lt[i*3-1], bty="n", x.intersp=xleg, y.intersp=yleg)
  
  # RMSE
  hEst <- hist(na.omit(rmseEstTMB[sdProIndex,1]),
               breaks=seq(min(rmseEstTMB[sdProIndex,1:2], na.rm=TRUE),
                          max(rmseEstTMB[sdProIndex,1:2], na.rm=TRUE), length.out=50), 
               plot=FALSE)
  hTrue <- hist(na.omit(rmseEstTMB[sdProIndex,2]),
                breaks=seq(min(rmseEstTMB[sdProIndex,1:2], na.rm=TRUE),
                           max(rmseEstTMB[sdProIndex,1:2], na.rm=TRUE), length.out=50),
                plot=FALSE)
  plot(hEst, 
       ylim=c(0,round(max(hEst$counts,hTrue$counts)*1.21)),
       xlab="", main="", ylab="", border=FALSE, col=grey(0.5))
  plot(hTrue, border=FALSE, col=rgb(0,0,1,0.5), add=TRUE)
  box()
  legend(placl, lt[i*4], bty="n", x.intersp=xleg, y.intersp=yleg)
}
title(xlab=substitute(widehat(sigma)[epsilon]), line=0.3, cex.lab=1.2,outer=TRUE, 
      adj=0.18)
title(xlab=substitute(widehat(sigma)[eta]), line=0.3, cex.lab=1.2,outer=TRUE, 
      adj=0.51)
title(xlab="RMSE", line=0.1, cex.lab=1.2,outer=TRUE, 
      adj=0.87)
