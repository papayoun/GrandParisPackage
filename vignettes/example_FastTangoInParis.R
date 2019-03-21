rm(list = ls())
library(GrandParisPackage)
library(parallel)
set.seed(122)
trueTheta <- pi/4; trueSigma2 <- 1;
n <- 1000; times <- seq(0, by = 1, length = n);
simulated_POD <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
observations <- simulated_POD[, "observations"]
estTheta = T; estSig = F; updateOrders <- rep(T, n)
firstStep <- 0.5
nStart <- floor(n / 3)
steps <- c(rep(firstStep, nStart), firstStep * (1:(n - nStart) )^(-0.500001))
gradientSteps   <- matrix(steps , ncol = 2, nrow = n)
allRes <- mclapply(1:detectCores(), function(seed){
  set.seed(seed)
  thetaStart <- runif(1, 0, 2 * pi)
  Res <- fastTangOR(observations, times,
                    thetaModel = thetaStart, sigma2 = trueSigma2,
                    updateOrders = updateOrders, gradientSteps = gradientSteps, 
                    all = T, estimateSigma2 = F, randomWalkParam = 1)
}, mc.cores = detectCores())

thetaEst <- sapply(allRes, function(x) x$Estimates[,1]) %% (2*pi)
smoothEst <- function(vecEst, n0 = 500){
  n <- length(vecEst)
  if(n0 >= n)
    n0 <- max(1, floor(n / 10))
  output <- vecEst
  tailEst <- vecEst[(n0 + 1) : n]
  output[(n0 + 1) : n] <- cumsum(tailEst) * 1 / ((n0 + 1) : n - n0)
  output
}
matplot(thetaEst, type = "b", lty = 1, pch = 20)
thetaMoy <- apply(thetaEst, 2, smoothEst)
matplot(thetaMoy, type = "b", lty = 1, pch = 20)
abline(h = c(trueTheta + (-5:5) * 2 * pi, trueTheta + pi), lwd = 2, 
       col = "orange")
