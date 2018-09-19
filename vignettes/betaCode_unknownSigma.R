rm(list = ls())
library(GrandParisPackage)
library(parallel)
set.seed(122)
trueTheta <- pi / 4; trueSigma2 <- 1;
n <- 50; times <- seq(0, by = 2, length = n);
SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
observations <- SINEprocess[, "observations"]
estTheta = T; estSig = F; updateOrders <- rep(T, n)
firstStep <- 0.5
nStart <- floor(n / 6)
stepsTheta <- c(rep(firstStep, nStart), firstStep * (1:(n - nStart))^(-0.5000001))
stepsSigma <- 0.1 * (1:n)^(-0.50001)
gradientSteps   <- cbind(stepsTheta, stepsSigma)
allRes <- mclapply(1:12, function(seed){
  set.seed(seed)
  thetaStart <- trueTheta
  Res <- fastTangOR(observations, times,
                    thetaModel = thetaStart, sigma2 = trueSigma2, particleSize = 5,
                    updateOrders = updateOrders, gradientSteps = gradientSteps, 
                    all = T, estimateSigma2 = T, randomWalkParam = 1)
}, mc.cores = detectCores())

thetaEst <- sapply(allRes, function(x) x$Estimates[,1]) %% (2*pi)
sigmaEst <- sapply(allRes, function(x) x$Estimates[,2])
par(mfrow = c(1, 2))
matplot(thetaEst, type = "b", lty = 1, pch = 20, main = expression(theta), ylim = c(0,6))
matplot(sigmaEst, type = "b", lty = 1, pch = 20, main = expression(sigma))
