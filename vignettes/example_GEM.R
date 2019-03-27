# library(GrandParisPackage)
library(parallel)
seed <- 122
set.seed(seed)
trueTheta <- pi / 4; trueSigma2 <- 1;
n <- 100; times <- seq(0, to = 50, length = n);
SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
observations <- SINEprocess[, "observations"]
n_case <- 12
allRes <- mclapply(1:n_case, function(i){
  thetaStart <- runif(1, 0, 2 * pi)
  sigma2Start <- runif(1, 0.2, 2)
  myTry <- GEM(observations = observations, observationTimes = times, thetaStart = thetaStart, 
               sigma2Start = sigma2Start, nIterations = 4, nModels = 20)
  colnames(myTry) = c("theta", "sigma2")
  myTry
}, mc.cores = min(n_case, detectCores()))

thetaEst <- sapply(allRes, function(x) x[,"theta"]) %% (2 * pi)
sigmaEst <- sapply(allRes, function(x) x[, "sigma2"])
par(mfrow = c(2, 1))
matplot(thetaEst, col = 1, type = "b", pch = 3, lty = 1, main = expression(theta))
abline(h = trueTheta, col = "red", lty = 2, lwd = 3)
matplot(sigmaEst, col = 1, type = "b", pch = 3, lty = 1, main = expression(sigma^2))
abline(h = trueSigma2, col = "red", lty = 2, lwd = 3)
