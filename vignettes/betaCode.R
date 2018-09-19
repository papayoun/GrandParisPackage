rm(list = ls())
library(GrandParisPackage)
library(parallel)
set.seed(122)
trueTheta <- pi / 4; trueSigma2 <- 1;
n <- 2000; times <- seq(0, by = 2, length = n);
SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
observations <- SINEprocess[, "observations"]
estTheta = T; estSig = F; updateOrders <- rep(T, n)
firstStep <- 0.5
nStart <- floor(n / 6)
steps <- c(rep(firstStep, nStart), firstStep * (1:(n - nStart) )^(-0.5000001))
gradientSteps   <- matrix(steps , ncol = 2, nrow = n)
allRes <- mclapply(1:12, function(seed){
  set.seed(seed)
  thetaStart <- runif(1, 0, 2 * pi)
  Res <- fastTangOR(observations, times,
                    thetaModel = thetaStart, sigma2 = trueSigma2,
                    updateOrders = updateOrders, gradientSteps = gradientSteps, 
                    all = T, estimateSigma2 = F, randomWalkParam = 1)
}, mc.cores = detectCores())

thetaEst <- sapply(allRes, function(x) x$Estimates[,1]) %% (2*pi)
smoothEst <- function(vecEst, n0 = 100){
  n <- length(vecEst)
  if(n0 >= n)
    n0 <- max(1, floor(n / 10))
  output <- vecEst
  tailEst <- vecEst[(n0 + 1) : n]
  output[(n0 + 1) : n] <- cumsum(tailEst) * 1 / ((n0 + 1) : n - n0 + 1)
}
matplot(thetaEst, type = "b", lty = 1, pch = 20)
thetaMoy <- apply(thetaEst, 2, smoothEst)
matplot(thetaMoy, type = "b", lty = 1, pch = 20)
abline(h = c(trueTheta + (-5:5) * 2 * pi, trueTheta + pi), lwd = 2, col = "orange")


# With unknown sigma ------------------------------------------------------

Res <- fastTangOR(observations, times,
                  thetaModel = thetaStart, sigma2 = trueSigma2,
                  updateOrders = updateOrders, gradientSteps = gradientSteps, 
                  all = T, estimateSigma2 = T, randomWalkParam = 1)

testedThetas <- seq(0, 2*pi, length.out = 100)
processValues <- SINEprocess[, "hiddenStates"]
logDensities <- sapply(testedThetas, function(theta){
  SINE_density(theta = theta, X = processValues, times = times, MC.size = 60, log = T)
})
gradLogDensities <- sapply(testedThetas, function(theta){
  mean(SINE_logDensityGradient(theta = theta, X = processValues, times = times, MC.size = 30))
})
partsRes <- sapply(allRes, function(x){
  colSums(x$Particles * x$Weights)
})
matplot(t(part1), type = "p", pch = 20, cex = 0.5, col = 1)
plot(SINEprocess[, 1], col = "blue")
lines(partsRes[,1], col = "green")
lines(partsRes[,4], col = "purple")

matplot(t(part2), type = "p", pch = 20, cex = 0.5, col = 2, add = T)


matplot(t(Res$Particles), pch = 20, type = "p", col = 1)
apply(Res$Particles, 2, var)
matplot(t(Res$Weights), pch = 20, type = "p", col = 1)

check <- function(i){
  parts <- Res$Particles[,i]
  observationVariance <- Res$Estimates[i - 1, 2]
  obs <- observations[i]
  myGrad <- function(particles, observation, observationVariance = 1){
    gaussTerm = dnorm(particles, observation, sqrt(observationVariance)) 
    quadrTerm = (particles - observation)^2;
    output =  gaussTerm * (quadrTerm / observationVariance - 1) / (2 * observationVariance)
    output
  }
  myLogGrad <- function(particles, observation, observationVariance = 1){
    gaussTerm = dnorm(particles, observation, sqrt(observationVariance)) 
    quadrTerm = (particles - observation)^2;
    output =  (quadrTerm / observationVariance - 1) / (2 * observationVariance)
    output
  }
  c(mean(dnorm(parts, obs, sqrt(observationVariance))),mean(myGrad(parts, obs, observationVariance)),
    mean(myLogGrad(parts, obs, observationVariance)))
}
