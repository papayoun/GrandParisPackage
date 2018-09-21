rm(list = ls())
library(GrandParisPackage)
library(parallel)
set.seed(1)
trueTheta <- pi / 4; trueSigma2 <- 0.2;
n <- 10; times <- seq(0, by = 2, length = n);
estTheta = T; estSig = F; updateOrders <- rep(T, n)
firstStep <- 0.5
nStart <- floor(n / 6)
stepsTheta <- c(rep(firstStep, nStart), firstStep * (1:(n - nStart))^(-0.5000001))
stepsSigma <- 1 * (1:n)^(-0.6)
gradientSteps   <- cbind(stepsTheta, stepsSigma)
SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
observations <- SINEprocess[, "observations"]
set.seed(1)
# indFin <- min(n, 19)
indFin <- n
RW <- 1
# getWeights <- function(indexDeb, weights = F){
#   obs <- observations[indexDeb + 1]
#   oldParticles <- Res$Particles[, indexDeb]
#   newParticles <- Res$Particles[, indexDeb + 1]
#   theta <- Res$Estimates[indexDeb, 1]; sigma2 <- Res$Estimates[indexDeb, 2]
#   delta = diff(times)[1]
#   quantities <- Prop_evalTrans(oldParticles, newParticles, obs, delta, RW, theta, sigma2)
#   if(weights)
#     return(quantities[, "trans"] * quantities[, "obs"] / quantities[, "prop"])
#   return(quantities)
# }
Res <- fastTangOR(observations[1:indFin], times[1:indFin], backwardSampleSize = 2,
                  thetaModel = pi, sigma2 = trueSigma2, particleSize = 250,
                  updateOrders = updateOrders, gradientSteps = gradientSteps, 
                  all = T, estimateSigma2 = T, estimateTheta = T, randomWalkParam = RW)
# oneRes <- Res
# matplot(t(oneRes$Weights),   type = "b", pch = 20, col = 1)
# matplot(t(oneRes$Particles), type = "b", pch = 20, col = 1)
# lines(observations, col = "blue")
# lines(SINEprocess[, "hiddenStates"], col = "red")
# thetaEst <- Res$Estimates[,1] %% (2*pi)
# sigmaEst <- Res$Estimates[,2]
# par(mfrow = c(1, 2))
# plot(thetaEst, type = "b", lty = 1, pch = 20, main = expression(theta), ylim = c(0,6))
# abline(h = c(trueTheta + (-5:5) * 2 * pi), lwd = 2, col = "orange")
# plot(sigmaEst, type = "b", lty = 1, pch = 20, main = expression(sigma))
# 
# 
# # Exploring particle Weights ----------------------------------------------
# 
# plotPartsWeights <- function(index) {
#   oneRes$Weights
#   parts <- oneRes$Particles[, index]
#   ws <- oneRes$Weights[, index]
#   trueState <- SINEprocess[index, "hiddenStates"]
#   obs <- SINEprocess[index, "observations"]
#   xlim = range(c(parts, trueState, obs))
#   plot(x = parts, y = ws, xlim = xlim, ylim = c(0, 1), pch = 20, main = paste("Time", index),
#        xaxt = "n", ylab = "", xlab = "") 
#   abline(v = c(trueState, obs), col = c("red", "blue"), lty = 2, lwd = 2)
# }
# par(mfrow = c(4, 4), mar = c(0.1, 3, 4, 0.1))
# sapply(1:n, plotPartsWeights)