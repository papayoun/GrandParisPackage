library(GrandParisPackage)

set.seed(123)
trueTheta <- pi; trueSigma2 <- 1;
n <- 1000; times <- seq(0, n - 1, length = n);
SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 0, times = times)
observations <- SINEprocess[, "observations"]
estTheta = T; estSig = F; updateOrders <- rep(T, n)
gradientSteps   <- matrix(1, ncol = 2, nrow = n)

foo <- function(seed){
  set.seed(seed)
  thetaStart <- runif(1, 0, 2 * pi)
  Res <- fastTangOR(observations, times,
                    thetaModel = thetaStart, sigma2 = trueSigma2,
                    updateOrders = updateOrders, gradientSteps = gradientSteps, all = F, estimateSigma2 = F)
  vecEst <- Res$Estimates[,1]
  save(vecEst, file = paste0("vignettes/simulatedResults/res_seed", seed, ".RData"))
  return(vecEst)
}
bigRes <- do.call(cbind, mclapply(1:48, foo, mc.cores = detectCores()))

library(parallel)
library(extrafont)
par(mfrow = c(2, 1), family = "LM Roman 10")
par(mar = c(mar=c(0, 4.1, 4.1, 0.5)))
matplot(times, SINEprocess[, 1:2], 
        xaxt = "n", type = c("p", "b"), pch = c(18, 20), col = c("blue", "red"),
        ylab = "Process values", cex = 0.5, xlab = "", cex.lab = 1.5)
mtext(side = 3, line = 2, expression(Estimation~of~theta~with~known~sigma^2), cex = 2)
legend(x = "topright", bty = "n", pch = c(18, 20), col = c("blue", "red"), cex = 1.5, 
       legend = c("Observed", "Hidden"))
par(mar = c(mar=c(4.1, 4.1, 0, 0.5)))
matplot(times, bigRes, lty = 1, type = "l", xlab = "", ylab = "Estimated value", xaxt = "n",
        cex.lab = 1.5)
abline (h = trueTheta, lwd = 2, col = "black", lty = 2)
axis(side = 1, at = seq(0, 1000, by = 100))
mtext(side = 1, line = 2.5, "Time", cex = 2)
