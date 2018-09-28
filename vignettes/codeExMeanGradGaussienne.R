rm(list = ls())
library(parallel)
echantillon <- rnorm(100)
gradlogObs <- function(x, variance){
  (x^2 / variance - 1) / (2 * variance)
}
gradEch <- function(ech, variance){
  mean(gradlogObs(ech, variance)) * dnorm(sqrt(mean(echantillon^2)), 0, sqrt(variance))
}
gradObs <- function(x, variance){
  gradlogObs(x, variance) * dnorm(x, 0, sqrt(variance))
}
variances <- seq(0.0001, 10, length.out = 500)

expectedLog <- do.call(c, mclapply(variances, function(var) mean(gradlogObs(echantillon, var)),
                              mc.cores = detectCores()))
expectedRaw <- do.call(c, mclapply(variances, function(var) mean(gradObs(echantillon, var)),
                                   mc.cores = detectCores()))
echRaw <- do.call(c, mclapply(variances, function(var) gradEch(echantillon, var)))

par(mfrow = c(3, 1), mar = c(3, 3, 1, 0.1))
plot(variances, expectedLog, ylim = c(-0.1, 1), type = "b", pch = 3, lwd = 2, cex = 0.1)
abline(h = 0); abline(v = 1)
plot(variances, expectedRaw,  ylim = c(-0.1, 1), type = "b", pch = 3, lwd = 2, cex = 0.1)
abline(h = 0); abline(v = 1)
plot(variances, echRaw,  ylim = c(-0.1, 1), type = "b", pch = 3, lwd = 2, cex = 0.1)
abline(h = 0); abline(v = 1)

