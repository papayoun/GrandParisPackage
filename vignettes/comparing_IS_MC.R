library(GrandParisPackage)
library(parallel)
# seed <- 123
# set.seed(seed)
# trueTheta <- pi / 4; trueSigma2 <- 1;
# n <- 5; times <- seq(0, to = 2, length = n);
# SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
# observations <- SINEprocess[, "observations"]
# n_case <- 6
# allRes <- mclapply(1:n_case, function(i){
#   thetaStart <- runif(1, 0, 2 * pi)
#   sigma2Start <- runif(1, 0.2, 2)
#   myTry <- GEM(observations = observations, observationTimes = times, thetaStart = thetaStart, 
#                particleSize = 5, backwardSampleSize = 50,
#                sigma2Start = sigma2Start, nIterations = 1, nModels = 1)
#   print("%%%%%%%%%%%% AUTRE %%%%%%%%%%%%")
#   myTry2 <- GEM_IS(observations = observations, observationTimes = times, thetaStart = thetaStart, 
#                particleSize = 5, backwardSampleSize = 50,
#                sigma2Start = sigma2Start, nIterations = 1, nModels = 1)
#   colnames(myTry) = c("theta", "sigma2")
#   myTry
# }, mc.cores = min(n_case, detectCores()))
# 
# thetaEst <- sapply(allRes, function(x) x[,"theta"]) %% (2 * pi)
# sigmaEst <- sapply(allRes, function(x) x[, "sigma2"])
# par(mfrow = c(2, 1))
# matplot(thetaEst, col = 1, type = "b", pch = 3, lty = 1, main = expression(theta))
# abline(h = trueTheta, col = "red", lty = 2, lwd = 3)
# matplot(sigmaEst, col = 1, type = "b", pch = 3, lty = 1, main = expression(sigma^2))
# abline(h = trueSigma2, col = "red", lty = 2, lwd = 3)


# Comparing computation time ----------------------------------------------

rm(list = ls())
# library(GrandParisPackage)
trueTheta <- pi / 4; trueSigma2 <- 1;
n <- 5; times <- seq(0, to = 2, length = n);
SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
observations <- SINEprocess[, "observations"]
thetaStart <- trueTheta
sigma2Start <- trueSigma2
# library(microbenchmark)

foo <- function(n_tilde){
  n_rep <- 100
  AR <- do.call(rbind, 
          lapply(1:n_rep, function(i){
            t1 <- Sys.time()
            res <- GrandParisPackage:::E_track(observations = observations, ind_tracked = 0, 
                                        observationTimes = times, thetaStart = thetaStart, 
                                        particleSize = 100, backwardSampleSize = n_tilde,
                                        sigma2Start = sigma2Start, nIterations = 1)
            dur <- as.numeric(Sys.time() - t1)
            data.frame(N = n_tilde, Xhat = res, Time = dur, 
                       Method = factor("AR", levels = c("IS", "AR")))
          }))
  IS <- do.call(rbind, 
                lapply(1:n_rep, function(i){
                  t1 <- Sys.time()
                  res <- GrandParisPackage:::E_track_IS(observations = observations, ind_tracked = 0, 
                                                     observationTimes = times, thetaStart = thetaStart, 
                                                     particleSize = 100, backwardSampleSize = n_tilde,
                                                     sigma2Start = sigma2Start, nIterations = 1)
                  dur <- as.numeric(Sys.time() - t1)
                  data.frame(N = n_tilde, Xhat = res, Time = dur, 
                             Method = factor("IS", levels = c("IS", "AR")))
                }))
  return(rbind.data.frame(AR, IS))
}
library(tidyverse)
seed <- 1
set.seed(seed)
res  <- mclapply(c(2, 5, 10, 15, 20, 30, 50, 60, 70, 80, 100), foo,
                 mc.cores = 11) 
res_df <- do.call(rbind.data.frame, res)
save(res_df, file = "vignettes/compare_IS_MC_ntilde.RData")


# Plotting results --------------------------------------------------------

load("vignettes/compare_IS_MC_ntilde.RData")

p1 <- ggplot(res_df, aes(x = factor(N), y = Time, col = Method)) + 
  geom_point() + scale_y_continuous(trans = "log") + 
  labs(x = expression(tilde(N)),
       y = "Comput. time (seconds, log scale)")
p2 <-  res_df %>% 
  ggplot(aes(x = factor(N), y = Xhat, col = Method, fill = Method)) + 
  geom_boxplot() +
  labs(x = expression(tilde(N)),
       y = expression(hat(E)(X[0]~"|"~Y)))
gridExtra::grid.arrange(p1, p2)

ggplot(res_df, aes(x = N, y = Time)) + 
  geom_point() + geom_smooth() +
  facet_wrap(~Method, scales = "free_y") +
  labs(x = expression(tilde(N)),
       y = "Comput. time (seconds)")


# Controlling N -----------------------------------------------------------

rm(list = ls())
# library(GrandParisPackage)
trueTheta <- pi / 4; trueSigma2 <- 1;
n <- 5; times <- seq(0, to = 2, length = n);
SINEprocess <- SINE_simulate(trueTheta, trueSigma2, 10, times = times)
observations <- SINEprocess[, "observations"]
thetaStart <- trueTheta
sigma2Start <- trueSigma2
# library(microbenchmark)

foo_npart <- function(n_particle){
  n_rep <- 100
  AR <- do.call(rbind, 
                lapply(1:n_rep, function(i){
                  t1 <- Sys.time()
                  res <- GrandParisPackage:::E_track(observations = observations, 
                                                     ind_tracked = 0, 
                                                     observationTimes = times, 
                                                     thetaStart = thetaStart, 
                                                     particleSize = n_particle, 
                                                     backwardSampleSize = 2,
                                                     sigma2Start = sigma2Start, 
                                                     nIterations = 1)
                  dur <- as.numeric(Sys.time() - t1)
                  data.frame(N = n_particle, Xhat = res, Time = dur, 
                             Method = factor("AR", levels = c("IS", "AR")))
                }))
  IS <- do.call(rbind, 
                lapply(1:n_rep, function(i){
                  t1 <- Sys.time()
                  res <- GrandParisPackage:::E_track_IS(observations = observations, 
                                                        ind_tracked = 0, 
                                                        observationTimes = times, 
                                                        thetaStart = thetaStart, 
                                                        particleSize = n_particle, 
                                                        backwardSampleSize = ceiling(n_particle * 0.1),
                                                        sigma2Start = sigma2Start, 
                                                        nIterations = 1)
                  dur <- as.numeric(Sys.time() - t1)
                  data.frame(N = n_particle, Xhat = res, Time = dur, 
                             Method = factor("IS", levels = c("IS", "AR")))
                }))
  return(rbind.data.frame(AR, IS))
}
library(tidyverse)
seed <- 1
set.seed(seed)
res  <- mclapply(floor(c(seq(50, 400, length = 8),
                         500, 1000, 2000)), foo_npart,
                 mc.cores = 11) 
res_df_npart <- do.call(rbind.data.frame, res)
save(res_df_npart, file = "vignettes/compare_IS_MC_npart_vary_ntilde.RData")


# Plotting results --------------------------------------------------------

load("vignettes/compare_IS_MC_npart_vary_ntilde.RData")

p1 <- ggplot(res_df_npart, aes(x = N, y = Time, col = Method)) + 
  geom_point() + scale_y_continuous(trans = "log") + 
  labs(x = "N",
       y = "Comput. time (seconds, log scale)")
p2 <-  res_df_npart %>% 
  ggplot(aes(x = factor(N), y = Xhat, col = Method, fill = Method)) + 
  geom_boxplot() +
  labs(x = expression(tilde(N)),
       y = expression(hat(E)(X[0]~"|"~Y)))
gridExtra::grid.arrange(p1, p2)

ggplot(res_df_npart, aes(x = N, y = Time)) + 
  geom_point() + geom_smooth() +
  facet_wrap(~Method, scales = "free_y") +
  labs(x = expression(tilde(N)),
       y = "Comput. time (seconds)")
