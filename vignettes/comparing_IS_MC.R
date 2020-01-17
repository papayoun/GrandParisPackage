
# Cleaning ----------------------------------------------------------------

rm(list = ls())

# Packages ----------------------------------------------------------------

library(GrandParisPackage)
library(parallel)
library(tidyverse)
library(tictoc)
# Simulating data ---------------------------------------------------------

my_seed <- 333 # For all experiments, the random seed

set.seed(my_seed) # For reproducibility 
trueTheta <- pi / 4; trueSigma2 <- 1;
n <- 11; times <- seq(from = 0, to = 5, length = n);
SINEprocess <- SINE_simulate(theta = trueTheta, sigma2 = trueSigma2, 
                             x0 = 10, times = times)
observations <- SINEprocess[, "observations"]


# Experiment_function -----------------------------------------------------

# A function that computes the approximated estimate of E[X_0 | Y_{0:n}]
# for different tilde(N)
# It does it 50 times

compare_on_ntilde <- function(n_particle = 100, n_tilde = NULL){
  if(is.null(n_tilde)){
    stop("You must provide ntilde!")
  }
  n_rep <- 400
  AR <- do.call(rbind, 
                mclapply(1:n_rep, function(i){
                  tic()
                  res <- GrandParisPackage:::E_track(observations = observations, 
                                                     ind_tracked = 0, 
                                                     observationTimes = times, 
                                                     thetaStart = trueTheta, 
                                                     particleSize = n_particle, 
                                                     backwardSampleSize = n_tilde,
                                                     sigma2Start = trueSigma2, 
                                                     nIterations = 1)
                  my_tictoc <- toc()
                  dur <- my_tictoc$toc - my_tictoc$tic
                  names(dur) <- NULL
                  data.frame(N = n_particle,
                             N_tilde = n_tilde, Xhat = res, Time = dur, 
                             Method = factor("AR", 
                                             levels = c("AR", "IS")))
                }, mc.cores = 10))
  IS <- do.call(rbind, 
                mclapply(1:n_rep, function(i){
                  tic()
                  res <- GrandParisPackage:::E_track_IS(observations = observations, 
                                                     ind_tracked = 0, 
                                                     observationTimes = times, 
                                                     thetaStart = trueTheta, 
                                                     particleSize = n_particle, 
                                                     backwardSampleSize = n_tilde,
                                                     sigma2Start = trueSigma2, 
                                                     nIterations = 1)
                  my_tictoc <- toc()
                  dur <- my_tictoc$toc - my_tictoc$tic
                  names(dur) <- NULL
                  data.frame(N = n_particle,
                             N_tilde = n_tilde, Xhat = res, Time = dur, 
                             Method = factor("IS", 
                                             levels = c("AR", "IS")))
                }, mc.cores = 10))
  return(rbind.data.frame(AR, IS))
}


# Obtaining estimates -----------------------------------------------------

set.seed(my_seed) # For reproducibility 
my_ntildes <- c(2, 5, 10, 20, 30)


# Run only once (experiment) ----------------------------------------------

res  <- lapply(my_ntildes, function(n_t)
  compare_on_ntilde(n_particle = 100, n_tilde = n_t))
res_df <- do.call(rbind.data.frame, res)
save(res_df, file = "vignettes/compare_IS_MC_ntilde_N100_400rep.RData")

# Plotting results --------------------------------------------------------

load("vignettes/compare_IS_MC_ntilde_N100_400rep.RData")

main_plot <- ggplot(res_df, aes(x = factor(N_tilde), fill = Method)) +
  labs(x = expression(tilde(N))) +
  theme(text = element_text(family = "Symbola", size = 24))

p1 <- main_plot + 
  geom_boxplot(aes(y = Time)) + 
  scale_y_continuous(trans = "log10") +
  labs(y = "Comput. time (seconds)")
p2 <-  main_plot + 
  geom_boxplot(aes(y = Xhat)) +
  labs(y = expression(hat("\U1d53c")~"["~X[0]~"|"~Y[0:n]~"]"))
my_plot <- gridExtra::grid.arrange(p1, p2)

# Then save it in pdf format 10 x 10 inches
# as comparing_IS_AR.pdf

# Controlling N -----------------------------------------------------------

compare_on_n <- function(n_particle = 100, alphas = NULL, frac = NULL){
  n_tilde_AR <-  2
  levels_IS <- paste0("IS_a_", alphas, "_p_", frac)
  n_rep <- 50
  AR <- do.call(rbind, 
                mclapply(1:n_rep, function(i){
                  tic()
                  res <- GrandParisPackage:::E_track(observations = observations, 
                                                     ind_tracked = 0, 
                                                     observationTimes = times, 
                                                     thetaStart = trueTheta, 
                                                     particleSize = n_particle, 
                                                     backwardSampleSize = n_tilde_AR,
                                                     sigma2Start = trueSigma2, 
                                                     nIterations = 1)
                  my_tictoc <- toc()
                  dur <- my_tictoc$toc - my_tictoc$tic
                  data.frame(N = n_particle,
                             N_tilde = n_tilde_AR, Xhat = res, Time = dur, 
                             Method = factor("AR", 
                                             levels = c("AR", levels_IS)))
                }, mc.cores = 10))
  IS <- do.call(rbind.data.frame,
                mapply(function(my_alpha, my_prop){
                  n_tilde_IS = ceiling(my_prop * n_particle^my_alpha)
                  do.call(rbind.data.frame,
                          mclapply(1:n_rep, function(i){
                    tic()
                    res <- GrandParisPackage:::E_track_IS(observations = observations,
                                                          ind_tracked = 0,
                                                          observationTimes = times,
                                                          thetaStart = trueTheta,
                                                          particleSize = n_particle,
                                                          backwardSampleSize = n_tilde_IS,
                                                          sigma2Start = trueSigma2,
                                                          nIterations = 1)
                    my_tictoc <- toc()
                    dur <- my_tictoc$toc - my_tictoc$tic
                    names(dur) = NULL
                    data.frame(N = n_particle,
                               N_tilde = n_tilde_IS, Xhat = res, Time = dur,
                               Method = factor(paste0("IS_a_", my_alpha, "_p_", my_prop),
                                               levels = c("AR", levels_IS)))
                  }, mc.cores = 10))
                }, alphas, frac, SIMPLIFY = F))
  return(rbind.data.frame(AR, IS))
}

my_n_particles <- c(50, 100, 200, 500, 1000, 2000)
set.seed(my_seed)
res_npart  <- lapply(my_n_particles, 
                 function(my_n) compare_on_n(n_particle = my_n, 
                                             alphas = c(0.5, 0.6, 1), 
                                             frac = c(1, 1, 0.1))) 
res_df_npart <- do.call(rbind.data.frame, res_npart)
save(res_df_npart, file = "vignettes/compare_IS_MC_vary_npart.RData")

# Plotting results --------------------------------------------------------

load("vignettes/compare_IS_MC_vary_npart.RData")

my_labels <- expression("AR,"~tilde(N)==2, "IS,"~tilde(N)== N^0.5,
                        "IS,"~tilde(N)== N^0.6, "IS,"~tilde(N)== N/10)
main_plot <- res_df_npart %>% 
  mutate(Method = factor(Method, labels = my_labels)) %>% 
  ggplot(aes(x = factor(N), fill = Method)) +
  labs(x = expression(N)) +
  theme(text = element_text(family = "Symbola", size = 24)) +
  scale_fill_discrete(labels = my_labels) +
  theme(legend.text.align = 0)

p1 <- main_plot + 
  geom_boxplot(aes(y = Time)) + 
  scale_y_continuous(trans = "log10") +
  labs(y = "Comput. time (seconds)")
p2 <-  main_plot + 
  geom_boxplot(aes(y = Xhat)) +
  labs(y = expression(hat("\U1d53c")~"["~X[0]~"|"~Y[0:n]~"]"))
my_plot <- gridExtra::grid.arrange(p1, p2)

ggplot(res_df_npart, aes(x = N, y = Time)) + 
  geom_point() + geom_smooth() +
  facet_wrap(~Method, scales = "free_y") +
  labs(x = expression(tilde(N)),
       y = "Comput. time (seconds)")
