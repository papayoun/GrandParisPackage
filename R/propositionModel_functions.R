loadModule("ProposalModel_Module", TRUE)


#' @title Simulate a SINE process
#' @name Prop_evalTrans
#' @description  function that simulates a SINE model trajectory using exact algorithm
#' @param theta Theta parameter of the model, default to pi
#' @param sigma2 observation process variance
#' @param oldParticles Starting particles
#' @param newParticles vector of simulation times
#' @return a vector same length as times
#' @export
Prop_evalTrans <- function(oldParticles, newParticles, obs, delta = 1, RW =1, theta = 3.14, sigma2 = 1){
  model <- new(ProposalSINEModel, RW, theta, sigma2, TRUE, TRUE)
  transDens <- model$transDens(oldParticles, newParticles, 0, delta, 30, TRUE)
  propDens <- model$propDens(oldParticles, newParticles, obs, delta, FALSE)
  obsDens <- model$obsDens(newParticles, obs)
  cbind(trans = transDens, prop = propDens, obs = obsDens)
}