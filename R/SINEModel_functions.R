loadModule("SINEModel_Module", TRUE)

#' @title Simulate a SINE process
#' @name SINE_simulate
#' @description  function that simulates a SINE model trajectory using exact algorithm
#' @param theta Theta parameter of the model, default to pi
#' @param sigma2 observation process variance
#' @param x0 starting value
#' @param times vector of simulation times
#' @return a vector same length as times
#' @export
SINE_simulate <- function(theta = pi, sigma2 = 1, x0 = 0, times = seq(0, 1, length.out = 10)){
  model <- new(SINE_POD, theta, sigma2)
  return(do.call(cbind, model$rSINE_POD(x0, times)))
}

#' @title compute density  in the SINE model
#' @name  SINE_density
#'
#' @param theta Theta parameter of the model, default to pi
#' @param X Process values
#' @param times vector of simulation times
#' @param log should the log be returned
#' @param MC.size Monte Carlo effort to compute density
#' @return a vector having length = length(X) - 1 giving the density of each segment
#' @export
SINE_density <- function(X, times, theta = pi, MC.size = 30, log = F){
  if(!is.null(dim(X)))
    stop("X should be a vector")
  n <- length(X)
  if(n < 2)
    stop("n should be of length at least 2")
  model <- new(SINE_POD, theta, 1)
  if(log){
    foo <- function(i){
      model$logDensity(X[i], X[i + 1], times[i], times[i + 1], MC.size)
    }
  }
  else{
    foo <- function(i){
      model$density(X[i], X[i + 1], times[i], times[i + 1], MC.size)
    }
  }
  return(sapply(1:(n - 1), foo))
}

#' @title compute logdensity gradient in the SINE model
#' @name  SINE_logDensityGradient 
#'
#' @param theta Theta parameter of the model, default to pi
#' @param X Process values
#' @param times vector of simulation times
#' @param MC.size Monte Carlo effort to compute density
#' @return a vector having length = length(X) - 1 giving the density of each segment
#' @export
SINE_logDensityGradient <- function(X, times, theta = pi, MC.size = 30){
  if(!is.null(dim(X)))
    stop("X should be a vector")
  n <- length(X)
  if(n < 2)
    stop("n should be of length at least 2")
  model <- new(SINE_POD, theta, 1)
  foo <- function(i){
    model$gradLogDensity(X[i], X[i + 1], times[i], times[i + 1], MC.size)
  }
  return(sapply(1:(n - 1), foo))
}
