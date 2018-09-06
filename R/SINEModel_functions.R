loadModule("SINEModel_Module", TRUE)

#' @title Simulate a SINE process
#' @name SINE_simulate
#' @description  function that simulates a SINE model trajectory using exact algorithm
#' @param theta Theta parameter of the model, default to pi
#' @param x0 starting value
#' @param times vector of simulation times
#' @return a vector same length as times
#' @export
SINE_simulate <- function(theta = pi, x0 = 0, times = seq(0, 1, length.out = 10)){
  model <- new(SINEModel, theta)
  return(model$rSINE(x0, times, 10000))
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
SINE_density <- function(theta = pi, X, times, MC.size = 30, log = F){
  if(!is.null(dim(X)))
    stop("X should be a vector")
  n <- length(X)
  if(n < 2)
    stop("n should be of length at least 2")
  model <- new(SINEModel, theta)
  if(log){
    foo <- function(i){
      model$logDensity(X[i], X[i + 1], times[i], times[i + 1], MC.size, 1000)
    }
  }
  else{
    foo <- function(i){
      model$density(X[i], X[i + 1], times[i], times[i + 1], MC.size, F)
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
SINE_logDensityGradient <- function(theta = pi, X, times, MC.size = 30){
  if(!is.null(dim(X)))
    stop("X should be a vector")
  n <- length(X)
  if(n < 2)
    stop("n should be of length at least 2")
  model <- new(SINEModel, theta)
  foo <- function(i){
    model$gradLogDensity(X[i], X[i + 1], times[i], times[i + 1], MC.size, 1000)
  }
  return(sapply(1:(n - 1), foo))
}
