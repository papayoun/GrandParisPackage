loadModule("MyModule", TRUE)

#' function that calls the class of Rcpp
#'
#' @export
createInstance <- function(class, init = 2){
  res <- new(class, init)
  return(res)
}