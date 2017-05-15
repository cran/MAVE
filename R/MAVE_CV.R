#' Select best direction using cross-validation
#'
#'This function selects the dimension of the central (mean)  space
#'based on the calculation of MAVE  using cross-validation method.
#'
#' @param dr the result of MAVE function
#' @param max.dim the maximum dimension for cross-validation.
#'
#' @return dr.dim contains all information in dr plus cross-validation values of corresponding
#' direction

#' @export
#'
#' @examples
#'  x <- matrix(rnorm(400*5),400,5)
#'  b1 <- matrix(c(1,1,0,0,0),5,1)
#'  b2 <- matrix(c(0,0,1,1,0),5,1)
#'  eps <- matrix(rnorm(400),400,1)
#'  y <- x%*%b1 + (x%*%b2)*eps
#'
#'  #seleted dimension of central space
#'  dr.cs <- mave(y~x,method='csmave')
#'  dr.cs.dim <- mave.dim(dr.cs)
#'
#'  #seleted dimension of central mean space
#'  dr.mean <- mave(y~x,method='meanmave')
#'  dr.mean.dim <- mave.dim(dr.mean)

mave.dim<-function(dr, max.dim=10){

  max.dim = min(max.dim, ncol(dr$x))
  which.dim <- 1:min(dr$max.dim, max.dim)
  if(length(which.dim)==0){
    warning('No CV is computed given the dr$max.dim and the argument max.dim')
    return(dr)
  }
  if(is.null(dr$cv)){
    dr$cv <- rep(Inf, max.dim)
  }
  for(dim in which.dim){
    dr$cv[dim] <- CVfastCpp(dr$x%*%mave.dir(dr,dim),dr$ky)
  }
  class(dr) <- 'mave.dim'
  dr$call<-match.call()

  return(dr)
}
