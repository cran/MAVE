#' Select best direction using cross-validation
#'
#'This function selects the dimension of the central (mean)  space
#'based on the calculation of MAVE  using cross-validation method
#'
#' @param rd the result of MAVE function
#'
#' @return dim contains cross-validation values of corresponding
#' direction

#' @export
#'
#' @examples
#'  x <- matrix(rnorm(200*5),200,5)
#'  b1 <- matrix(c(1,1,0,0,0),5,1)
#'  b2 <- matrix(c(0,0,1,1,0),5,1)
#'  eps <- matrix(rnorm(200),200,1)
#'  y <- x%*%b1 + (x%*%b2)*eps
#'
#'  #seleted dimension of central space
#'  rd.cs <- MAVE(x,y,'csmave')
#'  dim.cs <- DIM(rd.cs)
#'  
#'  #seleted dimension of central mean space
#'  rd.mean <- MAVE(x,y,'meanmave')
#'  dim.mean <- DIM(rd.mean)

DIM<-function(rd){
  dir1D=list2M1d(rd$dir);
  rd$cv<-CVfastCpp(rd$x,rd$ky,dir1D)
  rd$call<-match.call()
  rd$best_dir=rd$dir[[which.min(rd$cv)]]
  return(rd)
}
