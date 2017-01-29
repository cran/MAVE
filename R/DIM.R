#' Select best direction using cross-validation
#'
#'This function selects the dimension of the central (mean)  space
#'based on the calculation of MAVE  using cross-validation method
#'
#' @param rd the result of MAVE function
#'
#' @return rd contains cross-validation values of corresponding
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
#'  rd1 <- MAVE(x,y,'csmave')
#'  rd2 <- MAVE(x,y,'meanmave')
#'
#'  rd1 <- DIM(rd1)
#'  rd2 <- DIM(rd2)

DIM<-function(rd){
  dir1D=list2M1d(rd$dir);
  rd$cv<-CVfast(rd$x,rd$ky,dir1D)
  rd$call<-match.call()
  return(rd)
}
