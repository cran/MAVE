#' Dimension reduction
#'
#' This function provides several methods to estimate the central space or central mean space of y on x.
#'  It returns the matrix of central space or central mean space for different dimensions and contains
#'  other information used for dimension selection by \code{\link{DIM}}.
#'
#' @param x The design matrix.
#' @param y The respond vector.
#' @param method This parameter specify which method will be used in dimension reduction. It provides
#' five methods, including "csMAVE","csOPG","meanOPG","meanMAVE","KSIR"
#' by default, method = 'csOPG'
#' \itemize{
#' \item 'meanOPG' and 'meanMAVE' estimate dimension reduction space
#'            for conditional mean
#' \item 'csMAVE' and 'csOPG' estimate the central dimension reduction
#'            space
#' \item 'kSIR' is a kernel version of sliced inverse regression (Li, 1991). It is fast, but
#'            with poor efficiency
#' }
#'
#' @return rd rd is a list which contains:
#' \itemize{
#' \item dir: dir[[d]] is the central space with d-dimension
#'         d = 1, 2, ..., p reduced direction of different dimensions
#' \item x: parameter used for DIM for selection
#' \item ky: parameter used for DIM for selection
#' }
#'
#' @export
#'
#' @seealso \code{\link{DIM}} for dimension selection
#'
#' @references Li K C. Sliced inverse regression for dimension reduction[J]. Journal of the American Statistical Association, 1991, 86(414): 316-327.
#' @references Xia Y, Tong H, Li W K, et al. An adaptive estimation of dimension reduction space[J]. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 2002, 64(3): 363-410.
#' @references Xia Y. A constructive approach to the estimation of dimension reduction directions[J]. The Annals of Statistics, 2007: 2654-2690.
#' @references Wang H, Xia Y. Sliced regression for dimension reduction[J]. Journal of the American Statistical Association, 2008, 103(482): 811-821.
#'
#' @examples
#'  x <- matrix(rnorm(200*5),200,5)
#'  b1 <- matrix(c(1,1,0,0,0),5,1)
#'  b2 <- matrix(c(0,0,1,1,0),5,1)
#'  eps <- matrix(rnorm(200),200,1)
#'  y <- x%*%b1 + (x%*%b2)*eps
#'
#'  rd1 <- MAVE(x,y,'csopg')
#'  rd2 <- MAVE(x,y,'meanopg')
#'  rd3 <- MAVE(x,y,'ksir')
#'
#' @useDynLib MAVE
#' @importFrom Rcpp evalCpp


MAVE<-function(x,y,method='CSOPG'){
  method=toupper(method);
  methodvec=c('CSOPG','CSMAVE','MEANOPG','MEANMAVE','KSIR')
  if(is.matrix(y)&&ncol(y)!=1){
    stop('Sorry, the package now couldn\'t work with with multi-dimension response. The function will appear later.')
  }
  if(!(method %in% methodvec)){
    stop('method should be one of CSMAVE, CSOPG, MEANOPG, MEANMAVE, KSIR')
  }
  if(nrow(x)!=length(y)){
    stop('the row of x and the row of y is not compatible')
  }
  if(length(x)==length(y)){
    stop('x is one dimensional, no need do dimension reduction')
  }
  rd=MAVEfast(x,y,method)
  rd$call=match.call()
  rd$method=method
  rd$dir=M3d2list(rd$BB)
  rd$BB=NULL
  attr(rd,'class')<-'mave'
  return(rd)
}
