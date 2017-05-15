#' Directions of CS or CMS of given dimension
#'
#' This function returns the basis matrix of CS or CMS of given dimension
#'
#' @param dr the output of \code{\link{mave}} or the output of \code{\link{mave.dim}}
#' @param dim the dimension of CS or CMS. The value of dim should be given when the class of the
#' argument dr is mave. When the class of the argument dr is mave.dim and dim is not given, the
#' function will return the basis matrix of CS or CMS of dimension selected by \code{\link{mave.dim}}
#'
#' @return dir the matrix of CS or CMS of given dimension
#'
#' @export
#'
#' @examples
#' x <- matrix(rnorm(400),100,4)
#' y <- x[,1]+x[,2]+as.matrix(rnorm(100))
#' dr <- mave(y~x)
#' dir3 <- mave.dir(dr,3)
#'
#' dr.dim <- mave.dim(dr)
#' dir3 <- mave.dir(dr.dim,3)
#' dir.best <- mave.dir(dr.dim)
#'

mave.dir<-function(dr, dim = NULL){

    if(class(dr)=='mave'){
      if(!is.null(dim)){
        if(dim<=dr$max.dim){
          return(as.matrix(dr$dir[[dim]]))
        }
        else{
            arg='central space'
            if(substr(dr$method,1,2)!='CS') arg='central mean space'
            stop(c('The ', arg, ' of dimension ',dim,' haven\'t been computed, please compute it first'))
        }
      }
      else{
        stop('dim should be given.')
      }
    }
    if(class(dr)=='mave.dim'){
      if(is.null(dim)){
        return(as.matrix((dr$dir[[which.min(dr$cv)]])))
      }
      else{
        return(as.matrix(dr$dir[[dim]]))
      }
    }

}
