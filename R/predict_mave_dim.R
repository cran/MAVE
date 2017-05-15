#' Make predictions based on a 'mave.dim' object
#'
#' This method make predictions based the best dimension of data using nadaraya-watson method. The best
#' dimension is selected by the \code{\link{mave.dim}} method or the dimension which is given.
#'
#' @param object the object of class 'mave'
#' @param newx Matrix of the new data to be predicted
#' @param dim the dimension of central space or central mean space. The matrix of the original data will be
#' multiplied by the matrix of dimension reduction directions of given dimension. Then the prediction will be
#' made based on the data of given dimensions. The default value is 'dim.min'. In this case, the dimension
#' selected by \code{\link{mave.dim}} is used. When dimension is given, the function is same with \code{\link{predict.mave}}
#' @param ... further arguments passed to or from other methods.
#' @return the prediced response of the new data
#' @examples
#'
#' X = matrix(rnorm(10000),1000,10)
#' beta1 = as.matrix(c(1,1,1,1,0,0,0,0,0,0))
#' beta2 = as.matrix(c(0,0,0,1,1,1,1,1,0,0))
#' err = as.matrix(rnorm(1000))
#' Y = X%*%beta1+X%*%beta2+err
#'
#' train = sample(1:1000)[1:500]
#' x.train = X[train,]
#' y.train = as.matrix(Y[train])
#' x.test = X[-train,]
#' y.test = as.matrix(Y[-train])
#'
#' dr = mave(y.train~x.train, method = 'meanopg')
#' dr.dim = mave.dim(dr)
#'
#' yp = predict(dr.dim,x.test)
#' #mean error
#' mean(abs(yp-y.test)/y.test)

#'@method predict mave.dim
#'@export


predict.mave.dim<-function(object,newx,dim='dim.min',...){
  if(dim=='dim.min'){
    dir <- mave.dir(object)
  }else{
    dir <- mave.dir(object,dim)
  }

  x <- object$x%*%dir
  newx <- newx%*%dir
  allx <-rbind(x,newx)
  allx <- scale(allx)
  n <- nrow(x)
  np <- nrow(newx)
  x <- as.matrix(allx[1:n,])
  newx <- as.matrix(allx[(n+1):(n+np),])
  yp <- PredictCpp(x,newx,object$y)
  return(yp)

}
