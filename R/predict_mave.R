#' Make predictions based on a 'mave' object
#'
#' This method make predictions based the reduced dimension of data using nadaraya-watson method.
#'
#' @param object the object of class 'mave'
#' @param newx Matrix of the new data to be predicted
#' @param dim the dimension of central space or central mean space. The matrix of the original data will be
#' multiplied by the matrix of dimension reduction directions of given dimension. Then the prediction will be
#' made based on the data of given dimensions.
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
#'
#' yp = predict(dr,x.test,3)
#' #mean error
#' mean(abs(yp-y.test)/y.test)
#'
#'@method predict mave
#'@export

predict.mave<-function(object, newx, dim,...){

  dir <- mave.dir(object,dim)
  x <- as.matrix(object$x%*%dir)
  newx <- as.matrix(newx%*%dir)
  allx <-rbind(x,newx)
  allx <- scale(allx)
  n <- nrow(x)
  np <- nrow(newx)
  x <- as.matrix(allx[1:n,])
  newx <- as.matrix(allx[(n+1):(n+np),])
  yp <- PredictCpp(x,newx,object$y)
  return(yp)

}

