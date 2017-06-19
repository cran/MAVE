#' Make predictions based on a 'mave' object
#'
#' This method make predictions based the reduced dimension of data using \code{\link{mars}} function.
#'
#' @param object the object of class 'mave'
#' @param newx Matrix of the new data to be predicted
#' @param dim the dimension of central space or central mean space. The matrix of the original data will be
#' multiplied by the matrix of dimension reduction directions of given dimension. Then the prediction will be
#' made based on the data of given dimensions.
#' @param ... further arguments passed to \code{\link{mars}} function such as degree.
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
#' yp = predict(dr,x.test,dim=3,degree=2)
#' #mean error
#' mean((yp-y.test)^2)
#'
#'@method predict mave
#'@export

predict.mave<-function(object, newx, dim, ...){

  n0 = nrow(object$x)
  n1 = nrow(newx)

  x.train.mave <- mave.data(object, object$x, dim)
  x.test.mave <- mave.data(object, newx, dim)

  thresh0 <- 0.5/sqrt(n0)
  thresh <- 0.5/sqrt(n1)
  y.pred <- matrix(Inf, n1, 1)
  y.ub <- max(object$y) + 0.5*(max(object$y)-min(object$y))
  y.lb <- min(object$y) - 0.5*(max(object$y)-min(object$y))

  idx <- which((y.pred>y.ub)|(y.pred<y.lb))
  counter <- 0
  while(length(idx)>0){
    if(any(substr(rep("thresh", 5), 1, 2:6) %in% names(list(...)))){
      model <- mda::mars(x.train.mave,object$y,...)
    }
    else{
      model <- mars(x.train.mave,object$y,thresh=thresh,...)
    }
    y.pred[idx] <- predict(model,x.test.mave[idx,,drop=F])
    idx <- which( (y.pred>y.ub) | (y.pred<y.lb) )
    thresh = thresh + thresh0

    counter <- counter+1
    if(counter>100){
      break
    }
  }
  return(as.matrix(y.pred))

}
