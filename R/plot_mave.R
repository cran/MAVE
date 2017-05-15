#' Plot of mave object
#'
#' Plot the scatterplot of given dimension directions and reponse variables.
#' @aliases plot.mave
#' @aliases plot.mave.dim
#' @rdname plot.mave
#' @param x the object returned by mave
#' @param dim the dimension
#' @param plot.method the method for plotting scatter plot. The default is 'pairs'
#' @param ... arguments passed to the plot.method.
#' @examples
#' x = matrix(rnorm(2000),400,5)
#' beta1 = as.matrix(c(1,1,0,0,0))
#' beta2 = as.matrix(c(0,0,1,1,0))
#' err = as.matrix(rnorm(400))
#' y = (x%*%beta1)^2+x%*%beta2+err
#' dr = mave(y~x, method = 'meanopg')
#' dr.dim = mave.dim(dr)
#' plot(dr,dim=3)
#' plot(dr.dim)
#' @export
#' @method plot mave



plot.mave<-function(x,dim=4,plot.method=pairs,...){
  dr <- x
  dir <- mave.dir(dr,dim)
  if(dr$call[[1]]=='mave'){
    y.name=all.vars(dr$call)[1]
  }
  else{
    y.name=all.vars(dr$call)[2]
  }

  plot.method(cbind(dr$y,dr$x%*%dir),label=c(y.name,colnames(dir)))
}

#' @rdname plot.mave
#' @export
plot.mave.dim<-function(x,dim='dim.min',plot.method=pairs,...){
  dr.dim <- x
  if(dim=='dim.min'){
    dir <- mave.dir(dr.dim)
  }else{
    dir <- mave.dir(dr.dim,dim)
  }

  if(dr.dim$call[[1]]=='mave'){
    y.name=all.vars(dr.dim$call)[1]
  }
  else{
    y.name=all.vars(dr.dim$call)[2]
  }
  plot.method(cbind(dr.dim$y,dr.dim$x%*%dir),label=c(y.name,colnames(dir)))
}
