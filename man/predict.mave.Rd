% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_mave.R
\name{predict.mave}
\alias{predict.mave}
\alias{predict.mave.dim}
\title{Make predictions based on the dimension reduction space}
\usage{
\method{predict}{mave}(object, newx, dim, ...)

\method{predict}{mave.dim}(object, newx, dim = "dim.min", ...)
}
\arguments{
\item{object}{the object of class 'mave'}

\item{newx}{Matrix of the new data to be predicted}

\item{dim}{the dimension of central space or central mean space. The matrix of the original data will be
multiplied by the matrix of dimension reduction directions of given dimension. Then the prediction will be
made based on the data of given dimensions. The value of dim should be given when the class of the
argument dr is mave. When the class of the argument dr is mave.dim and dim is not given, the
function will return the basis matrix of CS or CMS of dimension selected by \code{\link{mave.dim}}}

\item{...}{further arguments passed to \code{\link{mars}} function such as degree.}
}
\value{
the prediced response of the new data
}
\description{
This method make predictions based the reduced dimension of data using \code{\link{mars}} function.
}
\examples{

X = matrix(rnorm(10000),1000,10)
beta1 = as.matrix(c(1,1,1,1,0,0,0,0,0,0))
beta2 = as.matrix(c(0,0,0,1,1,1,1,1,0,0))
err = as.matrix(rnorm(1000))
Y = X\%*\%beta1+X\%*\%beta2+err

train = sample(1:1000)[1:500]
x.train = X[train,]
y.train = as.matrix(Y[train])
x.test = X[-train,]
y.test = as.matrix(Y[-train])

dr = mave(y.train~x.train, method = 'meanopg')

yp = predict(dr,x.test,dim=3,degree=2)
#mean error
mean((yp-y.test)^2)

dr.dim = mave.dim(dr)

yp = predict(dr.dim,x.test,degree=2)
#mean error
mean((yp-y.test)^2)

}
\seealso{
\code{\link{mave}} for computing the dimension reduction space and \code{\link{mave.dim}} for
estimating the dimension of the dimension reduction space
}
