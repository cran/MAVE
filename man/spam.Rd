% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spam.R
\docType{data}
\name{spam}
\alias{spam}
\title{4601 email record}
\format{
A data frame with 4601 rows and 57 variables
}
\description{
A dataset containing 4601 record of email with 57 features. These features are the
relative frequency of most commonly used phrases and punctions. The data of these
features are recorded 1 to 57 columns of the spam data. The outcome is spam or email which is
denoted as 1 or 0, recorded in the 58th column of the data.
}
\examples{
data(spam)
train = sample(1:4601)[1:1000]
x.train <- as.matrix(spam[train,1:57])
y.train <- as.matrix(spam[train,58])
x.test <- as.matrix(spam[-train,1:57])
y.test <- as.matrix(spam[-train,58])
x.train <- sqrt(x.train)
x.test <- sqrt(x.test)
}
\keyword{datasets}
