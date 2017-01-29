### R code from vignette source 'Short_introduction_to_package_MAVE.Rnw'

###################################################
### code chunk number 1: Short_introduction_to_package_MAVE.Rnw:37-42
###################################################
X <- matrix(rnorm(400),100,4)
eps <- matrix(rnorm(100),100,1)
Y <- as.matrix(X[,1]+X[,2]+(X[,3]+X[,4])*eps)
library(MAVE)
rd <- MAVE(X,Y)


###################################################
### code chunk number 2: Short_introduction_to_package_MAVE.Rnw:47-48
###################################################
rd


###################################################
### code chunk number 3: Short_introduction_to_package_MAVE.Rnw:52-54
###################################################
rd <- DIM(rd)
rd


###################################################
### code chunk number 4: Short_introduction_to_package_MAVE.Rnw:57-58
###################################################
names(rd)


###################################################
### code chunk number 5: Short_introduction_to_package_MAVE.Rnw:61-62
###################################################
rd$dir


