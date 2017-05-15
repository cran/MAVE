### R code from vignette source 'Short_introduction_to_package_MAVE.Rnw'

###################################################
### code chunk number 1: Short_introduction_to_package_MAVE.Rnw:39-47
###################################################
set.seed(12345)
library(MAVE)
x <- matrix(rnorm(400*5),400,5)
b1 <- matrix(c(1,1,0,0,0),5,1)
b2 <- matrix(c(0,0,1,1,0),5,1)
eps <- matrix(rnorm(400),400,1)
y <- x%*%b1 + (x%*%b2)*eps
dr <- mave(y~x, method='csmave')


###################################################
### code chunk number 2: Short_introduction_to_package_MAVE.Rnw:52-53
###################################################
dr


###################################################
### code chunk number 3: Short_introduction_to_package_MAVE.Rnw:56-58
###################################################
dr.dim <- mave.dim(dr)
dr.dim


###################################################
### code chunk number 4: Short_introduction_to_package_MAVE.Rnw:61-62
###################################################
names(dr.dim)


###################################################
### code chunk number 5: Short_introduction_to_package_MAVE.Rnw:65-66
###################################################
dr.dim$dir


###################################################
### code chunk number 6: Short_introduction_to_package_MAVE.Rnw:69-70
###################################################
mave.dir(dr.dim)


