M3d2list<-function(BB){
	p=dim(BB)[1]
	B=list()
	B[[1]]=matrix(BB[,1,1])
	for(i in 2:p){
		B[[i]] = BB[,1:i,i]
	}
	return(B)
}

