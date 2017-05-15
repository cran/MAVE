M3d2list<-function(BB,xnames){
	p=dim(BB)[1]
	B=list()
	for(i in 1:p){
		B[[i]] = as.matrix(BB[,1:i,i])
		colnames(B[[i]])<-paste('dir',1:i,sep='')
    rownames(B[[i]])<-xnames
	}
	return(B)
}

