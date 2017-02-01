#'@method print mave
#'@export
print.mave<-function(x,...){
  cat('\nCall:\n',paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
  arg='central space'
  if(substr(x$method,1,2)!='CS') arg='central mean space'
  p=length(x$dir)
  if(sum('cv'==names(x))==0){
    if(p<4){
      print(x$dir)
    }else{
      print(x$dir[[1]])
      print(x$dir[[2]])
      print(x$dir[[3]])
    }
    cat('(only the first 3 ',arg,' are displayed,
        to display the space of dimension k, call object$dir[[k]])\n')
  }
  else{
    pp=1
    while(pp<=p){
      np=min(pp+10,p)
      cat('Dimension\t')
      for(i in pp:np) cat(i,'\t')
      cat('\n')
      cat('CV-value\t')
      for(i in pp:np) cat(round(x$cv[i],2),'\t')
      cat('\n')
      pp=np+1
    }
    cat('\n')
    d=which(x$cv==min(x$cv))
    cat('The selected dimension is ',d)
    cat('\n\n')
    #cat(paste('The matrix of the best',arg,'selected by cross-validation is of',d,'dimensions,','\nwhich is given below\n',sep=' '))
    #print(x$dir[[d]])
  }
}
