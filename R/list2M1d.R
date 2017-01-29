list2M1d<-function(BB){
  p=length(BB)
  B=rep(0,p^3)
  for(i in 1:p){
    B[((i-1)*p^2+1):((i-1)*p^2+i*p)]=BB[[i]]
  }
  return(B)
}
