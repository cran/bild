gradlogL.bin0<- function(parameters, X,data, trace)
{
  gradient1<- function(param, X, y)
  {
    npar <- as.integer(length(param)+1)
    beta <- as.double(param[1:(npar-1)])
    log.psi <- as.double(0)
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    theta <- work <- as.double(rep(0,n))
    gbeta <- dbeta <- dbeta1<- as.double(rep(0,npar-1))
    glpsi<- as.double(0)
    db <- matrix(as.double(0),nrow=3,ncol=npar-1)
    x <- matrix(as.double(X),nrow=n, ncol=npar-1)
    der<-as.double(rep(0,npar-1))
    
    results <- .Fortran("mbgd1",gbeta,glpsi,beta,log.psi, 
                        npar,x,y,theta,work,der,db,dbeta,dbeta1,n,PACKAGE="bild")
    
    
    return(results[[1]])}
  
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  counts<-data[[3]]
  gradlogL<-0
  k1<-1
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<- gradient1(param=parameters, X=X[k1:k2,],y=y[k1:k2])
    gradlogL<-gradlogL+counts[i]*z
    k1<-k2+1
  }
  return(-gradlogL)}	
