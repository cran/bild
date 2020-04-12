gradlogL.bin2<- function(parameters, X,data, trace)
{
  gradient2 <-  function(param,X,y)
  {
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    npar <- as.integer(length(param))
    beta <- as.double(param[1:(npar-2)])
    lpsi <- as.double(param[(npar-1):npar])
    theta <- work <- as.double(rep(0,n))
    g.beta <- d.beta <- d.beta1<-d.beta2 <- der <- as.double(rep(0,npar-2))
    g.lpsi1<-g.lpsi2 <- as.double(0)
    db <- matrix(as.double(0),3,npar-2)
    db1 <- matrix(as.double(0),4,npar-2)
    db2 <- matrix(as.double(0),5,npar-2)
    x <- matrix(as.double(X),nrow=n, ncol=npar-2)
    
    result <- .Fortran("bgd2m",g.beta,g.lpsi1,g.lpsi2,beta,
                       lpsi,npar,x,y,theta,work,d.beta,d.beta1,d.beta2,n,der,db,db1,db2,PACKAGE="bild")
    
    return(c(result[[1]],result[[2]],result[[3]]))
  }
  nparam <- as.integer(length(parameters))
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
    z<- gradient2(param=parameters, X=X[k1:k2,],y=y[k1:k2])
    gradlogL<-gradlogL+counts[i]*z
    k1<-k2+1
  }
  return(-gradlogL)}