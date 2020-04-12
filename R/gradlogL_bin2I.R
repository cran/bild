
gradlogL.bin2I <- function(parameters, X,data,integrate,trace)
{
  gradient2 <-  function(param,X,y,integrate)
  {
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    npar <- as.integer(length(param)-1)
    beta <- as.double(param[1:(npar-2)])
    bt <- as.double(param[1:(npar-2)])
    lpsi <- as.double(param[(npar-1):npar])
    omega<-as.double(param[npar+1])
    theta <- work <- as.double(rep(0,n))
    g.beta <- d.beta <- d.beta1<-d.beta2 <- der <- as.double(rep(0,npar-2))
    g.lpsi1<-g.lpsi2 <- as.double(0)
    d.beta <- matrix(as.double(0),3,npar-2)
    gvar<-as.double(0)
    x <- matrix(as.double(X),nrow=n, ncol=npar-2)
    der <- as.double(rep(0,npar-2))
    db <- matrix(as.double(0),3,npar-2)
    db1<- matrix(as.double(0),4,npar-2)
    db2<- matrix(as.double(0),5,npar-2)
    
    li<-as.double(integrate$lig)
    ls<-as.double(integrate$lsg)
    epsabs<-as.double(integrate$epsabs)
    epsrel<-as.double(integrate$epsrel)
    limit<-as.integer(integrate$limit)
    key<-as.integer(integrate$key)
    
    result <- .Fortran("gint",g.beta,g.lpsi1,g.lpsi2,gvar,
                       x,theta,work,y,lpsi,beta,bt,d.beta,d.beta1,d.beta2,
                       der,db,db1,db2,omega, npar,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")
    
    return(c(result[[1]],result[[2]],result[[3]],result[[4]]))}
  
  loglik2 <- function(param, X, y,integrate)
  {
    npar <- as.integer(length(param)-1)
    beta <- as.double(param[1:(npar-2)])
    bt <- as.double(param[1:(npar-2)])
    lpsi <- as.double(param[(npar-1):npar])
    omega<-as.double(param[npar+1])
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    x <-matrix(as.double(X),nrow=n,ncol=npar-2)
    theta <- work <-as.double(rep(0,n))
    logL <- as.double(0)
    
    li<-as.double(integrate$lig)
    ls<-as.double(integrate$lsg)
    epsabs<-as.double(integrate$epsabs)
    epsrel<-as.double(integrate$epsrel)
    limit<-as.integer(integrate$limit)
    key<-as.integer(integrate$key)
    
    results <- .Fortran("integ",logL,bt,beta,lpsi,omega,npar,
                        x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")
    
    
    return(results[[1]])}
  
  nparam <- as.integer(length(parameters)-1)
  omega1<-parameters[nparam+1]
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  counts<-data[[3]]
  dbeta<-as.double(rep(0,nparam-2))
  dlog.psi1<-0
  dlog.psi2<-0
  dvar<-0
  k1<-1
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<-loglik2(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
    grad<-gradient2(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
    
    for (j in 1:(nparam-2))
    {
      dbeta[j]<-dbeta[j]+counts[i]*(grad[j]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
    }
    dlog.psi1<-dlog.psi1+counts[i]*(grad[nparam-1]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
    dlog.psi2<-dlog.psi2+counts[i]*(grad[nparam]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
    #using the chain rule
    dvar<-dvar+counts[i]*(grad[nparam+1]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))*exp(omega1)
    k1<-k2+1
  }
  gr<-c(dbeta,dlog.psi1,dlog.psi2,dvar)
  return(-gr)}
