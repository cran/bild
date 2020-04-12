gradlogL.bin1I <- function(parameters, X,data,integrate,trace)
{
  gradient1 <-  function(param,X,y,integrate)
  {
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    npar <- as.integer(length(param)-1)
    beta <- as.double(param[1:(npar-1)])
    bt <- as.double(param[1:(npar-1)])
    lpsi <- as.double(param[npar])
    omega<-as.double(param[npar+1])
    theta <- work <- as.double(rep(0,n))
    g.beta <- d.beta <- d.beta1<- der <- as.double(rep(0,npar-1))
    g.lpsi1<-g.lpsi2 <- as.double(0)
    d.beta <- matrix(as.double(0),3,npar-1)
    gvar<-as.double(0)
    x <- matrix(as.double(X),nrow=n, ncol=npar-1)
    db <- matrix(as.double(0),3,npar-1)
    
    li<-as.double(integrate$lig)
    ls<-as.double(integrate$lsg)
    epsabs<-as.double(integrate$epsabs)
    epsrel<-as.double(integrate$epsrel)
    limit<-as.integer(integrate$limit)
    key<-as.integer(integrate$key)
    
    
    result <- .Fortran("gint1",g.beta,g.lpsi1,gvar,bt,beta,lpsi,omega,npar,x,y,theta,work,n,
                       d.beta,d.beta1,der,db,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")
    
    
    return(c(result[[1]],result[[2]],result[[3]]))}
  
  loglik1 <- function(param, X, y,integrate)
  {
    npar <- as.integer(length(param)-1)
    beta <- as.double(param[1:(npar-1)])
    bt <- as.double(param[1:(npar-1)])
    lpsi <- as.double(param[npar])
    omega<-as.double(param[npar+1])
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    x <-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta <- work <-as.double(rep(0,n))
    logL <- as.double(0)
    
    li<-as.double(integrate$lig)
    ls<-as.double(integrate$lsg)
    epsabs<-as.double(integrate$epsabs)
    epsrel<-as.double(integrate$epsrel)
    limit<-as.integer(integrate$limit)
    key<-as.integer(integrate$key)
    
    results <- .Fortran("integ1",logL,bt,beta,lpsi,omega,npar,
                        x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")
    
    return(results[[1]])}
  
  nparam <- as.integer(length(parameters)-1)
  omega1<-parameters[nparam+1]
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  counts<-data[[3]]
  dbeta<-as.double(rep(0,nparam-1))
  dlog.psi1<-0
  dvar<-0
  k1<-1
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
    grad<-gradient1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
    
    for (j in 1:(nparam-1))
    {
      dbeta[j]<-dbeta[j]+counts[i]*(grad[j]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))	
    }
    dlog.psi1<-dlog.psi1+counts[i]*(grad[nparam]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))		
    #using the chain rule
    dvar<-dvar+counts[i]*(grad[nparam+1]/(exp(z)*sqrt(2*pi)*exp(omega1/2)))*exp(omega1)
    k1<-k2+1
  }
  gr<-c(dbeta,dlog.psi1,dvar)
  return(-gr)}
