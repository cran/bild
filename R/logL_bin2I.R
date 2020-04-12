logL.bin2I<- function(parameters, X, data,integrate,trace)
{
  loglik2 <- function(param, X, y,integrate,trace)
  {
    npar <-as.integer(length(param)-1)
    beta<- as.double(param[1:(npar-2)])
    bt<- as.double(param[1:(npar-2)])
    log.psi<-as.double(param[(npar-1):npar])
    omega<-as.double(param[npar+1])
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    n <- as.integer(length(y)) 
    x<-matrix(as.double(X),nrow=n,ncol=npar-2)
    theta<- work<- as.double(rep(0,n))
    logL <- as.double(0)
    
    li<-as.double(integrate$li)
    ls<-as.double(integrate$ls)
    epsabs<-as.double(integrate$epsabs)
    epsrel<-as.double(integrate$epsrel)
    limit<-as.integer(integrate$limit)
    key<-as.integer(integrate$key)
    
    results <- .Fortran("integ",logL,bt,beta,log.psi,omega,npar,
                        x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="bild")
    
    
    return(results[[1]])}
  
  if(trace)	cat(paste(format(parameters[length(parameters)-2], digit=4), collapse=" "), "\t")
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)-1], digit=4)), collapse=" "), "\t")
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")
  
  nparam<-length(parameters)
  omega1<-as.double(parameters[nparam])
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  counts<-data[[3]]
  mt<-as.double(rep(0,length(y)))
  logL1<-as.double(0)
  k1<-1
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<-loglik2(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate,trace=trace)
    #logL1 gives the log-likelihood
    logL1<-logL1+counts[i]*z
    k1<-k2+1
  }
  if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
  
  return(nloglik=-logL1)}
