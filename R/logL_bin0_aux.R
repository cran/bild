
logL.bin0.aux<- function(parameters, X, data, trace)
{
  loglik1<- function(param, X, y, trace)
  {#calculate logLi(beta/bi) for each individual
    npar <- as.integer(length(param)+1)
    beta <- as.double(param[1:(npar-1)])
    log.psi <- as.double(0)
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    x <-matrix(as.double(X),nrow=n,ncol=(npar-1))
    theta <- work <- pij<-as.double(rep(0,n))
    logL <- L<-as.double(0)
    
    results <- .Fortran("mblik1",logL,pij,beta,log.psi, 
                        npar,x,y,theta,work,n,PACKAGE="bild")
    
    return(list(loglik=results[[1]],pij=results[[2]]))
  }
  
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  prob<-as.double(rep(0,length(y)))
  counts<-data[[3]]
  logL<-as.double(0)
  k1<-1
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
    logL<-logL+counts[i]*z$loglik
    prob[k1:k2]<-z$pij
    k1<-k2+1
  }
  if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
  
  return(list(nloglik=-logL, pij=prob))}
