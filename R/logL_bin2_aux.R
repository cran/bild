logL.bin2.aux <- function(parameters, X, data, trace)
{
  loglik2 <- function(param, X, y,trace)
  {
    npar <-as.integer(length(param))
    beta <- as.double(param[1:(npar-2)])
    lpsi <- as.double(param[(npar-1):npar])
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y)) 
    x <-matrix(as.double(X),nrow=n,ncol=npar-2)
    theta <- work <- pij<-as.double(rep(0,n))
    logL <- L<-as.double(0)
    
    results <- .Fortran("blik2m",logL,pij,beta,lpsi,npar,
                        x,y,theta,work,n,PACKAGE="bild")
    
    return(list(loglik=results[[1]],pij=results[[2]]))}
  
  if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=3), collapse=" "), "\t")
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=3)), collapse=" "), "\t")
  
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
    z<- loglik2(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
    logL<-logL+counts[i]*z$loglik
    prob[k1:k2]<-z$pij
    k1<-k2+1
  }
  if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
  
  return(list(nloglik=-logL, pij=prob))}
