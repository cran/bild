########to compute pij

logL.bin2I.aux<- function(parameters, X, data,trace)
{
  loglik2 <- function(param, X, y,trace)
  {
    npar <-as.integer(length(param)-1)
    beta<- as.double(param[1:(npar-2)])
    bt<- as.double(param[1:(npar-2)])
    lpsi<-as.double(param[(npar-1):npar])
    omega<-as.double(param[npar+1])
    y<- as.integer(y)
    n <- as.integer(length(y)) 
    x<-matrix(as.double(X),nrow=n,ncol=npar-2)
    theta<- work<- pij<-prob<-as.double(rep(0,n))
    logL <- as.double(0)
    eta<-i.fit<-fit<-as.vector(npar-2)
    eta<-x%*%beta
    
##### Estimate bi 
    m<-glm(as.numeric(y)~offset(eta), family=binomial)

    # bi<-coef(m) #27/03/2020
    bi<-coef(m)[1]
    beta[1]<-beta[1]+bi
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    
    results <- .Fortran("blik2m",logL,pij,beta,lpsi, 
                        npar,x,y,theta,work,n,PACKAGE="bild")
    
    i.fit<-x%*%beta

    return(list(loglik=results[[1]],pij=results[[2]],fit=i.fit, bi.est=bi ))}
  
  if(trace)	cat(paste(format(parameters[length(parameters)-2], digit=4), collapse=" "), "\t")
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)-1], digit=4)), collapse=" "), "\t")
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")
  
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  prob<-fitted<-as.double(rep(0,length(y)))
  bi.estimate<-as.double(rep(0,length(n.cases)))
  counts<-data[[3]]
  logL<-as.double(0)
  k1<-1
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<- loglik2(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
    prob[k1:k2]<-z$pij
    fitted[k1:k2]<-z$fit
    bi.estimate[i]<-z$bi.est
    k1<-k2+1
  }
  
  for (j in 1:length(y))
    {
    if (prob[j] > 0.999999999999) prob[j]<- 0.999999999999
#    else if  (prob[j] < 1e-11) prob[j]<- 1e-11
    }

  
  return(list(pij=prob, fit=fitted,bi.est=bi.estimate))}
