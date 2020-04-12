setMethod("summary",
    signature(object = "bild"),
    function(object,cov=FALSE,cor=FALSE)
	{
	cat("\nCall:\n")
	print(object@call)
	cat("\nNumber of profiles in the dataset: ", object@ni.cases)
	cat("\nNumber of profiles used in the fit: ", object@n.cases)
	cat("\nDependence structure:", object@call$dependence)
	cat("\nLog likelihood: ", round(object@log.likelihood, 4))
	cat("\nAIC: ", round(object@aic, 4),"\n")
  
	coef <- object@coefficients
	nas <- is.na(coef[, 1])
	cnames <- names(coef[, 1][!nas])
	coef <- matrix(rep(coef[, 1][!nas], 4), ncol = 4)
	coef.aux<-matrix(rep(coef[, 1][!nas], 4), ncol = 4)
	#	coef[, 1] <- 1:dim(coef)[[1]]
	coef[, 2] <- object@se[, 1][!nas]
	coef[, 3] <- round(coef[, 1]/coef[, 2], 3)
	coef.aux[,1]<-pnorm(coef[, 1]/coef[, 2])
	coef.aux[,2]<-1-pnorm(coef[, 1]/coef[, 2])
	for (i in 1:dim(coef)[[1]])
	{coef.aux[i,3]<-min(coef.aux[i, 1],coef.aux[i, 2])}
	coef[, 4] <- round(2*coef.aux[,3],6)
	dimnames(coef) <- list(cnames, c("Estimate", "Std. Error", "z value", "p-value"))
	
	# for dependence="ind" 
	if(all(is.na(match(cnames, "log.psi1"))) && all(is.na(match(cnames, "omega"))))
	   
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[ ,  ])
	}
	
	# for dependence="MC1" 
	else if(!all(is.na(match(cnames, "log.psi1"))) && all(is.na(match(cnames, "log.psi2")))  &&
	        all(is.na(match(cnames, "omega"))))
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:(dim(coef)[[1]]-1),  ])
	  
	  cat("\nDependence parameter: \t\n") 
	  table.r0<-matrix(coef[dim(coef)[[1]], ], nrow=1)
	  dimnames(table.r0)<-list("log.psi1", c("Estimate", "Std. Error", "z value", "p-value"))
	  print(table.r0)
	}
	
	# for dependence="MC2" 
	else if(!all(is.na(match(cnames, "log.psi1"))) && !all(is.na(match(cnames, "log.psi2")))  &&
	        all(is.na(match(cnames, "omega"))))
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:(dim(coef)[[1]]-2) ,  ])
	  
	  cat("\nDependence parameter: \t\n") 
	  table.r2<-data.frame(coef[(dim(coef)[[1]]-1):dim(coef)[[1]], ])
	  dimnames(table.r2)<-list(c("log.psi1","log.psi2"),  c("Estimate", "Std. Error", "z value", "p-value"))
	  print(table.r2)
	}
	
	
	# for dependence="indR" 
	else if (all(is.na(match(cnames, "log.psi1"))) && !all(is.na(match(cnames, "omega"))))
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[-dim(coef)[[1]],  ])
	  
	  coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
	  cat("\nRandom effects:\t\n") 
	  table.r1<-data.frame(coef[dim(coef)[[1]], 1])
	  dimnames(table.r1)<-list("(Intercept)", "Variance")
	  print(table.r1)
	}
	
	# for dependence="MC1R" 
	else if (!all(is.na(match(cnames, "log.psi1"))) && all(is.na(match(cnames, "log.psi2")))  &&
	        !all(is.na(match(cnames, "omega"))))
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:(dim(coef)[[1]]-2),  ])
	  
	  cat("\nDependence parameter: \t\n") 
	  table.r0<-matrix(coef[dim(coef)[[1]]-1, ], nrow=1)
	  dimnames(table.r0)<-list("log.psi1", c("Estimate", "Std. Error", "z value", "p-value"))
	  print(table.r0)
	  
	  ### for random coef  
	  coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
	  cat("\nRandom effect:\t\n") 
	  table.r1<-data.frame(coef[dim(coef)[[1]], 1])
	  dimnames(table.r1)<-list("(Intercept)", "Variance")
	  print(table.r1)
	}
	
	# for dependence="MC2R" 
	else if (!all(is.na(match(cnames, "log.psi1"))) && !all(is.na(match(cnames, "log.psi2")))  &&
	         !all(is.na(match(cnames, "omega"))))
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:(dim(coef)[[1]]-3),  ])
	  
	  cat("\nDependence parameter: \t\n") 
	  table.r2<-data.frame(coef[(dim(coef)[[1]]-2):(dim(coef)[[1]]-1), ])
	  dimnames(table.r2)<-list(c("log.psi1","log.psi2"),  
	                           c("Estimate", "Std. Error", "z value", "p-value"))
	  print(table.r2)
	  
	  ### for random coef  
	  coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
	  cat("\nRandom effect:\t\n") 
	  table.r1<-data.frame(coef[dim(coef)[[1]], 1])
	  dimnames(table.r1)<-list("(Intercept)", "Variance")
	  print(table.r1)
	}
	
	if (cov)
	{
	  cat("\nCovariance of Coefficients: \n")
	  print(object@covariance, digits = 2)}
	if (cor){
	  cat("\nCorrelation of Coefficients: \n")
	  print(object@correlation, digits = 2)}
	cat("\nMessage: ", object@message,"\n")
    })


