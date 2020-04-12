
setClass("summary.bild", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",call="language"))

setClass("bild", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric", residuals="numeric", 
		s.residuals="numeric",ind.probability="numeric", prob.matrix="matrix", Fitted="numeric", bi.estimate="matrix",
		Fitted.av="numeric", Time="numeric", model.matrix= "matrix", y.matrix="matrix",
		subset.data="data.frame", y.av="numeric", f.value="factor",call="language"))


setGeneric("getAIC",def=function(object) standardGeneric("getAIC"))

setGeneric("getLogLik",def=function(object) standardGeneric("getLogLik"))

setGeneric("getcoef",def=function(object) standardGeneric("getcoef"))

setGeneric("getvcov",def=function(object) standardGeneric("getvcov"))

setGeneric("randeff",def=function(object) standardGeneric("randeff"))

setGeneric("fixeff",def=function(object) standardGeneric("fixeff"))

setGeneric("vareff",def=function(object) standardGeneric("vareff"))

setGeneric("model.mat",def=function(object) standardGeneric("model.mat"))


bild<-function(formula = formula(data), data, time,id, subSET, aggregate=FALSE, start = NULL, trace = FALSE,
dependence="ind", method="BFGS", control=bildControl(), 
integrate=bildIntegrate())
{
# *****************DEFINITION OF INTERNAL FUNCTIONS ******************
# na.action for binary families

na.discrete.replace1 <- function(frame,  n.times, ti.repl)
	{
	vars <- names(frame)
	names(vars) <- vars
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	badlines<-NULL
	bad.ind<-NULL
	for(j in 1:length(vars)) 
	{k1<-1
	for (i in 1:n.cases)
	{k2<-cumti.repl[i]
	x <- frame[[j]][k1:k2]
	pos <- is.na(x)
	if(any(pos))
			if(j == 1){distance.between.na <- diff(seq(1, n.times)[!pos])
					if (any(distance.between.na > 2 ))
					{badlines<-c(badlines,c(k1:k2))
					bad.ind<-c(bad.ind,i)
					x[pos] <- -1}
				x[pos] <- -1}
			else stop("NA's on covariates not allowed")
	frame[[j]][k1:k2]<-x
	k1<-k2+1
	} 
	}
		if (length(bad.ind)>=1)
		cat("Warning Message: Condition on NA's not respected:\nresults might be inaccurate\n")
		return(list(data=frame, badlines=badlines,bad.ind=bad.ind))
	}


na.discrete.replace2 <- function(frame,  n.times, ti.repl)
	{
	vars <- names(frame)
	names(vars) <- vars
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	badlines<-NULL
	bad.ind<-NULL
	for(j in 1:length(vars)) 
	{k1<-1
	for (i in 1:n.cases)
	{k2<-cumti.repl[i]
	x <- frame[[j]][k1:k2]
	pos <- is.na(x)
	if(any(pos))
			if(j == 1){distance.between.na <- diff(seq(1, n.times)[!pos])
				if (any(distance.between.na > 2 ))
				{badlines<-c(badlines,c(k1:k2))
				bad.ind<-c(bad.ind,i)
				x[pos] <- -1}
				if (any(distance.between.na == 2 ) & length(distance.between.na)>1)
				{for (i2 in 1:(length(distance.between.na)-1))
				{if (distance.between.na[i2]==2 & distance.between.na[i2+1]==2)
				{badlines<-c(badlines,c(k1:k2))
				bad.ind<-c(bad.ind,i)
				x[pos] <- -1}}}
			x[pos] <- -1}
			else stop("NA's on covariates not allowed")
	frame[[j]][k1:k2]<-x
	k1<-k2+1
	} 
	}
		if (length(bad.ind)>=1)
		cat("Warning Message: Condition on NA's not respected:\nresults might be inaccurate\n")
		return(list(data=frame, badlines=badlines,bad.ind=bad.ind))
	}

#
# compute the fitted values for variuos families of distributions
#
	various.fitted <- function(Fitted)
	{
		#return(1/(1 + care.exp( - Fitted)))
		return(1/(1 +( - Fitted)))
	}
#

# ******************* MAIN PROGRAM *******************************
#
	call <- match.call()
#	vect.time <- F
	if(missing(data) || !is.data.frame(data))
		stop("a data.frame must be supplied")
	if(is.null(names(data)))
		stop("objects in data.frame must have a name")
	expr1 <- terms(formula, data=data)
	expr <- attr(expr1, "variables")
	var.names <- all.vars(expr)
	response <- all.vars(expr)[1]
	
if(any(is.na(match(var.names, names(data)))))
		stop("Variables in formula not contained in the data.frame")

if(!missing(time)) {Time<-as.vector(data[[time]])}
if (missing(time)) { if (all (is.na(match(names(data), "time")))) stop ("time must be defined")
			else  Time<-as.vector(data$time)}
				
if(!missing(id)) {id<-as.vector(data[[id]])}
if (missing(id)){ if (all(is.na(match(names(data), "id")))) stop ("id must be defined")
		else id<-as.vector(data$id)}

# select subset if necessary
	if(!missing(subSET)) {id1 <- eval(substitute(subSET), data)
			      data<-subset(data, id1)}

	if(!missing(aggregate)) {f.name <- deparse(substitute(aggregate))
  			f.value <- as.factor(data[[f.name]])}
	

	if(missing(aggregate)) {f.value<-as.factor(0)}

#returns data of a subset
subset.data<-data

	ti.repl<-as.vector(0)
	i1<-1
	i2<-1
	for (i in 1:(length(data[[response]])-1))
	{ 
		if (id[i]==id[i+1])
		{ 	i2<-i2+1
			ti.repl[i1]<-i2}
		else {	ti.repl[i1]<-i2
			i1<-i1+1
			i2<-1}
	}

	n.cases <- length(ti.repl)
	n.tot<-cumsum(ti.repl)[n.cases]
	n.time<-length(unique(Time))
	ni.cases <- length(ti.repl)
	pos.ind<-cumsum(ti.repl)

	counts<-as.vector(0)

	if(is.na(match("counts", names(data))))
			counts <- data$counts <- rep(1, n.cases)
		else   {for (i in 1:n.cases)
			{counts[i]<-data$counts[pos.ind[i]]}
		}

	final.data <- data
	var.names <- c(var.names, "counts")
	data <- data[var.names]
	n.var <- length(data)
	Y.resp <- as.vector(data[[response]])
	if((all(Y.resp[!is.na(Y.resp)] == 1 | Y.resp[!is.na(Y.resp)] == 0)) == FALSE)
		stop("Unfeasible values of response variable: must be 0,1,NA")

# ********** creation of individual profile according to NA patterns *******************

	data2<-data
		
	if (dependence=="MC2"|| dependence=="MC2R")

	{final.data <- na.discrete.replace2(frame=data,  n.times=n.time, ti.repl=ti.repl)}

	else if (dependence=="ind"|| dependence=="indR"||dependence=="MC1"|| dependence=="MC1R")

	{final.data <- na.discrete.replace1(frame=data,  n.times=n.time, ti.repl=ti.repl)}


if (length(final.data$bad.ind)>=1)
	{data<-final.data$data
	data<-data[-final.data$badlines,]
	data2<-data2[-final.data$badlines,]
	counts<-counts[-final.data$bad.ind]
	ti.repl<-ti.repl[-final.data$bad.ind]
	n.cases <- length(ti.repl)
	n.tot<-cumsum(ti.repl)[n.cases]}
else {data<-final.data$data}


# ********** design matrices creation *******************
	# define a plausible starting point for the optimizer if not given
		data1 <- na.omit(data2)
		data1.resp <- data1[, response]
		data1.resp <- log((data1.resp + 0.5)/(1.5 - data1.resp))
		data1[, c(response)] <- data1.resp

		if (dependence=="MC1")	 init<-0
		else  if (dependence=="MC1R")  init<-c(0,0)
		else  if (dependence=="MC2")    init<-c(0,0)
		else  if (dependence=="MC2R")   init<-c(0,0,0)
		else  if (dependence=="indR")  init<-0

		if(is.null(start) && dependence!="ind")
			start <- c(lm(formula, data1, weights = counts)$coefficients, init)
		else if(!is.null(start) && dependence!="ind")
			start <- c(lm(formula, data1, weights = counts)$coefficients, start)
		else if (dependence=="ind") start <- c(lm(formula, data1, weights = counts)$coefficients)

if (any(is.na(start))) stop("starting values produced by lm contains NA")

	id.not.na<-rep(TRUE,n.tot)
	X <- model.matrix(expr1, data)
	names.output <- dimnames(X)[[2]]
	sum.ti <- sum(ti.repl)
	data <- list(ti.repl, data[[response]], counts)
	data2<-list(ti.repl, data2[[response]], counts)
	p <- dim(X)[2] + 1
	pr<-prob<-aux<-as.double(rep(0,length(data[[2]])))

	if (dependence=="ind")
	{	if(trace)	cat("\t log.likelihood\n")
	temp <-optim(par= start, fn = logL.bin0,  gr= gradlogL.bin0, method=method, 
	data = data, X = X, trace=trace,control=control)}

	else  if (dependence=="indR")
	{	if(trace)	cat("\n omega\t log.likelihood\n")
	temp <- optim(par = start, fn =logL.bin0I,gr = gradlogL.bin0I,  method=method, 
	data = data, X = X, integrate=integrate, trace=trace,control=control)}

	else if (dependence=="MC1")
	{	if(trace)	cat("\n log.psi1\t log.likelihood\n")
	temp <-optim(par= start, fn = logL.bin1 ,  gr= gradlogL.bin1, method=method, 
	data = data, X = X, trace=trace,control=control)}
	
	else  if (dependence=="MC1R")
	{	if(trace)	cat("\n log.psi1\t omega\t log.likelihood\n")
	temp <- optim(par = start, fn =logL.bin1I,gr = gradlogL.bin1I,  method=method, 
	data = data, X = X, integrate=integrate, trace=trace,control=control)}
	
	else  if (dependence=="MC2")
	{	if(trace)	cat("\n log.psi1\t log.psi2\t log.likelihood\n")
	temp <- optim(par= start, fn =logL.bin2, gr =gradlogL.bin2, method=method,
	 data = data, X = X, trace = trace,control=control)}
	
	else  if (dependence=="MC2R")
	{	if(trace)	cat("\n log.psi1\t log.psi2\t omega\t log.likelihood\n")
	temp <- optim(par = start, fn= logL.bin2I, gr = gradlogL.bin2I,method=method, 
	data = data, X = X,integrate=integrate, trace = trace,control=control)}

	coefficients <- temp$par
	log.lik <-  - temp$value
	if (trace) 
	cat("Convergence reached. Computing the information matrix now\n")

	 if (dependence=="ind")
	Info <- num.info(coefficients, "gradlogL.bin0", X, data)
	else  if (dependence=="indR")
	Info <- num.infoI(coefficients, "gradlogL.bin0I", X, data, integrate=integrate)
	else if (dependence=="MC1")
	Info <- num.info(coefficients, "gradlogL.bin1", X, data)
	else  if (dependence=="MC1R")
	Info <- num.infoI(coefficients, "gradlogL.bin1I", X, data, integrate=integrate)
	else  if (dependence=="MC2")
	Info <- num.info(coefficients, "gradlogL.bin2", X, data)
	else  if (dependence=="MC2R")
	Info <- num.infoI(coefficients, "gradlogL.bin2I", X, data, integrate=integrate)

	se <- matrix(sqrt(diag(solve(Info))), ncol = 1)
	coefficients <- matrix(coefficients, ncol = 1)
	if (dependence=="ind")
	dimnames(coefficients) <- dimnames(se) <- list(names.output, " ")
	else  if (dependence=="indR")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "omega"), " ")
	else if (dependence=="MC1")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "log.psi1"), " ")
	else  if (dependence=="MC1R")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "log.psi1","omega"), " ")
	else if (dependence=="MC2")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "log.psi1","log.psi2"), " ")
	else if (dependence=="MC2R")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "log.psi1","log.psi2","omega" ), " ")

	covariance <- solve(Info)
	cr<- diag(1/sqrt(diag(covariance)))
	correlation <- cr %*% covariance %*% cr

	if (dependence=="ind")
	{dimnames(covariance) <- list(names.output, names.output)
	dimnames(correlation) <- list(names.output, names.output)}

	else  if (dependence=="indR")
	{dimnames(covariance) <- list(c(names.output, "omega"), c(names.output, "omega"))
	dimnames(correlation) <- list(c(names.output, "omega"), c(names.output, "omega"))}

	else if (dependence=="MC1")
	{dimnames(covariance) <- list(c(names.output, "log.psi1"), c(names.output, "log.psi1"))
	dimnames(correlation) <- list(c(names.output, "log.psi1"), c(names.output, "log.psi1"))}
	else  if (dependence=="MC1R")
	{dimnames(covariance) <- list(c(names.output, "log.psi1","omega"), c(names.output, 	"log.psi1","omega"))
	dimnames(correlation) <- list(c(names.output, "log.psi1","omega"), c(names.output, 	"log.psi1","omega"))}
	else if (dependence=="MC2")
	{dimnames(covariance) <- list(c(names.output, "log.psi1","log.psi2"), c(names.output,"log.psi1","log.psi2"))
	dimnames(correlation) <- list(c(names.output, "log.psi1","log.psi2"), c(names.output,"log.psi1","log.psi2"))}
	else if (dependence=="MC2R")
	{dimnames(covariance) <- list(c(names.output, "log.psi1","log.psi2","omega"), c(names.output, 		"log.psi1","log.psi2","omega"))
	dimnames(correlation) <- list(c(names.output, "log.psi1","log.psi2","omega"), c(names.output, 		"log.psi1","log.psi2","omega"))}

#### To compute fitted values 
#### To compute estimated transition probabilities 
	Fitted <- rep(NA, n.tot)
  if (dependence=="ind")
	{aux <- logL.bin0.aux (parameters=coefficients, X=X, data=data, trace=trace)
	Fitted[id.not.na] <- X %*% coefficients[1:(p - 1)]	 
	bi.estimate<-matrix(NaN, ncol = 1)
	prob<-aux$pij}
	else  if (dependence=="indR")
	{aux<- logL.bin0I.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	Fitted<-aux$fit
	prob<-aux$pij
	bi.estimate<- aux$bi.est
	bi.estimate<-matrix(bi.estimate, ncol = 1)
	colnames(bi.estimate)<-names.output[1]}
	else if (dependence=="MC1")
	{aux <- logL.bin1.aux (parameters=coefficients, X=X, data=data, trace=trace)
	Fitted[id.not.na] <- X %*% coefficients[1:(p - 1)]	 
	bi.estimate<-matrix(NaN, ncol = 1)
	prob<-aux$pij}
	else  if (dependence=="MC1R")
	{aux<- logL.bin1I.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	Fitted<-aux$fit
	prob<-aux$pij
	bi.estimate<- aux$bi.est
	bi.estimate<-matrix(bi.estimate, ncol = 1)
	colnames(bi.estimate)<-names.output[1]}
	else  if (dependence=="MC2")
	{aux <- logL.bin2.aux(parameters=coefficients, X=X, data=data, trace=trace)
	Fitted[id.not.na] <- X %*% coefficients[1:(p - 1)]	 
	bi.estimate<-matrix(NaN, ncol = 1)
	prob<-aux$pij}
	else  if (dependence=="MC2R")
	{aux <- logL.bin2I.aux(parameters=coefficients, X=X, data=data2, trace=trace)
	Fitted<-aux$fit
	prob<-aux$pij
	bi.estimate<- aux$bi.est
	bi.estimate<-matrix(bi.estimate, ncol = 1)
	colnames(bi.estimate)<-names.output[1]}
###
	
	ncoef<-length(coefficients)
	aic<-(2*temp$value+2*ncoef)
	ti.repl<-data[[1]]
	cumti.repl<-cumsum(ti.repl)
	n.cases<- length(ti.repl)
	y<-data2[[2]]
	counts<-data[[3]]
	residuals<-wts<-freq<-as.double(rep(0,length(y)))
	ncounts<-sum(counts[1:n.cases])

	for (j in 1: length(y))
	{
	 if (is.na(y[j])) 
	   residuals[j]<- NA
	 else
	   residuals[j]<-(y[j]-prob[j])/sqrt(prob[j]*(1-prob[j]))
	}


	Fitted <- plogis(Fitted)
	Fitted[is.na(y)] <- NA
	res<-as.double(rep(0,n.time))
	
	for (j in 1: n.time)
	{ k2<-0
	soma.n<-soma.d<-0
     	for ( i in 1:n.cases)
	{ 
	k3<-k2+j
 	if(!is.na(y[k3]))
     {
  		soma.n<-soma.n+(y[k3]-prob[k3])
  		soma.d<-soma.d+(prob[k3]*(1-prob[k3]))
      }
	k2<-cumti.repl[i]
  	}
	res[j]<-soma.n/sqrt(soma.d)
	}

	y.matrix<-matrix(y,ncol=n.time,byrow=TRUE)
	y.av<-apply(y.matrix,2,mean,na.rm=TRUE)
	Fitted.matrix<-matrix(Fitted,ncol=n.time,byrow=TRUE)
	Fitted.av<-apply(Fitted.matrix,2,mean,na.rm=TRUE)
	prob.matrix<-matrix(prob,ncol=n.time,byrow=TRUE)


bl<- new("bild", coefficients = coefficients, se = se, covariance =covariance, correlation=correlation, 
	log.likelihood=- temp$value, message = temp$convergence, n.cases=n.cases, ni.cases=ni.cases, aic=aic,    
  residuals=residuals, s.residuals=res, ind.probability=prob, prob.matrix=prob.matrix, Fitted=Fitted, 
	bi.estimate=bi.estimate, Fitted.av=Fitted.av, Time=Time, model.matrix=X, y.matrix=y.matrix, 
	subset.data=subset.data, y.av=y.av, f.value=f.value, call=call)

}


