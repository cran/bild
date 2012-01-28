\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,bild,missing-method}
\title{Methods for function plot in package "bild"}
\description{Six plots (selectable by which) are currently available: a plot of residuals against fitted values (which=1), 
a plot of standardized residuals against time (which=2), a plot of the autocorrelation function of the residuals (which=3),  
a plot of the partial autocorrelation function of the residuals (which=4), a plot for the fitted model (which=5) 
and a plot for the individual mean profile (which=6). By default, the first five are provided. 
}
\usage{
\S4method{plot}{bild,missing}(x,which=c(1:5),ylab=NULL,main=NULL, 
        ask=prod(par("mfcol"))<length(which)&&dev.interactive(), 
        subSET,add.unadjusted=FALSE,ident=FALSE, 
        caption=c("Residuals vs Fitted", "Residuals vs Time", "ACF residuals", 
        "PACF residuals", "Individual mean profiles"), cex.caption=1)}
\arguments{ 
  \item{x}{an object of class \code{\link{bild}}.}
  \item{which}{if a subset of the plots is required, specify a subset of the numbers 1:6.}
  \item{ylab}{label to some plots (which=5 and 6).}  
  \item{main}{title to some plots in addition to the caption (which=5 and 6).} 
  \item{ask}{logical expression; if TRUE, the user is asked before each plot.}   
  \item{subSET}{logical expression indicating elements to keep in individual mean profile plots: missing values are taken as FALSE. 
  The \code{subSET} argument is only for option which=6.}
  \item{add.unadjusted}{logical expression indicating whether or not to add the unadjusted fit for plot in which=5.}     
  \item{ident}{logical expression indicating whether or not to add the number of the subject to individual mean profile plots.
  The \code{ident} argument is only for option which=6.} 
  \item{caption}{captions to appear above the plots.}   
  \item{cex.caption}{controls the size of caption.}   
    }
    \details{   
    The option which=5 provides the parametric fitted model if the dependence structure is  
    "ind" (independence), "MC1" (first order Markov Chain) or "MC2" (second order Markov Chain). 
    When the dependence structure is 
    "MC1R" (first order Markov Chain with random intercept) or 
    "MC2R" (second order Markov Chain with random intercept) the parametric adjusted fit is provided and     
    the user can set \code{add.unadjusted}=TRUE to provide the unadjusted fitted. 
    The option which=6 is used only if the random intercept is present and provides individual mean profile.
    }
  \section{Methods}{
\describe{
 \item{\code{signature(x="ANY", y="ANY")}:}{Generic function.}
    \item{\code{signature(x="bild", y="missing")}:}{Plot diagnostics for \code{\link{bild}} object.}
}}
\examples{\donttest{ 
str(locust)

####  dependence="MC2R"
Integ <- bildIntegrate(li=-2.5,ls=2.5, lig=-2.5, lsg=2.5)   
locust2r_feed1 <- bild(move~(time+I(time^2))*sex, data=locust, 
    start=NULL, subSET=feed=="1",aggregate=sex, dependence="MC2R",  
    integrate=Integ)
    


summary(locust2r_feed1)

plot(locust2r_feed1, which=5, ylab="probability of locomoting",
     add.unadjusted=TRUE)

plot(locust2r_feed1, which=6, subSET=sex=="1", main="sex==1 & Feed=1", 
    ident=TRUE)

locust2r <- bild(move~(time+I(time^2))*feed,data=locust, start=NULL, 
    trace=TRUE, aggregate=feed, dependence="MC2R", integrate=Integ)

par(mfrow=c(2,2))
plot(locust2r, which=1)
plot(locust2r, which=2)
plot(locust2r, which=3)
plot(locust2r, which=4)
par(mfrow=c(1,1))

plot(locust2r, which=5, ylab="probability of locomoting", 
    main="Feed & Unfeed groups", add.unadjusted=TRUE)

plot(locust2r, which=6, ylab="probability of locomoting",
    main="Fed & Unfed groups", ident=TRUE)
}}
\keyword{methods}
