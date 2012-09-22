\name{anova-methods}
\docType{methods}
\alias{anova-methods}
\alias{anova,bild-method}
\title{Methods for function anova in package "bild"}
\description{ Compute an analysis deviance table for two fitted model objects. }

\usage{
\S4method{anova}{bild}(object, ..., test = TRUE, correct = FALSE)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=bild-class]{bild}}.}
  \item{...}{an object of class \code{\link[=bild-class]{bild}}.}
  \item{test}{an optional logical value controlling whether likelihood ratio tests 
  should be used to compare the fitted models represented by \code{object} and by \code{y}. The default is TRUE. }  
  \item{correct}{an optional logical value controlling whether the p-value of the likelihood ratio test 
  must be corrected. The default is FALSE.}  
   }  
\details{ 
\code{correct} = TRUE is used to test the presence of a random intercept term and the solution proposed by Self and Liang (1987) is adopted 
only to the p-value. }


\section{Warning}{ 
The comparison between two models by anova will only be valid if they are fitted to the same dataset.}

\section{Methods}{
\describe{
\item{\code{signature(object = "ANY")}:}{Generic function.}
\item{\code{signature(object="bild")}:}{Anova for \code{\link{bild}} object.} 
}}

\references{Self, Steven G. and Liang, Kung-Yee (1987). Asymptotic properties of maximum likelihood estimators and likelihood 
ratio tests under nonstandard conditions. \emph{Journal of the American Statistical Association}, 82, 605-610. }


\examples{\donttest{

#####  data= locust

locust1 <- bild(move~(time+I(time^2))*feed*sex, data=locust, aggregate=feed, start=NULL, dependence="MC1")

locust2 <- bild(move~(time+I(time^2))*feed, data=locust, aggregate=feed, start=NULL, dependence="MC1")

anova(locust1,locust2)

locust3 <- bild(move~(time+I(time^2))*feed, data=locust, aggregate=feed, start=NULL, dependence="MC2")

anova(locust3,locust2)

#####  data= muscatine

# we decompose the time effect in orthogonal components
muscatine$time1 <- c(-1, 0, 1)
muscatine$time2 <- c(1, -2, 1)

musc1 <- bild(obese ~ time1 , data = muscatine, time = "time1", dependence = "MC1")

musc1r <- bild(obese ~ time1 , data = muscatine, time = "time1", dependence = "MC1R")

anova(musc1,musc1r,correct=TRUE)
}}

\keyword{methods}

