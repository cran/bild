\name{vareff-methods}
\docType{methods}
\alias{vareff-methods}
\alias{vareff,bild-method}
\title{Methods for function \code{vareff}}
\description{Methods for function \code{vareff } extracting the variance estimates of random effects of a fitted model object from class  \code{\link[=bild-class]{bild}}. }

\usage{
\S4method{vareff}{bild}(object)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=bild-class]{bild}}.}
    }  

\section{Methods}{
\describe{
\item{\code{signature(object="bild")}:}{vareff for \code{\link{bild}} object.} 
}}

\examples{
str(airpollution)

#####  dependence="indR"
air0R <- bild(wheeze~age+smoking, data=airpollution, time="age", dependence="indR")

vareff(air0R)
}

\keyword{methods}

