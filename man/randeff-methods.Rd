\name{randeff-methods}
\docType{methods}
\alias{randeff-methods}
\alias{randeff,bild-method}
\title{Methods for function \code{randeff}}
\description{Methods for function \code{randeff} extracting conditional random effects of a fitted model object from class  \code{\link[=bild-class]{bild}}. }

\usage{
\S4method{randeff}{bild}(object)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=bild-class]{bild}}.}
    }  

\section{Methods}{
\describe{
\item{\code{signature(object="bild")}:}{fitted for \code{\link{bild}} object.} 
}}

\examples{
str(airpollution)

#####  dependence="indR"
air0R <- bild(wheeze~age+smoking, data=airpollution, time="age", dependence="indR")

randeff(air0R)
}

\keyword{methods}

